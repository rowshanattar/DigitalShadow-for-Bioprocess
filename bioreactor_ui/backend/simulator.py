"""
simulator.py — stepwise, streamable integrator.

Advances the 7-state plant one frame (default 1 sim-minute) at a time so the
frontend can draw growing curves. Each frame also doubles as the discrete pH
control tick (pumps are sampled at the frame start and held fixed), which is
exactly the Cell-11 mechanism. The EKF is stepped alongside and corrected on
noisy ethanol samples at its cadence.

Fault injection is a plain method: faults can be scheduled before the run or
appended live, always with a guard that t_fault has not already passed.

Two estimator architectures share this loop:
  * Design A (default) — the 8-state healthy-shadow EKF8. Its state is advanced
    through `_advance` with an EMPTY fault list, i.e. the exact code path that
    produces the healthy counterfactual, then corrected on ethanol. A CUSUM runs
    on the normalised innovation.
  * Design B (`ekf.old_design`) — the legacy 5-state EKF, fed the measured
    T/pH/DO through φ.
"""
from __future__ import annotations

import numpy as np

from .core import absolute_modifier
from .plant import resolve_params, build_modifiers, make_rhs, instantaneous_rates, thermal_steady_state
from .faults import apply_faults, normalize_fault
from .ekf import EKF
from .ekf8 import EKF8
from .cusum import cusum_step


class Simulation:
    def __init__(self, cfg: dict):
        self.params = resolve_params(cfg)
        self.f_T, self.f_pH, self.f_DO = build_modifiers(self.params)

        ini, th, ph, ox = (self.params["initial"], self.params["thermal"],
                           self.params["ph"], self.params["oxygen"])
        self.y = np.array([ini["X0"], ini["G0"], ini["E0"], ini["muG0"],
                           ini["muE0"], th["T_set"], ph["pH_set"], ox["DO0"]],
                          dtype=float)
        self.t = 0.0
        self.t_final = self.params["sim"]["t_final"]
        self.dt_h = self.params["sim"]["dt_min"] / 60.0

        self.active_faults: list[dict] = []

        # PI pH-controller integral state — one per trajectory (the healthy
        # counterfactual is controlled independently of the faulted run).
        self.pi_int = 0.0
        self.pi_int_healthy = 0.0

        # EKF — Design A (8-state healthy shadow) unless old_design is set, in
        # which case the legacy 5-state Design-B filter runs instead.
        ek = self.params["ekf"]
        self.ekf_on = bool(ek["on"])
        self.old_design = bool(ek.get("old_design", False))
        if not self.ekf_on:
            self.ekf = None
        else:
            self.ekf = EKF(self.params) if self.old_design else EKF8(self.params)
        # the REAL sensor noise — the filter's own assumed sigma_E lives in its R
        self.sigma_E = ek["sigma_E_sensor"]
        self.cadence_h = max(ek["cadence_min"], 1e-6) / 60.0
        self.next_meas_t = self.cadence_h
        self.rng = np.random.default_rng(0)

        # ── CUSUM detector state (Design A only): two floats + a latch ──────
        self.cusum_k, self.cusum_h = ek["cusum_k"], ek["cusum_h"]
        self.sp = self.sn = 0.0
        self.alarm_t: float | None = None
        self.last_nu = self.last_z = 0.0

        # running biomass NRMSE accumulators (EKF estimate vs. true X)
        self._se_X, self._n_X, self._Xmax = 0.0, 0, 0.0
        # Peak ethanol the HEALTHY counterfactual has reached. If that stays
        # under a few sigma the detector is structurally blind: there is no
        # signal to lose, so 'quiet' means 'cannot see', not 'nothing wrong'.
        self._hE_peak = 0.0
        # The specific reason the ethanol channel may be deaf, computed once.
        # A proportional heater needs a standing error to make power at all, so
        # the vessel sits ~2 degC BELOW setpoint all batch; cardinals chosen
        # against T_set rather than against THIS number silently cripple the
        # healthy reference, and a crippled culture respires everything it takes
        # up and ferments nothing.
        self.T_ss = thermal_steady_state(self.params)
        self.fT_ss = float(self.f_T(self.T_ss))
        # ...and the pH twin of it. Reported alongside fT_ss because the BLIND
        # badge used to name the T window as the cause unconditionally, which
        # points the user at the wrong knob whenever pH is the one that is off.
        # This is the ABSOLUTE modifier (relative to the organism's optimum), not
        # f_pH(pH_set) — that is 1 by construction, since f_pH is anchored there.
        self.fpH_set = float(absolute_modifier(
            self.params["ph"]["pH_set"], self.params["ph"]["pH_min"],
            self.params["ph"]["pH_opt"], self.params["ph"]["pH_max"], n=1))

        self.rows: list[dict] = []       # true trajectory — CSV / PNG
        self.ekf_rows: list[dict] = []    # EKF estimates — PNG overlay
        self.det_rows: list[dict] = []    # innovation + CUSUM arms (Design A)
        # One entry per ethanol SAMPLE (not per frame): the sensor value the
        # rig actually returned, plus the shadow's response to it. The MMAE
        # bank replays this exact sequence into every hypothesis filter, so it
        # has to be kept — nothing can be re-drawn from the RNG after the fact.
        self.meas_rows: list[dict] = []
        self.alarm_idx: int | None = None    # index into meas_rows
        self.mmae: dict | None = None        # last root-cause result
        # Counterfactual "no-fault" baseline, integrated in lockstep so the
        # user can compare the faulted run against what would have happened.
        self.healthy_y = self.y.copy()
        self.healthy_rows: list[dict] = []

    # ── fault scheduling ─────────────────────────────────────────────────
    def inject_fault(self, fault: dict):
        f = normalize_fault(fault)
        t_fault = float(f.get("t_fault", self.t))
        if t_fault < self.t - 1e-9:
            return False, (f"t_fault = {t_fault:.2f} h has already passed "
                           f"(sim is at {self.t:.2f} h).")
        f["t_fault"] = t_fault
        self.active_faults.append(f)
        return True, f"{f['fault']} scheduled at t = {t_fault:.2f} h."

    @property
    def done(self) -> bool:
        return self.t >= self.t_final - 1e-9

    # ── RK4 step of the plant over dt with fixed actuator context ────────
    def _rk4(self, y, t, dt, F_A, F_B, heater_frac, UA_eff, kLa_eff):
        rhs = make_rhs(self.params, self.f_T, self.f_pH, self.f_DO,
                       F_A, F_B, heater_frac, UA_eff, kLa_eff)
        k1 = np.array(rhs(t, y))
        k2 = np.array(rhs(t + dt / 2, y + dt / 2 * k1))
        k3 = np.array(rhs(t + dt / 2, y + dt / 2 * k2))
        k4 = np.array(rhs(t + dt, y + dt * k3))
        y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        y[0] = max(y[0], 0.0); y[1] = max(y[1], 0.0); y[2] = max(y[2], 0.0)
        y[7] = max(y[7], 0.0)   # DO can't go negative
        return y

    # ── integrate one frame with PI pH control on a short tick ────────────
    def _advance(self, y, t0, t1, faults, integ):
        ph = self.params["ph"]
        pH_set = ph["pH_set"]
        Kc, Ki = ph["Kc_pump"], ph["Ki_pump"]
        F_A_max, F_B_max = ph["F_A"], ph["F_B"]
        dt_ctrl = max(ph["DT_CTRL_s"], 0.5) / 3600.0
        y = np.array(y, dtype=float)
        t = t0
        while t < t1 - 1e-12:
            h = min(dt_ctrl, t1 - t)
            pH = y[6]
            # PI (throttled) control: signed command = proportional + integral.
            # The integral term is what "slowly pushes pH back up" to the
            # setpoint after a fault dips it — a pure-P law would leave an offset.
            err = pH_set - pH                       # >0 → too acidic, add base
            u = Kc * err + Ki * integ               # signed pump command [mL/h]
            if u >= 0.0:
                base_F_B, base_F_A = min(u, F_B_max), 0.0    # base raises pH
            else:
                base_F_A, base_F_B = min(-u, F_A_max), 0.0   # acid lowers pH
            # Anti-windup: don't accumulate error while the pump is saturated
            # and the error would only drive it further into saturation.
            saturated = u > F_B_max or u < -F_A_max
            if not (saturated and (err > 0.0) == (u > 0.0)):
                integ += err * h
            F_A, F_B, heater_frac, UA_eff, kLa_eff = apply_faults(
                faults, t, self.params, base_F_A, base_F_B)
            y = self._rk4(y, t, h, F_A, F_B, heater_frac, UA_eff, kLa_eff)
            t += h
        return y, integ

    # ── RHS closure at a given filter state (Design A / MMAE Jacobian) ────
    def _rhs_for(self, y, integ, faults, t):
        """The plant RHS in the actuator context state `y` implies at time `t`.

        Same pump law as `_advance` and the same `apply_faults` resolution, so
        the Jacobian is taken of exactly the dynamics the filter is being
        propagated with. `faults=[]` gives the healthy shadow; an MMAE
        hypothesis passes its own fault list here (`mmae.py`).
        """
        ph = self.params["ph"]
        u = ph["Kc_pump"] * (ph["pH_set"] - y[6]) + ph["Ki_pump"] * integ
        if u >= 0.0:
            F_B, F_A = min(u, ph["F_B"]), 0.0
        else:
            F_A, F_B = min(-u, ph["F_A"]), 0.0
        F_A, F_B, heater_frac, UA_eff, kLa_eff = apply_faults(
            faults, t, self.params, F_A, F_B)
        return make_rhs(self.params, self.f_T, self.f_pH, self.f_DO,
                        F_A, F_B, heater_frac, UA_eff, kLa_eff)

    def _healthy_rhs(self, y, integ):
        """The fault-free RHS — `_rhs_for` with no faults (the healthy shadow)."""
        return self._rhs_for(y, integ, [], 0.0)

    # ── one frame ────────────────────────────────────────────────────────
    def step(self) -> dict:
        p = self.params
        t0 = self.t
        t_end = min(self.t + self.dt_h, self.t_final)

        # Faulted (real) trajectory and the no-fault counterfactual, in lockstep
        self.y, self.pi_int = self._advance(
            self.y, t0, t_end, self.active_faults, self.pi_int)
        self.healthy_y, self.pi_int_healthy = self._advance(
            self.healthy_y, t0, t_end, [], self.pi_int_healthy)
        self.t = t_end

        X, G, E, m1, m2, T, pHv, DOv = self.y
        hX, hG, hE, _, _, hT, hpH, hDO = self.healthy_y
        muG, muE = instantaneous_rates(p, self.f_T, self.f_pH, self.f_DO, self.y)

        ekf_frame = cusum_frame = None
        if self.ekf_on:
            if self.old_design:
                # Design B: measured T, pH, DO enter as the modifier φ, so the
                # filter compensates for them and follows the fault.
                # Two numbers, not one: phi_g scales the rate, f_do caps the
                # respiratory branch (Crabtree). See core.carbon_fluxes.
                self.ekf.predict(self.dt_h,
                                 self.f_T(T) * self.f_pH(pHv), self.f_DO(DOv))
            else:
                # Design A: predict with the HEALTHY dynamics — the same call
                # that produces the no-fault counterfactual above. Nothing about
                # the faulted run enters here; the ethanol correction below is
                # the filter's only contact with reality.
                self.ekf.x, self.ekf.pi_int = self._advance(
                    self.ekf.x, t0, t_end, [], self.ekf.pi_int)
                self.ekf.predict_cov(
                    self.dt_h, self._healthy_rhs(self.ekf.x, self.ekf.pi_int))

            if self.t >= self.next_meas_t - 1e-9:
                z = E + self.rng.normal(0.0, self.sigma_E)
                nu, S = self.ekf.correct(z)
                self.next_meas_t += self.cadence_h
                if not self.old_design:
                    # DETECT, right here in the correction step. z_norm is the
                    # innovation normalised by the filter's own predicted spread
                    # sqrt(S), so a fixed threshold means the same thing at every
                    # step. The arms are never reset — the alarm latches.
                    z_norm = nu / max(S ** 0.5, 1e-12)
                    self.sp, self.sn, fired = cusum_step(
                        self.sp, self.sn, z_norm, self.cusum_k, self.cusum_h)
                    self.last_nu, self.last_z = nu, z_norm
                    self.meas_rows.append({
                        "time_h": self.t, "frame": len(self.rows), "z": float(z),
                        "nu": nu, "S": S, "Sp": self.sp, "Sn": self.sn,
                    })
                    if fired and self.alarm_t is None:
                        self.alarm_t = self.t
                        self.alarm_idx = len(self.meas_rows) - 1

            sx = self.ekf.sigma()
            ex = self.ekf.x
            ekf_frame = {
                "X": float(ex[0]), "G": float(ex[1]), "E": float(ex[2]),
                "muG": float(ex[3]), "muE": float(ex[4]),
                "sX": float(sx[0]), "sG": float(sx[1]), "sE": float(sx[2]),
            }
            if not self.old_design:
                # The shadow's T/pH/DO are estimated states, so they get their
                # own overlay. They sit almost on the healthy curve by design:
                # ethanol cannot observe them, and they move only through the
                # corrected X/G/E feeding back into metabolic heat and
                # acidogenesis. That is correct behaviour, not a stuck trace.
                ekf_frame.update({"T": float(ex[5]), "pH": float(ex[6]),
                                  "DO": float(ex[7])})
                cusum_frame = {"nu": self.last_nu, "z": self.last_z,
                               "blind": bool(self._hE_peak < 3.0 * self.sigma_E),
                               "hE_peak": float(self._hE_peak),
                               "T_ss": float(self.T_ss), "fT_ss": self.fT_ss,
                               "fpH_set": self.fpH_set,
                               "Sp": self.sp, "Sn": self.sn,
                               "h": self.cusum_h, "alarm_t": self.alarm_t}
                self.det_rows.append({"time_h": self.t, **cusum_frame})
            self.ekf_rows.append({"time_h": self.t, **ekf_frame})

        self.rows.append({
            "time_h": self.t, "T": float(T), "pH": float(pHv),
            "DO": float(DOv), "G": float(G), "E": float(E), "X": float(X),
            "mu_G": float(muG), "mu_E": float(muE),
        })
        self._hE_peak = max(self._hE_peak, float(hE))
        self.healthy_rows.append({
            "time_h": self.t, "X": float(hX), "G": float(hG), "E": float(hE),
            "T": float(hT), "pH": float(hpH), "DO": float(hDO),
        })

        # Running estimator error on BIOMASS — the hidden state the filter exists
        # to estimate. (Ethanol would only measure how well it fits its own
        # input.) Normalised by peak true X so it reads as a % and stays
        # comparable across runs; this is the number to watch while sweeping
        # ekf.sigma_E_sensor.
        if ekf_frame is not None:
            self._se_X += (ekf_frame["X"] - float(X)) ** 2
            self._n_X += 1
            self._Xmax = max(self._Xmax, float(X))
            nrmse_X = 100.0 * (self._se_X / self._n_X) ** 0.5 / max(self._Xmax, 1e-9)
        else:
            nrmse_X = None

        return {
            "t_h": self.t,
            "true": {"X": float(X), "G": float(G), "E": float(E),
                     "T": float(T), "pH": float(pHv), "DO": float(DOv),
                     "muG": float(muG), "muE": float(muE)},
            "healthy": {"X": float(hX), "G": float(hG), "E": float(hE),
                        "T": float(hT), "pH": float(hpH), "DO": float(hDO)},
            "ekf": ekf_frame,
            "cusum": cusum_frame,      # None in Design B — no detector there
            "nrmse_X": nrmse_X,
            "faults": [{"fault": f["fault"], "t_fault": f["t_fault"]}
                       for f in self.active_faults],
        }
