"""
plant.py — 7-state ground-truth plant.

State: [X, G, E, µmax_G, µmax_E, T, pH]

Ported from fault-injection.ipynb Cell 9 (pH-extended plant) but fully
parameterised. The pump flows (F_A, F_B), heater state and effective UA are
passed in per control tick so the simulator can (a) reproduce the discrete
sawtooth pH controller and (b) inject actuator faults by overriding them.
"""
from __future__ import annotations

import copy
import math
import warnings
from .core import (DEFAULTS, absolute_modifier, carbon_fluxes, growth_rates,
                   make_cardinal_modifier, make_do_modifier)


# ── Van't Riet kLa correlation constants (non-coalescing broth, Rushton) ──
# Fixed vessel geometry / power characteristics; the tunable operating inputs
# (N_rpm, Q_air_Lmin) live in DEFAULTS["oxygen"] so the UI and the aeration
# fault can vary them. kLa = VR_C · (Pg/V)^VR_a · v_s^VR_b   [1/s → ×3600 → 1/h]
_D_VESSEL_M = 0.090   # vessel inner diameter [m], DN 90  — *confirm from manual*
_D_IMP_FRAC = 1.0 / 3.0   # impeller/vessel diameter ratio (Rushton)
_N_POWER    = 5.0     # power number per 6-blade Rushton turbine, turbulent
_N_IMPELLERS = 2      # two 6-blade Rushton turbines on the shaft (Minifors manual)
_RHO_BROTH  = 1000.0  # broth density [kg/m³] ~ water
_PG_P0      = 0.5     # gassed/ungassed power ratio — *confirm*
_VR_C, _VR_A, _VR_B = 0.002, 0.7, 0.2


def compute_kLa(N_rpm: float, Q_air_Lmin: float, V_L: float) -> float:
    """Volumetric O₂ transfer coefficient kLa [1/h] from agitation + aeration.

    Van't Riet correlation. Reducing N_rpm or Q_air_Lmin (agitation/aeration
    fault) lowers kLa, which starves the DO supply term.
    """
    V_m3 = V_L / 1000.0
    A_cross = math.pi / 4.0 * _D_VESSEL_M ** 2
    D_imp = _D_VESSEL_M * _D_IMP_FRAC
    N_rps = N_rpm / 60.0
    v_s = (Q_air_Lmin / 1000.0 / 60.0) / A_cross          # superficial gas vel [m/s]
    # ungassed power [W]; two Rushtons draw ≈ the sum of their individual powers
    P0 = _N_IMPELLERS * _N_POWER * _RHO_BROTH * N_rps ** 3 * D_imp ** 5
    Pg_V = (_PG_P0 * P0) / V_m3                            # gassed power/vol [W/m³]
    if Pg_V <= 0.0 or v_s <= 0.0:
        return 0.0
    return _VR_C * (Pg_V ** _VR_A) * (v_s ** _VR_B) * 3600.0   # [1/h]


def resolve_params(user_cfg: dict | None) -> dict:
    """Deep-merge a (partial) user config over DEFAULTS."""
    params = copy.deepcopy(DEFAULTS)
    if user_cfg:
        for group, vals in user_cfg.items():
            if group in params and isinstance(vals, dict):
                params[group].update(vals)
            else:
                params[group] = vals
    return params


def _assert_monotone(f, lo, hi, name, n=300):
    """Every modifier must rise monotonically towards its optimum.

    A Rosso denominator that changes sign inside the window is invisible in the
    output — the clip presents the pole as a clean 0/1 step — so nothing but a
    shape check catches it. This guard is the regression test for the f_pH step
    function that made a stuck acid pump indistinguishable from an oxygen crash
    (see the note above `make_cardinal_modifier` in core.py). It previously did
    not exist, which is exactly why that sat undetected.
    """
    step = (hi - lo) / (n - 1)
    prev = f(lo)
    for i in range(1, n):
        v = lo + i * step
        cur = f(v)
        if cur - prev < -1e-9:
            raise ValueError(
                f"{name} is not monotone rising on [{lo}, {hi}]: falls at "
                f"{v:.3f}. A Rosso denominator has changed sign and the clip is "
                f"hiding the pole as a step function.")
        prev = cur


def thermal_steady_state(params: dict) -> float:
    """The temperature the healthy P-controller actually settles at [degC].

    This is NOT T_set. A proportional heater needs a standing error to produce
    any power at all, so at equilibrium (ignoring metabolic heat)

        K_p (T_set - T) = UA (T - T_amb)   ->   T = (K_p T_set + UA T_amb)/(K_p + UA)

    With the defaults that is (50*30 + 18*22)/68 = 27.88 degC, i.e. the vessel
    runs 2.1 degC BELOW setpoint all batch. Cardinal temperatures have to be
    chosen against THIS number, not against T_set: a T_min of 27 looks like a
    3 degC margin and is really a 0.9 degC margin.
    """
    th = params["thermal"]
    return ((th["K_p"] * th["T_set"] + th["UA"] * th["T_amb"])
            / (th["K_p"] + th["UA"]))


def build_modifiers(params: dict):
    """Return (f_T, f_pH, f_DO) callables from the current cardinal parameters."""
    th, ph, ox = params["thermal"], params["ph"], params["oxygen"]

    n_T = 2 if (th["T_opt"] - th["T_min"]) >= (th["T_max"] - th["T_opt"]) else 1
    f_T = make_cardinal_modifier(th["T_min"], th["T_opt"], th["T_max"],
                                 anchor=th["T_set"], n=n_T)
    # pH arms are asymmetric — n = 1, the pole-free form. Anchored at the
    # SETPOINT, not the optimum, to match the notebook: when the controller
    # holds the reactor closer to the biological optimum than the anchor is,
    # f_pH > 1 is the physically correct answer and must not be clipped away.
    f_pH = make_cardinal_modifier(ph["pH_min"], ph["pH_opt"], ph["pH_max"],
                                  anchor=ph["pH_set"], n=1)
    f_DO = make_do_modifier(ox["K_O"], ox["DO_star"])

    _assert_monotone(f_T, th["T_min"] + 1e-6, th["T_opt"], "f_T")
    _assert_monotone(f_pH, ph["pH_min"] + 1e-6, ph["pH_opt"], "f_pH")
    _assert_monotone(f_DO, 0.0, ox["DO_star"], "f_DO")

  
    # HOW GOOD IS THE SETPOINT ITSELF? Both modifiers are re-anchored AT their
    # setpoint, so f_T(T_set) and f_pH(pH_set) are identically 1 and can never
    # answer that — they measure distance from the setpoint, not the quality of
    # it. `absolute_modifier` is the un-anchored value relative to the organism's
    # true optimum, which is what says whether the HEALTHY reference is any good.
 
    abs_pH = absolute_modifier(ph["pH_set"], ph["pH_min"], ph["pH_opt"],
                               ph["pH_max"], n=1)
    if abs_pH < 0.5:
        warnings.warn(
            f"cardinal pH {(ph['pH_min'], ph['pH_opt'], ph['pH_max'])} leaves the "
            f"culture at {100 * abs_pH:.0f} % of its peak rate at the pH the "
            f"controller holds (pH_set = {ph['pH_set']}, pH_opt = {ph['pH_opt']}). "
            f"The HEALTHY reference is itself crippled, so a fault that moves pH "
            f"back towards the optimum will make the culture grow BETTER than "
            f"'healthy' and the detector will read it backwards. Move pH_set "
            f"towards pH_opt, or widen the cardinal window.",
            RuntimeWarning, stacklevel=2)

    T_ss = thermal_steady_state(params)
    if f_T(T_ss) < 0.5:
        warnings.warn(
            f"cardinal T {(th['T_min'], th['T_opt'], th['T_max'])} leaves "
            f"f_T = {f_T(T_ss):.2f} at the steady state the P-controller "
            f"actually reaches ({T_ss:.2f} degC, not T_set = {th['T_set']}). "
            f"The HEALTHY culture is already growing at {100 * f_T(T_ss):.0f} % "
            f"of its rate, so it barely ferments and the ethanol sensor has "
            f"almost nothing to see, so faults that only SLOW growth (T, pH) "
            f"cannot raise the alarm at any threshold. Oxygen faults still can, "
            f"because they redirect carbon into fermentation rather than "
            f"slowing it. Widen the cardinal window until f_T at {T_ss:.1f} degC "
            f"is comfortably above 0.5 - and note BOTH arms matter: "
            f"(27, 30, 32) gives f_T = 0.28, (27, 30, 33) gives 0.50.",
            RuntimeWarning, stacklevel=2)
    return f_T, f_pH, f_DO


def make_rhs(params: dict, f_T, f_pH, f_DO, F_A: float, F_B: float,
             heater_frac: float, UA_eff: float, kLa_eff: float):
    """Build the 8-state RHS closure for one control tick.

    F_A, F_B are held fixed over the tick (discrete pH control). heater_frac
    (1.0 = healthy, 0.0 = dead, in between = derated) and
    UA_eff encode heater-failure and vessel-leak faults; kLa_eff encodes the
    agitation/aeration fault (reduced O₂ transfer).
    """
    kin, th, ph, ox = (params["kinetics"], params["thermal"],
                       params["ph"], params["oxygen"])
    Yex = kin["Yex"]
    V, rhoCp, Y_QX = th["V"], th["rhoCp"], th["Y_QX"]
    T_set, T_amb, K_p, Q_max = th["T_set"], th["T_amb"], th["K_p"], th["Q_max"]
    alpha_meta, alpha_pump = ph["alpha_meta"], ph["alpha_pump"]
    DO_star, DO_sat_mgL = ox["DO_star"], ox["DO_sat_mgL"]

    def rhs(t, y):
        X, G, E, m1, m2, T, pH, DO = y
        X = max(X, 0.0); G = max(G, 0.0); E = max(E, 0.0)
        pH = max(pH, 1.0); DO = max(DO, 0.0)

        # phi_g scales the RATE (T and pH only). Oxygen sets carbon FATE, so a
        # DO fault pushes ethanol UP while a T or pH fault pushes it DOWN.
        phi_g = f_T(T) * f_pH(pH)
        f_do = f_DO(DO)
        q_ox, q_fm, muE, mu = carbon_fluxes(G, E, m1, m2, kin, phi_g, f_do)
        q_Eox = muE / Yex                    # ethanol consumed [g E / g X / h]

        dX = mu * X
        dG = -(q_ox + q_fm) * X
        dE = (q_fm * kin["YGE_FERM"] - q_Eox) * X
        dm1 = 0.0
        dm2 = 0.0

        # Thermal ODE (W → J/h via ×3600; metabolic heat already J/h)
        Q_heat = heater_frac * max(0.0, min(Q_max, K_p * (T_set - T)))
        Q_loss = UA_eff * (T - T_amb)
        Q_metab = Y_QX * mu * X * V
        dT = (3600.0 * (Q_heat - Q_loss) + Q_metab) / (V * rhoCp)

        # pH ODE: yeast acidogenesis vs. (fixed) pump flows this tick
        yeast_acid = alpha_meta * mu * X
        dpH = -yeast_acid + alpha_pump * F_B - alpha_pump * F_A
        if pH <= 1.0 and dpH < 0:      # physical buffering floor
            dpH = 0.0

        # DO ODE (%-sat/h): OTR supply − respiratory OUR demand.
        # Demand tracks the OXIDATIVE flux only, so oxygen is not consumed by
        # the fermentative branch. It previously tracked total growth, which
        # double-counted: fermentation would have drawn oxygen it never uses.
        supply = kLa_eff * (DO_star - DO)
        demand = kin["K_O2"] * (q_ox + kin["YO_E_REL"] * q_Eox) * X
        dDO = supply - demand / DO_sat_mgL * 100.0
        if DO <= 0.0 and dDO < 0:       # physical floor: DO can't go negative
            dDO = 0.0

        return [dX, dG, dE, dm1, dm2, dT, dpH, dDO]

    return rhs


def instantaneous_rates(params: dict, f_T, f_pH, f_DO, y):
    """µ_G, µ_E at a state y (for CSV / display columns).

    µ_G is now the growth carried by glucose (respired + fermented) and µ_E the
    growth carried by ethanol — the same two display channels as before, but
    read off the carbon split rather than off two Monod terms.
    """
    X, G, E, m1, m2, T, pH, DO = y
    phi_g = f_T(T) * f_pH(max(pH, 1.0))
    f_do = f_DO(max(DO, 0.0))
    return growth_rates(max(X, 0.0), max(G, 0.0), max(E, 0.0), m1, m2,
                        params["kinetics"], phi_g, f_do)
