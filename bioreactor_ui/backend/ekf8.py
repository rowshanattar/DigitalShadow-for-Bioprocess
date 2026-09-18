"""
ekf8.py — Design-A 8-state "healthy shadow" EKF.

State: [X, G, E, µmax_G, µmax_E, T, pH, DO]   — T, pH and DO are ESTIMATED
states here, integrated from the HEALTHY ODEs. They are never fed in. That is
the whole difference from Design B (`ekf.py`), where T/pH/DO are measured inputs
folded into φ and the filter therefore *follows* every actuator fault.

Because the filter only ever integrates healthy dynamics and is corrected by the
ethanol sensor alone, it answers "how should this reactor be behaving?" The
innovation ν = z_E − Ê is then the fault signal, and the CUSUM in `cusum.py`
runs on ν/√S.

Ported from fault_detection_CU.ipynb:
  * matrices P0 / Q / R / H         — Cell 4
  * forward-difference Jacobian     — Cell 12
  * correction + physical clamps    — Cell 18 (`run_ekf`)

State PREDICTION does not live here: the simulator advances this filter through
the very same integrator that produces the healthy counterfactual trace
(`Simulation._advance` with an empty fault list), so "the EKF predicts with
healthy dynamics" is true by construction instead of by a second copy of the
ODEs that could drift out of sync. This class owns the covariance and the
correction.
"""
from __future__ import annotations

import numpy as np

N_STATE = 8


P0 = np.diag([0.1, 0.02, 0.02, 1e-5, 1e-5, 1e-4, 1e-6, 1e-4])
Q = np.diag([1e-3, 1e-3, 1e-3, 1e-9, 1e-9, 1e-4, 1e-6, 1e-4])
H = np.zeros((1, N_STATE))
H[0, 2] = 1.0                      # ethanol-only observation


def compute_jacobian(rhs, x, eps=1e-6):
    """Forward-difference Jacobian (8x8) of the HEALTHY dynamics at x."""
    f0 = np.array(rhs(0.0, x), dtype=float)
    J = np.zeros((N_STATE, N_STATE))
    for j in range(N_STATE):
        xp = np.array(x, dtype=float)
        xp[j] += eps
        J[:, j] = (np.array(rhs(0.0, xp), dtype=float) - f0) / eps
    return J


class EKF8:
    """Continuous-discrete healthy-shadow EKF, stepped one frame at a time."""

    def __init__(self, params: dict):
        ini, th, ph, ox = (params["initial"], params["thermal"],
                           params["ph"], params["oxygen"])
        self.x = np.array([ini["X0"], ini["G0"], ini["E0"],
                           ini["muG0"], ini["muE0"],
                           th["T_set"], ph["pH_set"], ox["DO0"]], dtype=float)
        self.P = P0.copy()
        self.Q = Q.copy()
        sigma_E = params["ekf"]["sigma_E"]
        self.R = np.array([[max(sigma_E ** 2, 1e-4)]])
        self.I8 = np.eye(N_STATE)

        # The shadow runs its own pH controller: it is a *healthy* reactor, so
        # its integral term must not see the faulted run's error history.
        self.pi_int = 0.0
        # last corrected state — needed for the mass-balance clamp below
        self.prev_x = self.x.copy()

    # ── covariance propagation over one frame ────────────────────────────
    def predict_cov(self, dt_h: float, rhs):
        """P <- Phi P Phi' + Q dt, with Phi = I + F dt at the current state.

        `rhs` is the healthy plant RHS closure for this frame. A frame is one
        minute (dt = 0.017 h) against dynamics of O(0.2 1/h), so the first-order
        Phi is accurate here and costs one Jacobian instead of a second ODE solve.
        """
        F = compute_jacobian(rhs, self.x)
        Phi = self.I8 + F * dt_h
        self.P = Phi @ self.P @ Phi.T + self.Q * dt_h

    # ── ethanol measurement update ───────────────────────────────────────
    def correct(self, z_ethanol: float):
        """Joseph-form update on ethanol. Returns (innovation, S)."""
        S = H @ self.P @ H.T + self.R
        K = self.P @ H.T @ np.linalg.inv(S)
        innovation = float(z_ethanol - self.x[2])
        self.x = self.x + (K @ np.array([[innovation]])).flatten()

        # keep concentrations and DO physical
        for idx in (0, 1, 2, 7):
            self.x[idx] = max(self.x[idx], 0.0)
        self.x[6] = max(self.x[6], 1.0)

        # Batch mass balance, enforced on the CORRECTION. The plant has no death
        # term (dX = (muG+muE)*X >= 0) and no feed (dG = -(muG/Ygx)*X <= 0), so a
        # correction that reverses either is inadmissible, not merely imprecise.
        # It matters because µmax is pinned: the filter cannot absorb a fault by
        # retuning growth, so without this it absorbs it into X and G instead —
        # dragging X below the faulted truth and CREATING glucose.
        self.x[0] = max(self.x[0], self.prev_x[0])   # biomass cannot fall
        self.x[1] = min(self.x[1], self.prev_x[1])   # glucose cannot rise

        IKH = self.I8 - K @ H
        self.P = IKH @ self.P @ IKH.T + K @ self.R @ K.T
        self.prev_x = self.x.copy()
        return innovation, float(S[0, 0])

    def sigma(self):
        """Standard deviations of [X, G, E] for the optional ±σ bands."""
        d = np.clip(np.diag(self.P)[:3], 0.0, None)
        return np.sqrt(d)
