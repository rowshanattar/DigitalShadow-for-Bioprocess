"""
ekf.py — Design-B 5-state Extended Kalman Filter.

State: [X, G, E, µmax_G, µmax_E]   (T, pH, DO are *measured inputs*, not states)

Ported from ekf_bioreactor.ipynb:
  * combined state+covariance ODE — Cell 13
  * predict / Joseph-form correct — Cell 15

The environmental modifiers are passed into each predict step, so the filter
compensates for known T/pH/DO changes — the mechanism by which it stays blind to
actuator faults on measured inputs (Design B).

TWO CHANGES from the original, both to follow the Crabtree plant:

  * The measured environment enters as TWO numbers, not one. `phi_g = f_T*f_pH`
    scales the glucose uptake rate; `f_do` caps the respiratory branch. Folding
    f_DO into a single phi (the previous behaviour) made a DO fault reduce
    ethanol, which is the wrong SIGN — see core.carbon_fluxes.
  * The Jacobian is now computed by FORWARD DIFFERENCES rather than by sympy.
    The carbon split contains `min(q_G, C_RESP*f_do)`, so the RHS is not
    differentiable at the respiro-fermentative switch — exactly where a DO fault
    lives — and a symbolic Jacobian would be wrong precisely there. This mirrors
    ekf8.compute_jacobian and the notebook's A.6.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .core import carbon_fluxes

N_STATE = 5
H = np.array([[0.0, 0.0, 1.0, 0.0, 0.0]])   # ethanol-only observation


class EKF:
    """Continuous-discrete EKF stepped one frame at a time."""

    def __init__(self, params: dict):
        self.kin = params["kinetics"]
        ini = params["initial"]
        self.Yex = self.kin["Yex"]

        self.x = np.array([ini["X0"], ini["G0"], ini["E0"],
                           ini["muG0"], ini["muE0"]], dtype=float)
        self.P = np.diag([0.1, 0.02, 0.02, 0.2, 0.02])
        self.Q = np.diag([1e-3, 1e-3, 1e-3, 1e-3, 1e-3])
        sigma_E = params["ekf"]["sigma_E"]
        self.R = np.array([[max(sigma_E ** 2, 1e-4)]])
        self.I5 = np.eye(N_STATE)
        self._phi_g, self._f_do = 1.0, 1.0

    # ── continuous dynamics ───────────────────────────────────────────────
    def _f(self, x):
        X, G, E, m1, m2 = x
        X = max(X, 0.0); G = max(G, 0.0); E = max(E, 0.0)
        q_ox, q_fm, muE, mu = carbon_fluxes(G, E, m1, m2, self.kin,
                                            self._phi_g, self._f_do)
        q_Eox = muE / self.Yex
        return np.array([
            mu * X,
            -(q_ox + q_fm) * X,
            (q_fm * self.kin["YGE_FERM"] - q_Eox) * X,
            0.0, 0.0,
        ])

    def _jac(self, x, eps=1e-6):
        """Forward-difference Jacobian — the carbon split is non-smooth."""
        f0 = self._f(x)
        J = np.zeros((N_STATE, N_STATE))
        for j in range(N_STATE):
            xp = np.array(x, dtype=float)
            xp[j] += eps
            J[:, j] = (self._f(xp) - f0) / eps
        return J

    def _combined_ode(self, t, s):
        x = s[:N_STATE]
        P = s[N_STATE:].reshape(N_STATE, N_STATE)
        F = self._jac(x)
        return np.concatenate([self._f(x), (F @ P + P @ F.T + self.Q).flatten()])

    def predict(self, dt_h: float, phi_g: float, f_do: float):
        """Advance state + covariance over dt_h hours at fixed (phi_g, f_do)."""
        self._phi_g, self._f_do = phi_g, f_do
        s0 = np.concatenate([self.x, self.P.flatten()])
        sol = solve_ivp(self._combined_ode, (0.0, dt_h), s0,
                        rtol=1e-6, atol=1e-9, method="RK45")
        s = sol.y[:, -1]
        self.x = s[:N_STATE]
        # Concentrations can't go negative — clamp X, G, E as the plant RK4 does,
        # so the glucose estimate doesn't dip below zero as G→0.
        for i in (0, 1, 2):
            self.x[i] = max(self.x[i], 0.0)
        self.P = s[N_STATE:].reshape(N_STATE, N_STATE)

    def correct(self, z_ethanol: float):
        """Ethanol measurement update (Joseph form). Returns (innovation, S).

        Same signature as EKF8.correct so the simulator can call either filter
        without branching — though S is only *used* in Design A, where the
        innovation feeds the CUSUM.
        """
        S = H @ self.P @ H.T + self.R
        K = self.P @ H.T @ np.linalg.inv(S)
        innovation = z_ethanol - self.x[2]
        self.x = self.x + (K @ np.array([[innovation]])).flatten()
        for i in (0, 1, 2):
            self.x[i] = max(self.x[i], 0.0)   # keep X, G, E physical (≥ 0)
        IKH = self.I5 - K @ H
        self.P = IKH @ self.P @ IKH.T + K @ self.R @ K.T
        return float(innovation), float(S[0, 0])

    def sigma(self):
        """Standard deviations of [X, G, E] for optional uncertainty bands."""
        d = np.clip(np.diag(self.P)[:3], 0.0, None)
        return np.sqrt(d)
