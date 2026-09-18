"""
cusum.py — the two-sided CUSUM detector.

Ported verbatim from fault_detection_CU.ipynb, Cell 14, so the UI and the
notebook visibly run the same recursion. It carries only two floats: no history
is needed, because the update depends on nothing older than the previous arm.

    S+ <- max(0, S+ + z - k)      ethanol ABOVE prediction
    S- <- max(0, S- - z - k)      ethanol BELOW prediction  <- the fault direction

`z` is the innovation normalised by the filter's own predicted spread sqrt(S),
so a fixed threshold means the same thing at every step.
"""
from __future__ import annotations

# Both in units of sqrt(S), which is what the detector sees (sqrt(S) = 0.329 g/L
# at R = 0.1).



K_SLACK = 0.15
H_CUSUM = 1.5


def cusum_step(sp: float, sn: float, z: float,
               k: float = K_SLACK, h: float = H_CUSUM):
    """ONE CUSUM update — this is the detector as it runs on the rig.

    Returns (S_plus, S_minus, fired). The arms are NOT reset after firing: this
    study asks when a fault was FIRST seen, so the arms stay a monotone evidence
    trace. Reset both to 0.0 here for continuous monitoring that can catch a
    second, later fault.
    """
    sp = max(0.0, sp + z - k)
    sn = max(0.0, sn - z - k)
    return sp, sn, (sp > h or sn > h)
