"""test_mmae_isolation.py — the MMAE bank must name the fault that was injected.

The regression this pins down: `f_pH` used the Rosso n = 2 form on the
asymmetric triple (3.0, 5.0, 12.0), whose denominator has a root at ~3.84. The
clip to [0, 1] turned that pole into a step — f_pH was exactly 0 below it and
exactly 1 above — so a stuck acid pump could only ever kill the culture
INSTANTLY, in the same shape as an oxygen crash. Both drive phi to exactly 0,
ethanol cannot tell them apart, and the bank read `acid_pump_stuck` as
`agitation_aer`. See the note above `make_cardinal_modifier` in core.py.

Only faults that are detectable at their hypothesis severity are asserted on:
the `*_off` pump faults do not trip the CUSUM in the robust default organism, so
there is no verdict to check. `acid_pump_off` never can — it is state-identical
to healthy, because a healthy batch never calls for acid (see faults.py).

Two scorings are asserted separately, because they answer different questions:
ETHANOL-ONLY is the honest soft-sensor question and the mode the f_pH bug
appeared in; PROBES adds the rig's own T/pH/DO instruments and is the only thing
that separates acid_pump_stuck from base_pump_stuck — those two saturate the same
kill, so ethanol is a provable coin flip (0.1 sigma apart).

Run: python3 -m backend.test_mmae_isolation
"""
from __future__ import annotations

from .plant import build_modifiers, resolve_params
from .simulator import Simulation
from .mmae import run_mmae


# fault id -> the injection that reaches the bank's own hypothesis severity
CASES = {
    "heater_off":         {"t_fault": 0.5, "Q_pct": 0.0},
    "heater_clog":        {"t_fault": 0.5, "Q_pct": 6.0},
    "acid_pump_stuck":    {"t_fault": 0.5, "F_A_stuck": 210.0},
    "agitation_aer_stop": {"t_fault": 0.5, "kla_pct": 0.0},
    "agitation_aer_clog": {"t_fault": 0.5, "kla_pct": 1.0},
}

# ethanol alone cannot separate these from their loop partner — see the module
# docstring. They are asserted in PROBES mode only.
PROBE_ONLY = {
    "base_pump_stuck": {"t_fault": 0.5, "F_B_stuck": 210.0},
}


def test_f_pH_has_a_gradient():
    """The shape check, independent of any filter: no 0 -> 1 jump."""
    _, f_pH, _ = build_modifiers(resolve_params(None))
    assert abs(f_pH(5.0) - 1.0) < 1e-12, "f_pH must be 1 at the optimum"
    # the old step sat at ~3.84; these must now be strictly ordered and interior
    lo, mid, hi = f_pH(3.0), f_pH(3.5), f_pH(4.0)
    assert 0.0 < lo < mid < hi, f"f_pH is not graded: {lo}, {mid}, {hi}"


def _verdict(fid, params, channel):
    sim = Simulation({})
    ok, msg = sim.inject_fault({"fault": fid, **params})
    assert ok, msg
    while not sim.done:
        sim.step()
    res = run_mmae(sim)
    assert res["state"] == "ok", f"{fid}: no verdict ({res.get('state')})"
    return res["scores"][channel]["verdict"]


def test_bank_names_the_injected_fault():
    for fid, params in CASES.items():
        v = _verdict(fid, params, "ethanol")
        assert v == fid, f"{fid} was diagnosed as {v} (ethanol-only)"


def test_probes_break_the_pH_pump_tie():
    for fid, params in PROBE_ONLY.items():
        v = _verdict(fid, params, "probes")
        assert v == fid, f"{fid} was diagnosed as {v} (with probes)"


def test_acid_pump_off_is_a_null_fault():
    """A healthy batch barely calls for acid, so killing that pump is invisible.

    In the notebook the pH law is pure PROPORTIONAL, so it never overshoots the
    setpoint, the acid pump is identically off, and acid_pump_off is null to
    0.0e+00 in every state. The app runs the PI law instead, whose integral term
    overshoots to pH 5.0017 — so the acid pump does fire, briefly, and the fault
    is not bit-identical. It is still null in the only sense that matters: the
    ethanol channel moves by 0.0000 sigma, i.e. nothing a sensor could ever see.
    """
    import numpy as np
    base = Simulation({"ekf": {"on": False}})
    while not base.done:
        base.step()
    hurt = Simulation({"ekf": {"on": False}})
    hurt.inject_fault({"fault": "acid_pump_off", "t_fault": 0.5})
    while not hurt.done:
        hurt.step()
    got = {k: np.abs(np.array([r[k] for r in base.rows])
                     - np.array([r[k] for r in hurt.rows])).max()
           for k in ("X", "G", "E", "T", "pH", "DO")}
    sigma_E = base.params["ekf"]["sigma_E_sensor"]
    assert got["E"] / sigma_E < 1e-3, (
        f"acid_pump_off moved ethanol by {got['E'] / sigma_E:.4f} sigma - "
        f"it is supposed to be undetectable")
    # pH is the one channel it touches at all, and only through PI overshoot
    assert got["pH"] < 1e-2, f"acid_pump_off moved pH by {got['pH']:.2e}"


if __name__ == "__main__":
    test_f_pH_has_a_gradient()
    print("f_pH gradient: ok")
    test_bank_names_the_injected_fault()
    print(f"ethanol-only isolation: {len(CASES)}/{len(CASES)} correct")
    test_probes_break_the_pH_pump_tie()
    print(f"probes break the pH-pump tie: {len(PROBE_ONLY)}/{len(PROBE_ONLY)} correct")
    test_acid_pump_off_is_a_null_fault()
    print("acid_pump_off is undetectable on ethanol (0.0000 sigma): ok")
