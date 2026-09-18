"""
faults.py — fault registry and control-context resolution.

Each *active* fault modifies the actuator context (heater on/off, effective UA,
pump flows) for t >= t_fault. The simulator computes the healthy context every
tick, then calls `apply_faults` to fold in whatever faults are scheduled.

The registry is FOUR ACTUATOR LOOPS x TWO FAILURE MODES, matching
`7_mmae_fault_isolation.ipynb` (A.4):

  loop      dead                     running-but-wrong
  ------    ---------------------    -------------------------------------
  heater    heater_off  Q_pct=0      heater_clog        Q_pct=6
  acid      acid_pump_off            acid_pump_stuck    F_A_stuck=210 mL/h
  base      base_pump_off            base_pump_stuck    F_B_stuck=210 mL/h
  oxygen    agitation_aer_stop       agitation_aer_clog kla_pct=1

WHAT THE *INJECTOR* OFFERS IS SHORTER THAN THAT. For the heater and the oxygen
loop, "dead" and "running-but-wrong" are the SAME code path: both resolve to a
single continuous severity (Q_pct, kla_pct), and 0 % is just the bottom rung of
it. Offering them as two catalogue entries with the same slider was asking the
user to pick a rung twice, so the catalogue carries one entry per loop —
"Heater degraded", "Agitation / aeration degraded" — and the % field says which
rung. The pump loops keep two entries because dead (F = 0) and stuck-ON
(F = saturation) are genuinely different code paths, not two values of one
knob.

The MMAE BANK still carries all eight: a hypothesis has to commit to a
severity, and heater_off vs heater_clog is exactly the "how dead" question the
bank answers. `_clog` ids therefore stay in FAULT_DEFAULTS and stay resolvable
in `apply_faults` — they are simply not offered as separate things to inject.
See `mmae.candidate_causes`.

AGITATION AND AERATION ARE ONE MECHANISM, NOT TWO. Both reach the biology only
through kLa, so any (rpm, air-flow) pair giving the same kLa produces a
bit-identical trajectory in all 8 states — the notebook proved agitation_stop
== aeration_stop to 1e-9. They are one hypothesis parameterised two ways, and
only an instrument on the actuator itself (a tachometer, an air mass-flow meter)
can separate them. Severity is therefore a kLa RATIO, not an rpm.

`vessel_leak` was removed: it is a passive thermal fault, and at the severe rung
it lands within 0.2 % of the same biomass loss as heater_off, so it added a
permanent tie to the bank without adding a mechanism the rig can act on.

Coming soon (declared for the UI, not executable):
  * whole Design and Ageing categories
"""
from __future__ import annotations

from .plant import compute_kLa


# Per-fault default parameters (user-overridable from the UI).

FAULT_DEFAULTS = {
    "heater_off":         {"Q_pct": 0.0},        # % of heater power still delivered
    "heater_clog":        {"Q_pct": 6.0},
    "acid_pump_off":      {},                     # F_A forced to 0 (dead)
    "acid_pump_stuck":    {"F_A_stuck": 210.0},  # mL/h, = pump saturation
    "base_pump_off":      {},                     # F_B forced to 0 (dead)
    "base_pump_stuck":    {"F_B_stuck": 210.0},
    "agitation_aer_stop": {"kla_pct": 0.0},      # % of healthy kLa still delivered
    "agitation_aer_clog": {"kla_pct": 1.0},
}


# Catalogue served to the frontend. `enabled` gates the UI control.
_T_FAULT = {"key": "t_fault", "label": "Fault time [h]", "default": 0.5}

FAULT_CATALOGUE = [
    {
        "category": "Operation",
        "enabled": True,
        "faults": [
            {"id": "heater_off", "label": "Heater degraded", "enabled": True,
             "params": [_T_FAULT,
                        {"key": "Q_pct", "label": "Power left [%]", "default": 0.0}]},
            {"id": "acid_pump_off", "label": "Acid pump OFF (dead)", "enabled": True,
             "params": [_T_FAULT]},
            {"id": "acid_pump_stuck", "label": "Acid pump stuck ON", "enabled": True,
             "params": [_T_FAULT,
                        {"key": "F_A_stuck", "label": "Stuck flow [mL/h]", "default": 210.0}]},
            {"id": "base_pump_off", "label": "Base pump OFF (dead)", "enabled": True,
             "params": [_T_FAULT]},
            {"id": "base_pump_stuck", "label": "Base pump stuck ON", "enabled": True,
             "params": [_T_FAULT,
                        {"key": "F_B_stuck", "label": "Stuck flow [mL/h]", "default": 210.0}]},
            {"id": "agitation_aer_stop", "label": "Agitation / aeration degraded", "enabled": True,
             "params": [_T_FAULT,
                        {"key": "kla_pct", "label": "kLa left [%]", "default": 0.0}]},
        ],
    },
    {"category": "Design",  "enabled": False, "faults": []},
    {"category": "Ageing",  "enabled": False, "faults": []},
]


def normalize_fault(fault: dict) -> dict:
    """Fill a fault request with defaults for its type."""
    fid = fault.get("fault") or fault.get("id")
    merged = dict(FAULT_DEFAULTS.get(fid, {}))
    merged.update({k: v for k, v in fault.items() if k not in ("fault", "id")})
    merged["fault"] = fid
    return merged


def apply_faults(active_faults: list[dict], t: float, params: dict,
                 base_F_A: float, base_F_B: float):
    """Resolve the actuator context at time t given the scheduled faults.

    Returns (F_A, F_B, heater_frac, UA_eff, kLa_eff).  heater_frac is the
    fraction of commanded heater power actually delivered (1.0 = healthy).
    """
    F_A, F_B = base_F_A, base_F_B
    heater_frac = 1.0
    UA_eff = params["thermal"]["UA"]
    ox = params["oxygen"]
    V = params["thermal"]["V"]
    kLa_eff = compute_kLa(ox["N_rpm"], ox["Q_air_Lmin"], V)

    for f in active_faults:
        if t < f.get("t_fault", 0.0):
            continue
        fid = f["fault"]
        if fid in ("heater_off", "heater_clog"):
            # 0 % = dead heater, 100 % = healthy, in between = underpowered.
            heater_frac = min(max(f.get("Q_pct", 0.0), 0.0), 100.0) / 100.0
        elif fid == "acid_pump_off":
            F_A = 0.0
        elif fid == "acid_pump_stuck":
            F_A = f.get("F_A_stuck", 210.0)
        elif fid == "base_pump_off":
            F_B = 0.0
        elif fid == "base_pump_stuck":
            F_B = f.get("F_B_stuck", 210.0)
        elif fid in ("agitation_aer_stop", "agitation_aer_clog"):
            # Severity is a RATIO of the healthy kLa, because kLa is the only
            # thing the plant sees — see the module docstring.
            kLa_eff *= min(max(f.get("kla_pct", 0.0), 0.0), 100.0) / 100.0

    return F_A, F_B, heater_frac, UA_eff, kLa_eff
