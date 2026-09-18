"""
mmae.py — Multiple-Model Adaptive Estimation: *which* fault happened?

The healthy-shadow EKF + CUSUM says **that** something is
wrong and roughly **when**. It says nothing about **what**. This module adds the
what, as the five steps of `MMAE.ipynb`:

  1. the shadow alarms on the ethanol innovation          (`simulator.step`)
  2. walk the firing CUSUM arm back to where it left zero  -> onset t0_hat
  3. build one EKF per candidate cause, each *believing* its own fault started
     at t0_hat
  4. feed all of them the SAME recorded ethanol samples
  5. compare: the Bayes recursion over their innovation sequences

A hypothesis filter is not a second EKF implementation. It is the healthy shadow
with a non-empty fault list handed to the very same `Simulation._advance` /
`Simulation._rhs_for` — so "the hypothesis believes fault f started at t0" is
true by construction. Before its own onset a hypothesis IS the healthy shadow,
which is what makes the prefix sharing below exact rather than an approximation.

Two rules carried over from the notebook, both measured there rather than
assumed:

  * **Onset ladder, backwards only.** The CUSUM onset is the weak link and its
    error is asymmetric: guessing early is nearly free (a hypothesis is just the
    healthy model before its onset), guessing late is fatal — a late onset always
    favours whichever cause acts fastest. So every cause is run at t0 and at
    30-minute steps back to t0-120 min; see LADDER_H for why the floor is -120
    and not -60.
  * **Profile, do not marginalise.** Collapse the onset axis with max, not sum.
    Summing over onsets rewards *slow* faults, which degrade gracefully as t0
    moves and so accumulate mass across many onsets, whatever is true. It
    produces confident wrong answers; the max does not.

Scoring uses the WHOLE batch, not a per-onset window: once candidates carry
different onsets, a per-onset window would reward the latest candidate simply
for having fewer samples to explain. Pre-onset samples are identical across
hypotheses and cancel, so this costs nothing.
"""
from __future__ import annotations

import math

import numpy as np

from .ekf8 import EKF8
from .faults import FAULT_CATALOGUE, FAULT_DEFAULTS, normalize_fault

# Onsets tried per cause, as offsets [h] from the CUSUM estimate. Backwards
# only — see the module docstring.
#
# DEPTH GOES TO -120 min, NOT -60. An earlier note here claimed depth saturated
# at -60; that was measured on a registry of faster faults. The oxygen-loop
# faults alarm late — their early innovations sit under the CUSUM slack k and
# keep resetting the arm, so the onset estimate lands up to +105 min after the
# truth. With a -60 min floor the ladder cannot reach back to the real onset and
# `agitation_aer_clog` is then read as `agitation_aer_stop` with p = 1.000 once
# the probes are in the likelihood — a confident error produced entirely by the
# floor. (7_mmae_fault_isolation.ipynb, C.2)
LADDER_H = (0.0, -0.5, -1.0, -1.5, -2.0)

# Two causes closer than this in total log-likelihood are not distinguishable by
# this data: 2.3 nats is 10:1 odds. An honest tie beats a confident error, so
# the verdict reports every cause inside the margin rather than the argmax.
TIE_MARGIN_NATS = 2.3

# Below this the two hypotheses produced *numerically identical* innovation
# sequences — a degeneracy in the plant (the notebook's provable
# agitation_stop == aeration_stop), not a shortfall of the estimator.
DEGENERATE_NATS = 1e-6

# Separation, in units of the filter's own spread sqrt(S), below which a
# hypothesis is indistinguishable from the null and its tie is an observability
# statement rather than a result. The notebook's rule: "anything below ~2 sigma
# is a tie no estimator can break."
SILENT_SIGMA = 2.0

# Minutes of batch AFTER the alarm before a verdict stops being flagged
# provisional. Below the threshold the UI withholds the cause NAME.
#
# THIS IS A RELEASE POINT, NOT A CERTAINTY GUARANTEE — do not describe it as
# "the verdict has settled". Measured over 30 trials (6 alarming scenarios x 5
# sensor-noise seeds, fault at t = 2 h, robust organism; see
# thesis_ekf_digital_twin/evaluation_accuracy_windows.py), POINTWISE rung
# accuracy of the ethanol scoring is:
#
#     offset   0     30    60    90   120   150   240   min after alarm
#     acc     40 %  60 %  77 %  83 %  83 %  97 % 100 %
#
# So at 120 one ethanol diagnosis in six is still wrong. 95 % needs 150 min and
# 100 % needs 240. 120 is KEPT anyway, deliberately: the panel always renders the
# tie set beside the name, so a wrong argmax inside a reported tie is a shortlist
# rather than a false diagnosis, and waiting to 150 would withhold a correct name
# from most faults for another 30 min.
#
# An earlier SINGLE-SEED measurement (8_evaluation.ipynb: mean 35, max 90) made
# this look ~3.4x conservative. That seed was lucky. The straggler across seeds is
# base_pump_stuck — 20 % correct from the alarm through 120 min.
#
# The accuracy curve is NOT MONOTONE (83 % at 90, 80 % at 105, 83 % at 120):
# verdicts flip near the tie margin, which is evidence that a fixed clock is the
# wrong instrument. Gating on verdict STABILITY (argmax unchanged across the last
# N re-scores) would measure the real thing. Open item, not implemented.
ETHANOL_SETTLE_MIN = 300.0

# The same, for the probe scoring. Probes separate causes by WHICH channel moved
# and in WHICH DIRECTION, not by accumulated curvature, so they need far less
# evidence: on the LOOP question (which actuator broke — the actionable one) they
# are at 100 % from the alarm itself, in every scenario and every seed.
#
# On the RUNG question they have a CEILING of ~83 % that no window ever lifts:
# the two oxygen rungs mutually confuse (agitation_aer_stop right in 2 of 5 seeds,
# agitation_aer_clog in 3 of 5) because a DO probe reads the same floor under
# kLa = 0 % and kLa = 1 % and is reduced to guessing between them. There is no
# "time to 95 %" for the probe rung — only a ceiling below it.
#
# 30 is therefore conservative for the loop question and irrelevant to the rung
# one. It is kept as a small settling allowance rather than tuned to either.
PROBE_SETTLE_MIN = 30.0

NULL = "healthy"

# The severity each hypothesis commits to. These are the notebook's registry
# values (7_mmae_fault_isolation.ipynb, A.4) and they now MATCH the injection
# defaults in faults.py, because the registry itself was rebuilt as four loops x
# two failure modes — the severity IS the mechanism here, so the bank and the
# injector no longer need to disagree.
#
# The two that still need care are the "running-but-wrong" rungs. Neither is a
# round number, and neither is mild:
#   heater_clog Q_pct=6   the P-controller draws only 144 W of 630 W at steady
#                         state, so anything above ~23 % is invisible
#   agitation_aer_clog    kla_pct=1: the inherited oxygen stoichiometry (see
#                         core.py) makes demand ~53x too small, so DO holds near
#                         saturation until kLa is nearly gone
# A hypothesis has to commit to a severity; the bank asks "which mechanism", not
# "which mechanism, how hard". A matching mechanism at a different severity still
# scores best, with a smaller margin.
HYPOTHESIS_SEVERITY = {
    "heater_off":         {"Q_pct": 0.0},
    "heater_clog":        {"Q_pct": 6.0},
    "acid_pump_off":      {},
    "acid_pump_stuck":    {"F_A_stuck": 210.0},
    "base_pump_off":      {},
    "base_pump_stuck":    {"F_B_stuck": 210.0},
    "agitation_aer_stop": {"kla_pct": 0.0},
    "agitation_aer_clog": {"kla_pct": 1.0},
}

# THE BANK IS WIDER THAN THE INJECTOR, AND ON PURPOSE. The catalogue in
# faults.py offers ONE entry per continuously-parameterised loop (heater,
# oxygen) because there the user sets the severity with a % field — asking them
# to also pick "dead" vs "underpowered" from a list was picking the same rung
# twice. The bank has no such field: each hypothesis must commit to a severity
# before it can be integrated, so both rungs stay here. Dropping them would not
# simplify the diagnosis, it would delete the answer "the heater is weak but not
# dead" from the set of answers the panel can give.
#
# Labels are spelled out here rather than looked up from the catalogue for the
# same reason: `heater_clog` no longer HAS a catalogue entry to borrow a label
# from.
HYPOTHESIS_LABELS = {
    "heater_off":         "Heater dead",
    "heater_clog":        "Heater underpowered",
    "acid_pump_off":      "Acid pump OFF (dead)",
    "acid_pump_stuck":    "Acid pump stuck ON",
    "base_pump_off":      "Base pump OFF (dead)",
    "base_pump_stuck":    "Base pump stuck ON",
    "agitation_aer_stop": "Agitation / aeration STOP",
    "agitation_aer_clog": "Agitation / aeration restricted",
}


def candidate_causes() -> dict[str, dict]:
    """The bank: every hypothesis, at its committed severity.

    Gated by the catalogue only at CATEGORY level (an "Operation" bank exists,
    the Design and Ageing categories do not yet), never per fault id — see the
    note above HYPOTHESIS_LABELS.
    """
    if not any(g.get("enabled") for g in FAULT_CATALOGUE
               if g["category"] == "Operation"):
        return {}
    causes: dict[str, dict] = {}
    for fid, sev in HYPOTHESIS_SEVERITY.items():
        params = dict(FAULT_DEFAULTS.get(fid, {}))
        params.update(sev)
        causes[fid] = {"label": HYPOTHESIS_LABELS.get(fid, fid), "params": params}
    return causes


# ── step 2: where did it start? ──────────────────────────────────────────────
def onset_estimate(meas_rows: list[dict], alarm_idx: int) -> float:
    """Last time the FIRING CUSUM arm was zero before the alarm.

    Unreliable in both directions, which is exactly why the ladder exists: a
    raised arm means nothing on its own (on healthy data an arm leaves zero
    ~11x per batch on noise alone), so the walk back can run through a pre-fault
    excursion and land early; a slow fault whose early innovations sit below the
    slack k lands late.
    """
    i = alarm_idx
    arm_key = "Sn" if meas_rows[i]["Sn"] >= meas_rows[i]["Sp"] else "Sp"
    j = i
    while j > 0 and meas_rows[j - 1][arm_key] > 0.0:
        j -= 1
    return float(meas_rows[j]["time_h"])


# ── steps 3-4: one filter per hypothesis, same measurements ─────────────────
class _Hyp:
    """A hypothesis filter: the shadow EKF plus the fault list it believes in."""

    def __init__(self, sim, faults: list[dict]):
        self.ekf = EKF8(sim.params)
        self.faults = faults

    def snapshot(self) -> dict:
        e = self.ekf
        return {"x": e.x.copy(), "P": e.P.copy(),
                "prev_x": e.prev_x.copy(), "pi_int": e.pi_int}

    def restore(self, s: dict):
        e = self.ekf
        e.x, e.P = s["x"].copy(), s["P"].copy()
        e.prev_x, e.pi_int = s["prev_x"].copy(), s["pi_int"]


def _replay(sim, hyp: _Hyp, grid, i0, i1, meas_by_frame, ll, tr, probes,
            snap_at=None):
    """Advance `hyp` over frames [i0, i1) of the run's own frame grid.

    This is `Simulation.step`'s Design-A branch with the hypothesis' fault list
    in place of the empty one, and the RECORDED ethanol sample in place of a
    fresh draw. `tr` receives this hypothesis' trace: the innovation that was
    scored, the ethanol it predicted, and the T/pH/DO it predicted — the
    channels ethanol cannot see, and therefore where the causes visibly differ.
    `snap_at` maps a frame index to a dict that receives a state snapshot taken
    BEFORE that frame is integrated.

    Both scores are accumulated in the SAME pass, because they differ only in
    the likelihood, not in the trajectory: `ll["E"][k]` is the ethanol-only
    log-likelihood of sample k, `ll["P"][k]` the extra term contributed by the
    rig's T/pH/DO instruments. Scoring both costs one replay, so the UI can
    report the soft-sensor answer and the instrumented answer side by side
    instead of making the user choose blind.
    """
    e = hyp.ekf
    for i in range(i0, i1):
        if snap_at is not None and i in snap_at:
            snap_at[i].update(hyp.snapshot())
        t0 = grid[i - 1] if i > 0 else 0.0
        t1 = grid[i]
        e.x, e.pi_int = sim._advance(e.x, t0, t1, hyp.faults, e.pi_int)
        e.predict_cov(t1 - t0, sim._rhs_for(e.x, e.pi_int, hyp.faults, t1))
        m = meas_by_frame.get(i)
        if m is not None:
            k = m["k"]
            tr["E"][k] = e.x[2]
            tr["T"][k], tr["pH"][k], tr["DO"][k] = e.x[5], e.x[6], e.x[7]
            nu, S = e.correct(m["z"])
            tr["nu"][k] = nu
            ll["E"][k] = -0.5 * nu * nu / S - 0.5 * math.log(2.0 * math.pi * S)
            # The rig's own T/pH/DO instruments. This is where the causes
            # actually differ — a dead heater and a jammed acid pump end in the
            # same stalled ethanol curve, but at 22 °C vs 28 °C — so it is what
            # breaks the ethanol-only tie.
            acc = 0.0
            for ch in ("T", "pH", "DO"):
                s = probes["sigma"][ch]
                d = probes[ch][k] - tr[ch][k]
                acc += -0.5 * d * d / (s * s) - 0.5 * math.log(2.0 * math.pi * s * s)
            ll["P"][k] = acc
    if snap_at is not None and i1 in snap_at:
        snap_at[i1].update(hyp.snapshot())


def _new_trace(n: int) -> dict:
    return {k: np.zeros(n) for k in ("nu", "E", "T", "pH", "DO")}


# ── step 5: Bayes over the innovation sequences ─────────────────────────────
def _score(names, runs, key, n):
    """Collapse one likelihood channel into a verdict: posterior, ties, margin.

    `key` selects which of the two scores accumulated in `_replay` is used —
    "E" for the ethanol-only bank the notebook studies, "EP" for the same bank
    with the rig's T/pH/DO instruments added to the likelihood. Everything else
    (the hypotheses, the onsets, the trajectories) is shared.
    """
    totals = np.array([runs[c][key]["total"] for c in names])
    order = [int(j) for j in np.argsort(totals)[::-1]]

    p = np.exp(totals - totals.max())          # uniform prior over causes
    p /= p.sum()

    # Running posterior — column k is the belief after sample k, each cause at
    # its own profiled-best onset.
    L = np.vstack([runs[c][key]["ll"] for c in names])
    C = np.cumsum(L, axis=1)
    C -= C.max(axis=0)
    P = np.exp(C)
    P /= P.sum(axis=0)

    margin = (float(totals[order[0]] - totals[order[1]])
              if len(order) > 1 else None)
    tied = [names[j] for j in order
            if totals[order[0]] - totals[j] < TIE_MARGIN_NATS]
    # Two hypotheses with numerically identical scores are degenerate ON THIS
    # DATA: no estimator and no threshold can separate them, because the plant
    # walked the same trajectory under both. The catalogue no longer CONTAINS
    # such a pair — agitation and aeration used to be two entries and were
    # provably equal to 1e-9, which is why they are now one — but the check
    # stays, because it is what proved that and would catch the next one.
    degenerate = [[a, b] for i, a in enumerate(names) for b in names[i + 1:]
                  if abs(runs[a][key]["total"] - runs[b][key]["total"])
                  < DEGENERATE_NATS]
    return {
        "verdict": names[order[0]],
        "verdict_label": runs[names[order[0]]]["label"],
        "tied": tied,
        "margin_nats": margin,
        "degenerate": degenerate,
        "order": [names[j] for j in order],
        "p": {c: float(p[j]) for j, c in enumerate(names)},
        "loglik": {c: float(totals[j]) for j, c in enumerate(names)},
        "post": {c: [float(v) for v in P[j]] for j, c in enumerate(names)},
        "onset": {c: runs[c][key]["onset"] for c in names},
    }


def run_mmae(sim, ladder_h=LADDER_H) -> dict:
    """The whole pipeline: onset -> bank -> posterior -> verdict.

    Runs against a *snapshot* of the run so far, so it is safe to call from a
    worker thread while the simulation keeps streaming. Returns a JSON-ready
    dict; `{"state": ...}` says why there is no verdict when there is none.

    Both scorings are returned, from one set of replays: `scores["ethanol"]` is
    the notebook's soft-sensor bank, `scores["probes"]` the same bank with the
    T/pH/DO instruments in the likelihood.
    """
    meas = list(sim.meas_rows)
    alarm_idx = sim.alarm_idx
    if sim.old_design or not sim.ekf_on:
        return {"state": "unavailable",
                "note": "Root-cause isolation needs the Design-A healthy shadow."}
    if alarm_idx is None:
        return {"state": "no_alarm",
                "note": "No alarm yet - nothing to isolate."}

    grid = [r["time_h"] for r in sim.rows]
    n_frames = len(grid)
    meas = [m for m in meas if m["frame"] < n_frames]
    n = len(meas)
    meas_by_frame = {m["frame"]: {"z": m["z"], "k": k} for k, m in enumerate(meas)}
    t_now = grid[-1]

    ek = sim.params["ekf"]
    probes = {"sigma": {"T": ek["sigma_T_probe"], "pH": ek["sigma_pH_probe"],
                        "DO": ek["sigma_DO_probe"]}}
    for ch in ("T", "pH", "DO"):
        probes[ch] = np.array([sim.rows[m["frame"]][ch] for m in meas])

    t0_hat = onset_estimate(meas, alarm_idx)
    # Ladder onsets, snapped DOWN to a frame boundary so a hypothesis can be
    # forked from a shared healthy prefix exactly rather than nearly.
    onsets, frames = [], []
    for off in ladder_h:
        t = max(0.0, t0_hat + off)
        fi = int(np.searchsorted(grid, t - 1e-12))
        fi = min(max(fi, 0), n_frames - 1)
        if fi not in frames:
            frames.append(fi)
            onsets.append(grid[fi - 1] if fi > 0 else 0.0)

    # ── the null hypothesis, which is also the shared prefix ────────────────
    # It is the healthy shadow, replayed. Every faulted hypothesis is identical
    # to it until its own onset, so one pass produces both the null's score and
    # the fork points for all the others.
    null = _Hyp(sim, [])
    ll_null = {"E": np.zeros(n), "P": np.zeros(n)}
    tr_null = _new_trace(n)
    snaps = {fi: {} for fi in frames}
    _replay(sim, null, grid, 0, n_frames, meas_by_frame, ll_null, tr_null,
            probes, snap_at=snaps)

    def _pack(label, onset, ll, tr):
        """One hypothesis run, scored both ways."""
        e, ep = ll["E"], ll["E"] + ll["P"]
        return {"label": label, "tr": tr,
                "E":  {"onset": onset, "ll": e,  "total": float(e.sum())},
                "EP": {"onset": onset, "ll": ep, "total": float(ep.sum())}}

    causes = candidate_causes()
    runs: dict[str, dict] = {NULL: _pack("No fault (null)", None, ll_null, tr_null)}
    for cid, meta in causes.items():
        best = None
        for fi, t_on in zip(frames, onsets):
            fault = normalize_fault({"fault": cid, **meta["params"],
                                     "t_fault": t_on})
            hyp = _Hyp(sim, [fault])
            hyp.restore(snaps[fi])
            # pre-onset samples: identical to the null, by construction
            ll = {k: v.copy() for k, v in ll_null.items()}
            tr = {k: v.copy() for k, v in tr_null.items()}
            _replay(sim, hyp, grid, fi, n_frames, meas_by_frame, ll, tr, probes)
            cand = _pack(meta["label"], t_on, ll, tr)
            # PROFILE over the onset axis (max), never marginalise (sum). The
            # two scorings profile independently, so each reports the onset that
            # best explains the evidence IT is allowed to see; the traces kept
            # are the ethanol-profiled ones, which is what the proof plots show.
            if best is None:
                best = cand
            else:
                for key in ("E", "EP"):
                    if cand[key]["total"] > best[key]["total"]:
                        best[key] = cand[key]
                if cand["E"]["total"] > best["E"]["total"]:
                    best["tr"] = cand["tr"]
        runs[cid] = best

    names = list(runs)
    scores = {"ethanol": _score(names, runs, "E", n),
              "probes":  _score(names, runs, "EP", n)}

    # ── how much batch is this verdict actually standing on? ────────────────
    # main.py fires the bank the instant the CUSUM latches, so by default this
    # runs on the [0, alarm] window. 7_mmae_fault_isolation.ipynb C.4 measures
    # what that window is worth, and the two channels answer differently:
    #
    #   ethanol only  the alarm fires as soon as the innovation is improbable,
    #                 which is well before the drift has a recognisable SHAPE.
    #                 Only 2 of 6 causes are right at that moment; the rest need
    #                 another 15-120 min of batch before the verdict settles, and
    #                 the acid/base pump pair never settles at all.
    #   with probes   every cause is already right at the alarm, at p = 1.000.
    #                 T/pH/DO differ in which channel moved and in which
    #                 direction, which is one sample's worth of evidence.
    #
    # So an ethanol verdict taken at the alarm is PROVISIONAL, and saying so is
    # the difference between a diagnosis and a guess that will quietly change.
    mins_after_alarm = (t_now - float(sim.alarm_t)) * 60.0
    scores["ethanol"]["provisional"] = bool(mins_after_alarm < ETHANOL_SETTLE_MIN)
    scores["probes"]["provisional"] = bool(mins_after_alarm < PROBE_SETTLE_MIN)
    for k in scores:
        scores[k]["evidence_min_after_alarm"] = float(mins_after_alarm)

    # Screening (the notebook's step 0): how far each hypothesis' predicted
    # ethanol departs from the null's, in units of the filter's own spread
    # sqrt(S). A cause that never separates from "no fault" cannot be isolated
    # from ethanol at any threshold — the signal is absent, not buried — so its
    # tie is a property of the plant, labelled as such rather than read as a
    # result. (It may still be wide open to the probes. The sharpest case is
    # acid_pump_stuck vs base_pump_stuck in a pH-sensitive strain: both saturate
    # the SAME kill, so ethanol is a coin flip — it sees only |dpH| through
    # f_pH, never its sign — while the pH probe separates them by 3.5 pH units
    # and the probes column returns p = 1.000.)
    root_S = math.sqrt(max(sim.ekf.R[0, 0], 1e-12))
    sep = {c: float(np.abs(runs[c]["tr"]["E"] - tr_null["E"]).max() / root_S)
           for c in names}
    sep[NULL] = 0.0
    silent = [c for c in names if c != NULL and sep[c] < SILENT_SIGMA]

    return {
        "state": "ok",
        "alarm_t": float(sim.alarm_t),
        "t_onset": t0_hat,
        "onsets": onsets,
        "t_evidence": float(t_now),
        "n_samples": n,
        "ethanol_settle_min": ETHANOL_SETTLE_MIN,
        "probe_settle_min": PROBE_SETTLE_MIN,
        "default_mode": "probes" if ek.get("mmae_probes") else "ethanol",
        "silent": silent,
        "scores": scores,
        "t_meas": [m["time_h"] for m in meas],
        "z_meas": [m["z"] for m in meas],
        # What the rig's own instruments read at those sample times. Plotted
        # against each hypothesis' prediction in both modes — with the ethanol
        # scoring it is the visual check on a verdict the bank reached without
        # them.
        "meas": {ch: [float(v) for v in probes[ch]] for ch in ("T", "pH", "DO")},
        "causes": [
            {"id": c, "label": runs[c]["label"], "sep": sep[c],
             "silent": c in silent,
             # the evidence itself: what this hypothesis predicted, and how
             # wrong it was on the one channel the soft sensor measures
             "nu": [float(v) for v in runs[c]["tr"]["nu"]],
             "E":  [float(v) for v in runs[c]["tr"]["E"]],
             "T":  [float(v) for v in runs[c]["tr"]["T"]],
             "pH": [float(v) for v in runs[c]["tr"]["pH"]],
             "DO": [float(v) for v in runs[c]["tr"]["DO"]]}
            for c in names
        ],
    }
