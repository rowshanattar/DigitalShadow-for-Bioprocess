# Bioreactor Digital-Twin Simulator (local web UI)

Interactive front-end for the baker's-yeast EKF digital twin. Edit every
parameter, watch the concentration curves grow in real time (1 simulated hour ≈
5 s), overlay the EKF estimate, inject faults live, get a root cause for the
alarm, and download CSV + PNG.

Built on the notebook models in the repository root: `1_ekf_baseline_yd.ipynb` (EKF),
`4_fault_injection.ipynb` (plant + faults), `5_fault_detection.ipynb` (healthy shadow +
CUSUM) and `6_mmae_fault_isolation.ipynb` (fault isolation). The backend is a cell-by-cell
port, so the app and the thesis cannot disagree.

## Run

Requires Python 3.12 (3.10+ should work; the pins are verified on 3.12).

```bash
cd bioreactor_ui
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn backend.main:app                                # add --reload while developing
```

Open <http://localhost:8000>.

> **Run it locally only.** `uvicorn` binds `127.0.0.1` by default and that is the
> intended deployment: there is no authentication, no rate limiting, and the
> simulation config is validated only for *physics* plausibility (`build_modifiers`),
> not as untrusted input. A long `t_final` or a fast cadence is a CPU-bound run this
> app will happily start for anyone who can reach it. Do not pass `--host 0.0.0.0`
> or expose the port without putting a real front end and auth in front of it.

> Plotly is loaded from a CDN, so the first load needs internet. Everything
> else (simulation, EKF, faults, exports) runs entirely on your machine.

## Using it

1. **Setup** — pick cell line (Yeast) and reactor (Minifors 2). Other options are
   placeholders ("coming soon").
2. Tweak any parameter in the left panels (kinetics, initial state, thermal,
   pH, EKF, run settings). Defaults match the current notebooks.
3. **▶ Run** — the six panels animate. Enable the **EKF** toggle to overlay the
   filter's estimate (dashed); enable **±σ band** for uncertainty ribbons.
4. **Fault injection** — choose the *Operation* category, set a fault's time and
   parameters, and click **＋ inject**. Before a run this queues the fault; during
   a run it injects live (rejected if the time has already passed). Design and
   Ageing categories, the DO panel, and extra pumps/agitation are "coming soon".
5. **Root cause** — when the CUSUM alarms, an MMAE filter bank runs one EKF per
   candidate fault (each believing its own onset, walked back from the alarm)
   and the panel ranks the causes by posterior probability, with the proof:
   the posterior as it built up, each hypothesis' innovation, and what each
   predicted for ethanol and for T/pH/DO. Two tabs: **ethanol only** (the
   soft-sensor question — often an honest tie) and **+ T/pH/DO probes** (the
   rig's own instruments in the likelihood, which breaks it).

   The panel **does not name a cause until the evidence window has closed**. A
   verdict taken at the alarm has not settled — on ethanol alone only half the
   causes are right at that moment and the rest change over the next 35 min on
   average (90 max), so naming one and quietly replacing it later reads as the
   method being unreliable. Until then it shows a progress bar and the reason,
   with a *Show the provisional answer anyway* button. The run **re-diagnoses
   itself** when each window matures, so the name appears on its own; **↻
   Re-diagnose** forces a re-score at any time.
6. **⤓ CSV / ⤓ PNG** — after a run, download the standardized
   `time_h,T,pH,DO,G,E,X,mu_G,mu_E` table and the multi-panel figure (detector
   and root-cause rows included).

## Layout

```
backend/   core.py plant.py faults.py ekf.py ekf8.py cusum.py mmae.py simulator.py main.py
frontend/  index.html app.js styles.css
```

The physics/EKF live in `backend/` (pure, parameter-driven, ported cell-by-cell
from the notebooks); `main.py` streams frames over a WebSocket paced to the
5 s/sim-hour cadence.


## How it works, in order

What happens where and when, from page load to CSV. `DOCUMENTATION.md` is the
reference for each piece; this is the sequence that connects them.

### Phase 0 — page load, once

`boot()` makes one HTTP call to `/api/defaults`, which returns `DEFAULTS` (every
parameter), `FAULT_CATALOGUE`, `CELLS` and `REACTORS`, then builds the DOM from
them.

**Why fetch rather than hard-code:** the frontend holds no physics. Every
default, unit and fault definition lives in `core.py` and `faults.py`, so
changing a parameter there changes the UI with no JS edit. `CELLS`/`REACTORS`
carry `enabled` flags, so CHO and E. coli are listed but greyed — the UI
advertises the roadmap without pretending to implement it.

The default root-cause tab is read from the backend (`ekf.mmae_probes`), so the
notebooks and the UI agree out of the box. Static assets are served
`Cache-Control: no-store` on purpose: this is edited while running, and a
browser holding a stale `app.js` looks exactly like a broken backend.

### Phase 1 — you configure

Parameters, faults, severities, fault times. All local DOM state; faults queue
in `pendingFaults`. Nothing is validated yet because nothing has been sent.

### Phase 2 — press Run

`run()` clears the plots *first*, then opens a WebSocket to `/ws/run` and sends
`{action:"start", config}`. One socket per run — closed at the end, reopened
next time — so a run can never inherit state from the previous one.

### Phase 3 — the backend builds the plant (the first thing that can fail)

```python
try:
    sim = Simulation(msg.get("config", {}))
except ValueError as exc:
    await ws.send_json({"type": "config_error", "note": str(exc)})
```

`Simulation.__init__` calls `build_modifiers`, and **this is where every guard
fires**: the anchor check raises if a setpoint sits outside its cardinal window,
`_assert_monotone` raises on a Rosso pole, and the `absolute_modifier` /
`f_T(T_ss)` checks warn on a crippled setpoint.

**Why here:** these are physics validations and the physics lives in the
backend — the frontend cannot know that `pH_max = 5.0` is fatal. **Why the
`try`:** without it the raise escapes and kills the socket with no message, which
is the same silence the guard exists to end.

On success `{"type":"started"}` goes out and the frontend flushes
`pendingFaults` as individual `inject_fault` messages, each answered by a
`fault_ack`.

### Phase 4 — the frame loop

Per iteration:

**4a. `sim.step()` — one simulated minute.** In order: the faulted trajectory
(RK4, pH PI controller on a 2 s inner tick, actuator context from
`apply_faults`); the healthy counterfactual (same integration, empty fault list,
its own PI integral); the EKF predict, which advances the shadow through *the
identical code path* as the counterfactual — so "the EKF predicts healthy
dynamics" is true by construction rather than by a second copy of the ODEs that
could drift out of sync; then, only at the ethanol cadence, a noisy sample, the
correction, `z = ν/√S`, and `cusum_step` in the same block. The detector lives
inside the correction step exactly as it would on the rig: two floats of state,
no stored history.

**4b. Send the frame** — `true`, `healthy`, `ekf`, `cusum`, `nrmse_X`. The
frontend `extendTraces` rather than redrawing, so a 600-frame run stays smooth.

**4c. Fire the bank on the alarm**, and schedule two re-diagnoses at
`alarm + PROBE_SETTLE_MIN` and `alarm + ETHANOL_SETTLE_MIN`.

**Why automatic:** the verdict at the alarm has not settled (measured mean
35 min, max 90). The UI withholds the name until the window closes, so the run
must produce a fresh verdict at that moment by itself — leaving it to a click
asks the operator to know when the evidence matured, the one thing they cannot
know.

**4d. Diagnosis runs off the event loop** (`asyncio.to_thread`). The bank is 8
hypotheses × up to 5 onsets replayed over the whole batch; `to_thread` keeps the
frame stream flowing. A request arriving while one is in flight is dropped
rather than queued, since newer evidence supersedes it.

**4e. Drain the inbox.** A separate `receiver()` task queues incoming messages;
the loop drains it non-blockingly, so `inject_fault` / `diagnose` / `stop` are
honoured *between* frames, never mid-integration.

**4f. Pace to wall-clock** — sleeps the *remainder* of the frame budget, so
heavy frames do not accumulate lag.

### Phase 5 — rendering, per frame

`pushFrame` extends the traces; `setAlarmBadge` sets the pill to **ALARM · t =
3.50 h**, **BLIND** (tooltip naming whichever of `f_T` / `f_pH` is binding), or
**quiet · 0.31 / 1.5**.

### Phase 6 — the verdict arrives

`onMmae` → `renderRca`. If the window has not closed, **no cause is named** — a
progress bar, the reason, and a reveal button. Otherwise: verdict, tie set,
margin in nats, posterior ranking, per-hypothesis traces.

Both scorings come back from **one** bank run, so switching the ethanol/probes
tab re-renders from the same result — instant, and the two answers are always
over identical evidence.

### Phase 7 — end of run

An in-flight bank is awaited, then a **final diagnosis on the full batch**: a
late alarm may never have had one kicked off, and one kicked mid-run scored less
evidence than now exists. CSV and PNG are built server-side and sent inline with
`{"type":"complete"}`; `endRun()` closes the socket and enables the downloads.

### The decisions worth naming

* **One socket per run** — no cross-run contamination, and disconnect
  unambiguously means stop.
* **Backend owns all state** — the frontend holds only what it has been told, and
  the config is captured once at `start`, so the UI and the simulation cannot
  drift apart mid-run.
* **The healthy counterfactual runs in lockstep, not afterwards** — you watch the
  divergence open up live, and the EKF's prediction *is* that same computation.
* **Diagnosis is automatic at three moments** — the alarm and each settle
  boundary — because those are exactly the moments the answer changes.
* **Faults are injectable mid-run**, guarded so `t_fault` cannot already have
  passed, so you can break the reactor while watching it.
