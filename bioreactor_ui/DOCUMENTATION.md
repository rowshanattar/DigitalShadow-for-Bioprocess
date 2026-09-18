# Bioreactor Digital-Twin Simulator — Documentation

Developer and scientific reference for the local web app that wraps the
baker's-yeast EKF digital twin. This document covers the architecture, the model
equations, the parameter schema, the client/server protocol, and how to extend
the app.

---

## Table of contents

1. [Overview](#1-overview)
2. [Architecture](#2-architecture)
3. [Installation & running](#3-installation--running)
4. [Scientific model](#4-scientific-model)
5. [Parameter reference](#5-parameter-reference)
6. [Fault catalogue](#6-fault-catalogue)
7. [The EKF and the CUSUM detector](#7-the-ekf-and-the-cusum-detector)
7b. [Root cause: the MMAE bank](#7b-root-cause-the-mmae-bank)
   — including [the configuration guards & BLIND badge](#47-the-configuration-guards-and-the-blind-badge) (§4.7)
   and [the settle windows](#7b4-how-much-batch-the-verdict-is-standing-on) (§7b.4)
8. [HTTP & WebSocket API](#8-http--websocket-api)
9. [Data exports](#9-data-exports)
10. [Extending the app](#10-extending-the-app)
11. [Troubleshooting](#11-troubleshooting)
12. [Provenance & references](#12-provenance--references)

---

## 1. Overview

The app turns the research notebooks into an interactive tool:

- **Edit every parameter** (kinetics, initial state, thermal, pH, oxygen, EKF,
  run settings) with the notebook defaults preloaded.
- **Run in real time**: the six concentration/environment panels grow as the
  batch simulates, paced at **5 wall-clock seconds per simulated hour**.
- **Compare three trajectories**: the true (faulted) run, the EKF *healthy
  shadow* — what the reactor should be doing — and the pure no-fault
  counterfactual. Optional ±σ uncertainty band.
- **Watch the detector**: the ethanol innovation and the two CUSUM arms build up
  live against the threshold, with the alarm time latched in the EKF card.
- **Get a root cause**: the moment the alarm fires, an MMAE filter bank runs one
  EKF per candidate fault and returns a ranked posterior over causes — with the
  evidence behind it, scored on ethanol alone and again with the T/pH/DO probes.
  The panel **withholds the cause name until the evidence window has closed**
  and re-diagnoses itself when it does ([§7b.5](#7b5-when-it-runs-and-when-it-is-allowed-to-name-a-cause)).
- **Inject faults**, before or live during a run, across the three Isermann
  categories (Operation is active; Design and Ageing are placeholders).
- **Get told when the configuration is unusable**: a setpoint outside the
  organism's viable window is *rejected* rather than run silently, and a batch
  whose healthy reference makes almost no ethanol is flagged **BLIND** instead
  of reading as quiet ([§4.7](#47-the-configuration-guards-and-the-blind-badge)).
- **Download** the run as a standardized CSV and a multi-panel PNG.

Runs entirely on `localhost`; no cloud, no auth, single user.

---

## 2. Architecture

```
bioreactor_ui/
├── backend/                 # Python — physics, EKF, server (pure & parameter-driven)
│   ├── core.py              # Monod kinetics, Crabtree carbon split, Rosso/O₂ modifiers, DEFAULTS
│   ├── plant.py             # 8-state plant RHS + Van't Riet kLa + build_modifiers guards
│   ├── faults.py            # fault registry + catalogue + apply_faults()
│   ├── ekf8.py              # 8-state Design-A healthy-shadow EKF  ← default
│   ├── ekf.py               # 5-state Design-B EKF (symbolic Jacobian) — legacy
│   ├── cusum.py             # two-sided CUSUM detector on the innovation
│   ├── mmae.py              # MMAE filter bank: WHICH fault (onset ladder + Bayes)
│   ├── simulator.py         # stepwise streaming integrator (1 sim-min per frame)
│   └── main.py              # FastAPI: /api/defaults, WS /ws/run, CSV/PNG builders
├── frontend/                # Browser — no build step
│   ├── index.html           # config panels + chart grid + toolbar
│   ├── app.js               # WS client, Plotly live plots, fault UI, downloads
│   └── styles.css           # pastel, minimal theme
├── requirements.txt
├── README.md                # quickstart
└── DOCUMENTATION.md         # this file
```

**Data flow**

```
                    WS /ws/run  (JSON frames, paced 5 s / sim-hour)
  ┌───────────┐    ───────────────────────────────►    ┌──────────────┐
  │  app.js   │   start / inject_fault / stop  ◄─────   │ Simulation   │
  │  (Plotly) │                                          │  step()      │
  └───────────┘    ◄───────────────────────────────    │  + EKF.step  │
       ▲            frame / fault_ack / complete         └──────────────┘
       │                                                        │
   GET /api/defaults  (params, faults, cells, reactors)         │
                                                     on complete: csv + png
```

The backend is deliberately UI-agnostic: `core`/`plant`/`faults`/`ekf`/`simulator`
have no FastAPI imports and can be driven from a script or notebook. `main.py` is
the only web layer.

**Key design choices**

- **Frame vs. control tick.** The simulator emits one display frame per `dt_min`
  (default 1 min), but the PI pH controller runs on its own faster tick
  (`DT_CTRL_s`, default 2 s): each frame is integrated as several fixed RK4
  sub-steps, and the pump flows (and the error integral) are re-evaluated at each
  sub-step. This keeps the pH integration accurate while the display samples once
  per frame.
- **The EKF is a healthy shadow (Design A).** It integrates the *fault-free* ODEs
  — through the very same `_advance` call that produces the no-fault
  counterfactual — and is corrected only by the noisy ethanol sensor. So it
  answers "how should this reactor be behaving?", the innovation becomes the
  fault signal, and a CUSUM on it is the detector. The legacy Design B, where
  T/pH/DO enter as measured inputs through `φ` and the filter follows the fault,
  is still available via `ekf.old_design`. See
  [§7](#7-the-ekf-and-the-cusum-detector).
- **Parameter-driven, not hard-coded.** `DEFAULTS` in `core.py` is the single
  source of truth; the frontend builds its inputs from `/api/defaults`, so adding
  a parameter in one place surfaces it in the UI automatically.
- **Physics validation lives in the backend, at `start`.** The frontend cannot
  know that `pH_max = 5.0` is fatal, so `Simulation.__init__` → `build_modifiers`
  is where every guard fires, and `main.py` wraps the construction in a `try` so
  the refusal comes back as a `config_error` message rather than a dead socket
  ([§4.7](#47-the-configuration-guards-and-the-blind-badge)).
- **Static assets are served `Cache-Control: no-store`** (`NoCacheStatic` in
  `main.py`). This is a tool that gets edited while it runs, and a browser
  holding a stale `app.js` looks exactly like a broken backend.

---

## 3. Installation & running

**Requirements:** Python ≥ 3.10 (developed on 3.12). Dependencies are pinned by
lower bound in `requirements.txt`.

```bash
cd bioreactor_ui
python -m venv .venv && source .venv/bin/activate     # optional but recommended
pip install -r requirements.txt
uvicorn backend.main:app                               # --reload while developing
```

Open <http://localhost:8000>. Change the port with `--port 8123`.

> **Internet note:** Plotly is loaded from a CDN (`cdn.plot.ly`), so the first
> page load needs a connection. The simulation, EKF, faults, and exports all run
> locally. To go fully offline, download `plotly-2.35.2.min.js` into `frontend/`
> and change the `<script src>` in `index.html` to a relative path.

Stop the server with `Ctrl-C` (foreground) or `pkill -f "uvicorn backend.main:app"`.

---

## 4. Scientific model

### 4.1 Plant state vector

`plant.py` integrates an 8-state ODE system:

```
y = [ X,  G,  E,  µmax_G,  µmax_E,  T,  pH,  DO   ]
      g/L g/L g/L  1/h      1/h    °C   –   % sat
```

`µmax_G`, `µmax_E` are constant in the plant (`dµ/dt = 0`); they are augmented
*states* only in the EKF, where they are estimated online. `DO` (dissolved
oxygen, in % saturation) is a plant state fed to the EKF as a *measured input*
via `f_DO` (like T and pH; see [§4.5b](#45b-dissolved-oxygen-ode-free-floating-no-controller) and [§7](#7-the-ekf-and-the-cusum-detector)).

### 4.2 Growth kinetics (diauxic Monod)

With the growth modifier `φ_g = f_T(T) · f_pH(pH)` scaling the glucose uptake
**rate**, and oxygen acting separately on carbon **fate** (see §4.5c):

```
µ_G  = µmax_G · G/(K_G + G) · φ
supp = 1 − µ_G/µmax_G                    (diauxic suppression, linear form)
µ_E  = µmax_E · E/(K_E + E) · supp · φ
```

### 4.3 Mass balances

```
dX/dt = (µ_G + µ_E) · X
dG/dt = −(µ_G / Y_X/G) · X
dE/dt =  (µ_G · Y_E/G / Y_X/G − µ_E / Y_X/E) · X
```

### 4.4 Thermal ODE (P-controlled jacket, no active cooling)

```
Q_heat  = heater_frac · clip(K_p · (T_set − T), 0, Q_max)   # heater only
                                                 # heater_frac = 1 healthy,
                                                 # = Q_pct/100 under heater_off
Q_loss  = UA_eff · (T − T_amb)                   # Newton cooling
Q_metab = Y_Q/X · (µ_G + µ_E) · X · V            # metabolic heat
dT/dt   = (3600·(Q_heat − Q_loss) + Q_metab) / (V·ρCp)
```

Thermal time constant τ = V·ρCp / UA ≈ 5.2 min at the defaults.

### 4.5 pH ODE (acidogenesis vs. acid/base pumps)

The acid/base pumps are **throttled PI (proportional–integral) actuators**: the
pump command is the sum of a proportional term (∝ pH error) and an integral term
(∝ accumulated error), saturated at the maximum flow. The command is
re-evaluated on a short tick (`DT_CTRL_s`, default 2 s) and held over the tick:

```
err   = pH_set − pH                               # >0 → too acidic, add base
I    += err · dt                                  # error integral (anti-windup)
u     = K_c·err + K_i·I                            # signed pump command [mL/h]
F_B   = clip( u, 0, F_B)   (u ≥ 0 → base raises pH)
F_A   = clip(−u, 0, F_A)   (u < 0 → acid lowers pH)
yeast_acid = α_meta · (µ_G + µ_E) · X
dpH/dt = −yeast_acid + α_pump·F_B − α_pump·F_A     (clamped: no drop below pH 2)
```

The PI law produces a **smooth** trajectory (no relay sawtooth) and, unlike a
pure-P controller, has **no steady-state offset** — the integral drives pH back
to the setpoint. In a healthy run pH holds at 5.0. Under a **partially-stuck
acid pump** (`F_A_stuck < F_B_max`) pH **dips** at the moment of the fault (the
proportional term is still small), then the integral accumulates and **slowly
pushes pH back up to the setpoint** — a dip-and-recover, not a fight. Because pH
is restored to the growth-favourable range, biomass keeps rising through the
fault (its slope eases only briefly during the dip). Only when a stuck flow
exceeds the opposing pump's saturation flow does pH run away and growth collapse
below the cardinal minimum. `K_c` sets the dip depth; `K_i` sets the recovery
speed. Anti-windup freezes the integral while a pump is saturated. The plant is
stepped with a fixed RK4 sub-step at the control-tick rate (`simulator._rk4`).

> **Note on fidelity:** the real Minifors pumps are digital ON/OFF hardware; the
> PI law here is a deliberate simplification chosen for a clean, readable pH
> trace and a realistic dip-and-recover fault response. An earlier revision used
> a bang-bang relay whose full-flow pulses limit-cycled under a constant
> disturbance.

### 4.5b Dissolved-oxygen ODE (free-floating, no controller)

DO is expressed in **% saturation** (matching a Minifors pO₂ probe). It is a
free-floating balance of oxygen supply (aeration) against respiratory demand;
there is no DO controller. With `kLa_eff` the effective transfer coefficient for
the tick (reduced by an aeration fault):

```
supply  = kLa_eff · (DO_star − DO)                          # OTR   [%-sat/h]
demand  = qO2_max · (µ_G + µ_E)/(µmax_G + µmax_E) · X        # OUR   [mg/(L·h)]
dDO/dt  = supply − demand / DO_sat_mgL · 100                # %-sat/h (clamped ≥ 0)
```

`kLa` is not a fixed constant: `plant.compute_kLa(N_rpm, Q_air_Lmin, V)` derives
it from agitation and aeration via the **Van't Riet** correlation, so reducing
rpm or airflow (the `agitation_aer` fault) physically starves the supply term. The
shaft carries **two 6-blade Rushton turbines** (Minifors manual), modelled with
additive power (`_N_IMPELLERS = 2`); healthy kLa(500 rpm, 3.5 L/min) ≈ 44.8 1/h.
At the defaults DO sits near saturation (≈92–100 %), so `f_DO ≈ 1` and growth is
essentially unaffected in a healthy run; DO only throttles growth if it is driven
down toward the O₂ half-saturation `K_O`. The OUR here is a specific-uptake
heuristic (scaled by growth activity), not a per-substrate O₂ yield.

### 4.5c Oxygen sets carbon FATE, not growth rate (Crabtree)

Oxygen used to multiply into `φ` alongside `f_T` and `f_pH`, so losing oxygen
*reduced* ethanol. That is the wrong sign, and it is why this app and
`7_mmae_fault_isolation.ipynb` used to disagree about every DO fault.

Respiration has a hard capacity. Glucose arriving faster than `C_RESP · f_DO`
**overflows into fermentation**, so losing oxygen makes **more** ethanol, not
less — the Crabtree effect. Ethanol is itself a respiratory substrate, so without
oxygen it cannot be consumed back either.

```
q_G  = (µmax_G · G/(K_G+G) · φ_g) / Ygx     total glucose uptake
q_ox = min(q_G, C_RESP · φ_g · f_DO)         respired (capacity- and O2-limited)
q_fm = q_G − q_ox                            overflow → fermented
µ    = q_ox·YGX_OX + q_fm·YGX_RED + µ_E
dE   = (q_fm·YGE_FERM − µ_E/Yex) · X
```

**Note the `φ_g` on the cap.** Respiration is enzymatic, so cold or acid slows it
exactly as they slow uptake. An earlier version of this port left the cap at its
30 °C value, which made a *chilled* culture **more** respiratory — backwards, and
with a bad consequence: slow the culture ~3.5× and `q_G` fell below a fixed cap,
fermentation stopped altogether, and ethanol went silent. A heater fault costing
33 % of the biomass then produced **0.6 σ** of ethanol signal and could not be
detected at any threshold, where the pre-Crabtree model saw it at **19 σ**.
Scaling the cap keeps the fermentative *fraction* invariant to T and pH; oxygen
still sets carbon fate on its own through `f_DO`, which is the point of the
Crabtree rework and is untouched by this.

`YGE_FERM`, `YGX_OX`, `YGX_RED`, `C_RESP` are a **decomposition** of the inherited
YD yields, not a replacement: those yields already are a ⅓-respiratory,
⅔-fermentative blend, and this writes the blend down so the halves can move apart
when oxygen goes away. The healthy batch is essentially unchanged (final X 3.99 vs
3.98); only fault behaviour flips.

**A DO fault now has the opposite ethanol sign to a T or pH fault**, and that sign
is what lets the MMAE bank separate the oxygen loop from the others.

> ⚠️ **Inherited inconsistency, not fixed.** `kLa`, `qO2_max` and `DO_sat_mgL` are
> not stoichiometrically consistent: respiring `C_RESP` g glucose/gX/h needs ~274 mg
> O₂/gX/h, ~53× this model's peak OUR. `K_O2` keeps the inherited demand scale and
> only redirects it onto the oxidative flux. The practical consequence is that DO
> holds near saturation until `kLa` is almost gone, so the `agitation_aer_clog` rung
> sits at **1 %** of healthy `kLa` and is not a physically meaningful "restriction".
> Re-anchoring the oxygen balance against a measured OUR is an open item.

### 4.5d `acid_pump_off` is a null fault

Yeast acidogenesis only ever pushes pH **down**, so the controller only ever calls
for base: `F_A ≈ 0` in a healthy batch. Killing the acid pump therefore changes
essentially nothing — the ethanol channel moves by **0.0000 σ**. (Under the
notebook's pure-P law it is null to `0.0e+00`; the app's PI law overshoots to
pH 5.0017, so the acid pump fires briefly and the deviation is 7.7e-4 pH instead of
exactly zero.) A dead acid pump is undiscoverable until something *else* drives pH
up. It is kept in the catalogue so the MMAE screening can *name* it rather than
silently mis-attribute the alarm.


### 4.6 Re-anchored Rosso cardinal modifier (and the O₂ Monod modifier)

`make_cardinal_modifier(v_min, v_opt, v_max, anchor, n)` builds `f(v)` from Rosso
(1993), re-normalized so `f(anchor) = 1`. Reused for both temperature (`f_T`) and
pH (`f_pH`). Outside `(v_min, v_max)` it returns `PHI_FLOOR = 1e-6`, never a hard
zero — a hard zero hands the EKF a zero Jacobian and the filter can never correct
its way back.

**`n` selects the shape of the window, not the variable.** `n = 2` is the CTMI
temperature form; its denominator picks up a root *inside* `(v_min, v_max)`
whenever the optimum sits below the window midpoint, and the clip serves that
pole up as a clean 0/1 step carrying no gradient at all. `n = 1` is pole-free and
handles asymmetric arms. So:

- `f_pH` always takes **n = 1**. (On the old `(3.0, 5.0, 12.0)` triple, `n = 2`
  put a pole at ≈3.84 and made a stuck acid pump indistinguishable from an oxygen
  crash in the MMAE bank.)
- `f_T` takes **n = 2** for the usual right-skewed window and falls back to
  **n = 1** on a left-skewed one such as `(27, 30, 40)`. `build_modifiers` picks
  this automatically: `n_T = 2 if (T_opt − T_min) ≥ (T_max − T_opt) else 1`.

**Clipping happens after re-anchoring, and only from below.** Clipping to `[0, 1]`
first truncated the physically correct `f > 1` that arises when the reactor sits
nearer the optimum than the anchor does.

`absolute_modifier(v, v_min, v_opt, v_max, n)` is the companion: the *un-anchored*
value, relative to the organism's true optimum. `make_cardinal_modifier` normalises
at the setpoint, so `f(setpoint) ≡ 1` and it can only ever answer "how far from
where we hold it", never "is where we hold it any good". `absolute_modifier`
answers the second question, and it is what the guards and the BLIND badge read.

Oxygen uses a **different shape**: `make_do_modifier(K_O, DO_star)` builds a
saturating **Monod** term `f_DO(DO) = DO/(K_O + DO)`, re-anchored so
`f_DO(DO_star) = 1` (at saturation there is no growth penalty; only the deficit
below saturation damps growth). It is clipped to `[0, 1]`. `f_DO` does **not**
multiply into `φ_g` — it caps the respiratory branch instead (§4.5c), which is
also why both filters now build their Jacobians by forward differences
(see [§7](#7-the-ekf-and-the-cusum-detector)).

### 4.7 The configuration guards, and the BLIND badge

Three checks stand between a parameter set and a run. All three exist because a
**crippled healthy reference silently inverts the whole method**: the shadow EKF
and the CUSUM both measure the faulted reactor against "healthy", so if healthy
is already dead, a fault that moves a variable *back* toward the optimum makes
the culture grow **better** than the reference and no threshold reads that
correctly.

| # | check | where | on failure |
|---|---|---|---|
| 1 | the anchor lies inside the viable window | `core.make_cardinal_modifier` | **raises** `ValueError` |
| 2 | every modifier rises monotonically toward its optimum | `plant._assert_monotone` | **raises** `ValueError` |
| 3 | the setpoint itself is any good (`absolute_modifier ≥ 0.5`) | `plant.build_modifiers` | `RuntimeWarning` |

**1 — the anchor guard.** The modifier is normalised at the setpoint, so a
setpoint at or outside a cardinal bound normalises by ≈0. This used to be a
silent `a = 1.0` fallback. Observed with `pH_set = pH_max = 5.0`: healthy
`X_end = 2.50` (no growth at all), base-pump-off `X_end = 4.00`, and **no alarm**.

**2 — the monotonicity guard.** A Rosso denominator that changes sign inside the
window is invisible in the output, because the clip presents the pole as a tidy
step. Nothing but a shape check catches it; this is the regression test for the
`f_pH` step function described in §4.6, and it did not exist when that bug shipped.

**3 — the setpoint-quality warnings.** Two of them, one per channel. The thermal
one is checked at `plant.thermal_steady_state(params)`, **not** at `T_set`: a
proportional heater needs a standing error to make any power, so at equilibrium

```
T_ss = (K_p·T_set + UA·T_amb) / (K_p + UA) = (50·30 + 18·22)/68 = 27.88 °C
```

— the vessel runs 2.1 °C *below* setpoint all batch, and cardinals chosen against
`T_set` are 2 °C more generous than they look. A `T_min` of 27 looks like a 3 °C
margin and is really 0.9 °C.

**Why a warning and not a raise:** an off-optimum setpoint is a legitimate thing
to simulate. It is only misleading if unreported.

**The BLIND badge** is the same story surfaced live. `Simulation` tracks
`_hE_peak`, the peak ethanol the *healthy counterfactual* has reached; when that
stays under `3·σ_E_sensor` the detector is structurally blind — there is no
signal to lose, so a quiet CUSUM means "cannot see", not "nothing wrong". The
frame carries `blind`, `hE_peak`, `T_ss`, `fT_ss` and `fpH_set`, and the UI names
**whichever of f_T / f_pH is the smaller** as the binding constraint. (It used to
blame the cardinal T window unconditionally, which points the user at the wrong
knob whenever pH is the one that is off.) A culture growing that slowly respires
everything it takes up and ferments nothing — so this blinds **T and pH faults
only**. Oxygen faults still alarm, because they redirect carbon into fermentation
rather than slowing it.


---

## 5. Parameter reference

Defined in `backend/core.py :: DEFAULTS`, grouped exactly as the UI panels.
All values are user-editable; the table gives the shipped defaults.

### `kinetics`
| Key | Default | Unit | Meaning |
|---|---|---|---|
| `KG` | 0.1 | g/L | glucose half-saturation |
| `KE` | 0.1 | g/L | ethanol half-saturation |
| `Ygx` | 0.15 | – | biomass/glucose yield (code set) |
| `Yge` | 0.34 | – | ethanol/glucose yield (code set) |
| `Yex` | 0.43 | – | biomass/ethanol yield (code set) |
| `YGE_FERM` | 0.51 | – | ethanol/glucose on the **fermentative** branch (theoretical max) |
| `YGX_OX` | 0.25 | – | biomass/glucose **respired** |
| `YGX_RED` | 0.10 | – | biomass/glucose **fermented** |
| `C_RESP` | 0.2573 | g G/g X/h | respiratory capacity — the Crabtree cap (§4.5c). *Stale calibration*: fitted for a 1.968 g/L ethanol peak, which now sits near 2.29 |
| `YO_E_REL` | 2.0 | – | O₂ per g ethanol respired, relative to per g glucose |
| `K_O2` | 9.050 | mg O₂/g | O₂ per g oxidised substrate — see the stoichiometry caveat in §4.5c |

The last five are a **decomposition** of the first three, not extra freedom; see
[§4.5c](#45c-oxygen-sets-carbon-fate-not-growth-rate-crabtree).

### `initial`
| Key | Default | Unit |
|---|---|---|
| `X0` | 2.5 | g/L |
| `G0` | 5.0 | g/L |
| `E0` | 0.0 | g/L |
| `muG0` | 0.15 | 1/h |
| `muE0` | 0.08 | 1/h |

### `thermal`
| Key | Default | Unit | Note |
|---|---|---|---|
| `V` | 1.35 | L | working volume (Minifors 2) |
| `rhoCp` | 4180 | J/L/K | water |
| `Y_QX` | 12000 | J/g | metabolic heat yield |
| `T_set` | 30 | °C | setpoint |
| `T_amb` | 22 | °C | ambient — *confirm* |
| `UA` | 18.0 | W/K | U·A = 242 × 0.07 |
| `K_p` | 50 | W/K | heater P-gain — *confirm* |
| `Q_max` | 630 | W | heater max (Minifors manual) |
| `T_min/opt/max` | 5 / 30 / 40 | °C | Rosso cardinal T. The **n = 2** form is used while `T_opt − T_min ≥ T_max − T_opt`, otherwise n = 1 (§4.6). Judge these against `T_ss ≈ 27.9 °C`, not `T_set` (§4.7) |

### `ph`
| Key | Default | Unit | Note |
|---|---|---|---|
| `pH_set` | 5.0 | – | setpoint |
| `pH_min/opt/max` | 2.5 / 4.7 / 8.0 | – | cardinal pH, the **`robust`** organism profile of `7_mmae_fault_isolation.ipynb` (A.2). `pH_opt ≠ pH_set` deliberately. Always evaluated with Rosso **n = 1** (§4.6). *Confirm*: `ph_modifier_handoff.md` traces the magnitude to Hinga (2002), a review of **marine phytoplankton** — it supports the shape, not yeast values. The notebook's `sensitive` profile is (4.6, 5.0, 6.0) with cardinal T (20, 30, 40); override both groups together |
| `alpha_meta` | 0.4 | pH/h per g/L·h⁻¹ | acidogenesis (uncontrolled drop ≈ α·ΔX ≈ 0.5 pH) |
| `alpha_pump` | 0.2563 | pH/mL | buffer coefficient |
| `F_A`, `F_B` | 210 | mL/h | acid/base pump **max (saturation)** flow — also the stuck-ON severity, so a stuck pump exactly saturates its opposite number |
| `Kc_pump` | 120 | mL/h per pH | PI **proportional** gain — sets the depth of the pH dip when a fault hits |
| `Ki_pump` | 250 | mL/h per pH·h | PI **integral** gain — sets how fast pH recovers to setpoint after the dip |
| `DT_CTRL_s` | 2 | s | controller / integration tick (RK4 sub-step) |

### `oxygen`
| Key | Default | Unit | Note |
|---|---|---|---|
| `N_rpm` | 500 | rpm | agitation — Van't Riet input & agitation-fault target |
| `Q_air_Lmin` | 3.5 | L/min | aeration — Van't Riet input & aeration-fault target |
| `DO_star` | 100 | % sat | O₂ saturation reference (probe = 100 %) |
| `DO0` | 100 | % sat | initial DO (saturated at inoculation) |
| `qO2_max` | 8.0 | mg O₂/g/h | specific O₂ uptake — *confirm* |
| `K_O` | 3.0 | % sat | O₂ Monod half-saturation (≈0.2 mg/L) |
| `DO_sat_mgL` | 6.7 | mg/L | mg/L ↔ %-sat conversion constant |

Van't Riet geometry/correlation constants (`D_vessel_m`, `Pg_P0`, `N_power`, …)
live as module constants in `plant.py`, not in the UI. Several are flagged
`***CONFIRM***` pending calibrated values.

### `ekf`
| Key | Default | Meaning |
|---|---|---|
| `on` | true | run the EKF alongside the plant |
| `old_design` | false | use the legacy 5-state Design-B filter instead of the healthy shadow |
| `cadence_min` | 5.0 | ethanol measurement cadence [min] |
| `sigma_E` | 0.3162 | ethanol σ the **filter assumes** → `R = σ²` [g/L] |
| `sigma_E_sensor` | 0.05 | the **real** ethanol sensor noise σ [g/L] |
| `cusum_k` | 0.15 | CUSUM slack [units of √S] |
| `cusum_h` | 1.5 | CUSUM threshold [units of √S] |
| `band` | false | show the ±σ uncertainty band |
| `mmae_probes` | false | which root-cause tab opens first (both are always computed) |
| `sigma_T_probe` | 0.1 | Pt100 σ used by the MMAE probe scoring [°C] |
| `sigma_pH_probe` | 0.02 | glass-electrode σ [–] |
| `sigma_DO_probe` | 1.0 | pO₂-probe σ [% sat] |

The two σ's are separate on purpose: tying them together would make a degrading
sensor invisible, since the filter would simply trust it less. Raise
`sigma_E_sensor` alone to sweep sensor quality and watch the biomass-NRMSE badge
climb. `sigma_E = √0.1` is the R used throughout the thesis notebooks, and
`cusum_k`/`cusum_h` were derived against it — see
[§7](#7-the-ekf-and-the-cusum-detector).

In Design A (default) T, pH and DO are **estimated states** predicted from the
healthy ODEs; in Design B they are measured inputs entering through
`φ = f_T·f_pH·f_DO`.


### `sim`
| Key | Default | Meaning |
|---|---|---|
| `t_final` | 10.0 | batch length [h] |
| `sec_per_hour` | 5.0 | wall-clock seconds per simulated hour (pacing) |
| `dt_min` | 1.0 | integration/frame step [min] |

> Items marked *confirm* were flagged `*** CONFIRM ***` in the notebooks
> (pending values from the Minifors manual / literature). They are safe editable
> defaults, not final calibrated values.

---

## 6. Fault catalogue

Defined in `backend/faults.py`. A fault modifies the actuator context
(heater power, effective UA, pump flows, effective kLa) for `t ≥ t_fault`.

The registry is **four actuator loops × two failure modes**, matching
`7_mmae_fault_isolation.ipynb` (A.4):

| loop | dead | running-but-wrong |
|---|---|---|
| heater | `heater_off` `Q_pct = 0` | `heater_clog` `Q_pct = 6` |
| acid | `acid_pump_off` | `acid_pump_stuck` `F_A_stuck = 210` |
| base | `base_pump_off` | `base_pump_stuck` `F_B_stuck = 210` |
| oxygen | `agitation_aer_stop` `kla_pct = 0` | `agitation_aer_clog` `kla_pct = 1` |

### 6.1 What the *injector* offers (`FAULT_CATALOGUE`)

Six entries, not eight. For the heater and the oxygen loop, "dead" and
"running-but-wrong" are the **same code path** — both resolve to one continuous
severity (`Q_pct`, `kla_pct`) and 0 % is just the bottom rung — so offering them
as two catalogue entries with the same slider asked the user to pick a rung
twice. The pump loops keep two entries because dead (`F = 0`) and stuck-ON
(`F = saturation`) are genuinely different code paths.

| Category | Fault `id` | UI label | Effect | Params (default) |
|---|---|---|---|---|
| Operation | `heater_off` | Heater degraded | `Q_heat → (Q_pct/100)·Q_heat`; 0 % = dead (T decays to T_amb) | `t_fault` (0.5), `Q_pct` (0) |
| Operation | `acid_pump_off` | Acid pump OFF (dead) | `F_A = 0` — **a null fault**, see §4.5d | `t_fault` |
| Operation | `acid_pump_stuck` | Acid pump stuck ON | `F_A` frozen at `F_A_stuck` | `t_fault`, `F_A_stuck` (210 mL/h) |
| Operation | `base_pump_off` | Base pump OFF (dead) | `F_B = 0` → pH drifts down | `t_fault` |
| Operation | `base_pump_stuck` | Base pump stuck ON | `F_B` frozen at `F_B_stuck` | `t_fault`, `F_B_stuck` (210 mL/h) |
| Operation | `agitation_aer_stop` | Agitation / aeration degraded | `kLa_eff = kLa · kla_pct/100`; 0 % → DO crashes | `t_fault`, `kla_pct` (0) |
| Design | — | 🔜 | — | — |
| Ageing | — | 🔜 | — | — |

`heater_clog` and `agitation_aer_clog` **stay in `FAULT_DEFAULTS` and stay
resolvable in `apply_faults()`** — the MMAE bank carries all eight, because a
hypothesis must commit to a severity before it can be integrated. They are simply
not offered as separate things to inject. See
[§7b.2](#7b2-hypothesis-severity--now-the-same-numbers-the-injector-uses).

### 6.2 Why the derate rungs are not round numbers

- **`heater_clog` `Q_pct = 6`.** The P-controller only draws `UA·ΔT = 18·8 = 144 W`
  of the 630 W available at steady state, so `Q_max` has to fall below ~23 %
  before the heater loses authority at all. `Q_pct = 15` is invisible (0.7 σ);
  6 % → 37.8 W.
- **`agitation_aer_clog` `kla_pct = 1`.** `kla_pct = 15` is invisible (0.4 σ): DO
  only falls to 73 % sat, where `f_DO = 0.989`. That is the inherited
  oxygen-stoichiometry caveat in §4.5c — demand is ~53× too small — so the DO
  ladder stays squeezed into `kla_pct ≤ 1` until the O₂ balance is re-anchored on
  a measured OUR. A "clog" needing a 99 % cut in kLa is **not** physically
  meaningful, and is labelled as such rather than presented as a mild fault.

### 6.3 Two structural facts about this registry

- **Agitation and aeration are one mechanism, not two.** Both reach the biology
  only through kLa, so any `(rpm, air-flow)` pair giving the same kLa produces a
  bit-identical trajectory in all 8 states — the notebook proved
  `agitation_stop == aeration_stop` to 1e-9. They are one hypothesis
  parameterised two ways, and only an instrument on the actuator itself (a
  tachometer, an air mass-flow meter) can separate them. Severity is therefore a
  **kLa ratio**, not an rpm.
- **`vessel_leak` was removed.** It is a passive thermal fault, and at the severe
  rung it lands within 0.2 % of the same biomass loss as `heater_off` — so it
  added a permanent tie to the bank without adding a mechanism the rig can act on.

**Scheduling semantics** (`Simulation.inject_fault`): a fault is accepted only if
`t_fault ≥` the current sim time (with a small epsilon). Injecting a past time
returns `(False, message)`, surfaced to the user as an error toast. Multiple
faults can be active simultaneously; `apply_faults()` folds them all in each tick.

---

## 7. The EKF and the CUSUM detector

Two estimator architectures ship side by side. **Design A is the default**;
Design B is kept behind the `old design (B)` switch so the two can be compared
on the same fault.

### 7.1 Design A — the 8-state healthy shadow (`backend/ekf8.py`)

- **State (8):** `[X, G, E, µmax_G, µmax_E, T, pH, DO]`. T, pH and DO are
  **estimated states here**, integrated from the *healthy* ODEs. Nothing about
  the faulted reactor is fed in — the noisy ethanol sample is the filter's only
  contact with reality.
- **Predict:** the simulator advances the filter's state through
  `Simulation._advance(..., faults=[], ...)` — the exact call that produces the
  no-fault counterfactual trace. "The EKF predicts with healthy dynamics" is
  therefore true by construction, not by a second copy of the ODEs. Covariance
  is propagated in `EKF8.predict_cov` as `P ← Φ P Φᵀ + Q·dt`, `Φ = I + F·dt`,
  with `F` the forward-difference Jacobian of the healthy RHS.
- **Correct:** Joseph-form, `H = e₃ᵀ` (ethanol only), plus physical clamps —
  non-negativity, `pH ≥ 1`, and the **batch mass balance**: `X` may not fall and
  `G` may not rise, because the plant has no death term and no feed. Without
  those the ethanol-only correction drags the estimate below the faulted truth
  and *creates* glucose.
- **Matrices** (ported from `fault_detection_CU.ipynb` Cell 4):
  `P0 = diag(0.1, 0.02, 0.02, 1e-5, 1e-5, 1e-4, 1e-6, 1e-4)`,
  `Q = diag(1e-3, 1e-3, 1e-3, 1e-9, 1e-9, 1e-4, 1e-6, 1e-4)`, `R = σ_E²`.
  µmax is pinned (1e-5/1e-9) on purpose: the filter must not be able to absorb a
  fault by retuning growth.

**Consequence (the architectural result):** the filter is a reference for *how
the reactor should be behaving*. Every run therefore shows three traces —

| trace | meaning |
|---|---|
| **true** | what is really happening (faulted plant) |
| **EKF** | what we expect to happen (healthy shadow, pulled by ethanol) |
| **no fault** | what should have happened if nothing were wrong |

and under a fault the EKF sits **between** the other two: healthy dynamics push
it toward the counterfactual, the ethanol correction pulls it toward truth. The
gap is the fault signal.

The shadow's T/pH/DO overlays sit almost on the healthy curve. That is correct:
ethanol cannot observe them, so they move only through the corrected X/G/E
feeding back into metabolic heat and acidogenesis.

### 7.2 The CUSUM detector (`backend/cusum.py`)

The innovation `ν = z_E − Ê` is normalised by the filter's own predicted spread,
`z = ν/√S` with `S = H·P·Hᵀ + R`, so a fixed threshold means the same thing at
every step. Two arms are then accumulated once per correction:

```
S⁺ ← max(0, S⁺ + z − k)      ethanol ABOVE prediction
S⁻ ← max(0, S⁻ − z − k)      ethanol BELOW prediction  ← the fault direction
```

An alarm fires the first time either arm exceeds `h`, and **latches** — the arms
are never reset, so they stay a monotone evidence trace and the UI reports when
the fault was *first* seen. The UI shows both arms against `h` in the right-hand
detector panel, the raw ν against the ±0.1 g/L "cultivation interrupted" band in
the left one, and the alarm time as a red pill in the EKF card.

`k = 0.15` and `h = 1.5` are both in units of √S (√S = 0.329 g/L at R = 0.1), and
they have **different provenance**:

- **`k` is derived.** `k = z₁/2` with `z₁ = 0.1 g/L / √S = 0.304`, the ethanol
  shortfall at which the cultivation counts as interrupted. Classical δ/2.
- **`h` is not.** The claim that h = 1.5 came from an MC over 20k healthy batches
  was audited against the notebook that ran it (`fault_detection_CU.ipynb`): the
  executable MC uses **4000** batches and at k = 0.15 returns **h = 0.93** for
  0.1 % per batch. 1.5 is ≈1.6× that — a conservative hand-picked threshold the
  MC *supports*, not one it produced. Measured false-alarm rate at 1.5: none in
  400k synthetic batches.

**The portable reading is in σ units of the channel's own healthy spread**, which
was measured twice independently (0.155 over 30 healthy batches in
`fault_detection_CU.ipynb`, 0.1506 over 6 in `8_evaluation.ipynb`):
**`k ≈ 1 σ_z`, `h ≈ 10 σ_z`**. Those ratios transfer to a new channel; the literal
0.15 / 1.5 do not, and copying them onto a correctly-normalised residual
false-alarms every batch. They also depend on R, which is why `sigma_E` defaults
to `√0.1 = 0.3162` — the R used throughout the thesis notebooks
(`ekf_bioreactor.ipynb` Cell 4). Change `sigma_E` and k/h must be re-derived.

Two caveats when comparing UI alarm times to the notebook:

- **The UI ships the `robust` organism profile, the notebook runs `sensitive`.**
  The shipped defaults are cardinal T (5, 30, 40) with cardinal pH
  (2.5, 4.7, 8.0); `7_mmae_fault_isolation.ipynb`'s `sensitive` profile is
  cardinal T (20, 30, 40) with cardinal pH (4.6, 5.0, 6.0) — **override both
  groups together**, they are one organism. With the literature `T_min = 5` a
  dead heater still leaves a large `f_T` and costs ~0 % biomass, so thermal
  faults are undetectable from ethanol at *any* threshold. Tightening the
  cardinal T window is what makes them visible — but tighten it against
  `T_ss ≈ 27.9 °C`, not `T_set`, and watch the BLIND badge: **both arms matter**,
  (27, 30, 32) gives `f_T = 0.28` and never alarms, (27, 30, 33) gives 0.50 and
  alarms at 1.17 h ([§4.7](#47-the-configuration-guards-and-the-blind-badge)).
- `h` was tuned over an 84-sample batch (7 h at 5 min). The UI's default batch is
  10 h ≈ 120 samples, so the false-alarm rate is slightly above 0.1 %.

### 7.3 Design B — the legacy 5-state `EKF` (`backend/ekf.py`)

Enabled with `ekf.old_design`. The detector row is hidden while it runs: its
innovation is not a fault signal.

- **State (5):** `[X, G, E, µmax_G, µmax_E]`. T, pH and **DO** are **not**
  estimated — they enter the predict step as measured inputs via
  `φ_g = f_T(T)·f_pH(pH)` plus `f_DO(DO)` as a separate cap on respiration. The
  Jacobian is now built by **forward differences**, because the carbon split
  contains `min(q_G, C_RESP·f_DO)` and is not differentiable at the
  respiro-fermentative switch — exactly where a DO fault lives. The old note here
  said folding `f_DO` into a single symbolic `φ` required no Jacobian rebuild; that
  was true and it was the problem — the DO change was
  one line in `simulator.step`.
- **Jacobian:** built once symbolically with SymPy (φ carried as a symbol),
  `lambdify`-ed, and evaluated per step. Parameters stay symbolic so UI edits take
  effect without rebuilding.
- **Predict:** integrates the combined 30-vector `[x; vec(P)]` over one frame with
  `dP/dt = F·P + P·Fᵀ + Q` (RK45, φ fixed over the tick).
- **Correct:** ethanol-only observation `H = [0,0,1,0,0]`, gain `K = P·Hᵀ·S⁻¹`,
  `S = H·P·Hᵀ + R`, `R = max(σ_E², 1e-4)`. Covariance uses the **Joseph form**
  for numerical stability. A noisy ethanol sample is drawn every `cadence_min`.
- Default noise/covariance (in `EKF.__init__`): `P0 = diag(0.1, 0.02, 0.02, 0.2,
  0.02)`, `Q = diag(1e-3 ×5)`.

**Consequence:** because the filter compensates for *known* T/pH/DO changes
through φ, it **follows** every actuator fault — a heater fault or an
`agitation_aer_stop` (kLa) fault leaves the ethanol innovation near zero. It is a
state estimator for the real reactor, which is useful, but it leaves nothing for
a detector to sit on. That is precisely why Design A exists.

---

## 7b. Root cause: the MMAE bank

The CUSUM says **that** something is wrong and roughly **when**. It says nothing
about **what**. `backend/mmae.py` adds the what, ported from
`7_mmae_fault_isolation.ipynb` as the same five steps:

| # | step | in code |
|---|---|---|
| 1 | the shadow alarms on the ethanol innovation | `simulator.step` |
| 2 | walk the firing CUSUM arm back to where it left zero → `t̂₀` | `mmae.onset_estimate` |
| 3 | one EKF per candidate cause, each *believing* its own fault started at `t̂₀` | `mmae._Hyp`, `_replay` |
| 4 | feed all of them the **same recorded** ethanol samples | `Simulation.meas_rows` |
| 5 | Bayes over their innovation sequences | `mmae._score` |

A hypothesis filter is **not** a second EKF implementation: it is the healthy
shadow with a non-empty fault list handed to the very same
`Simulation._advance` / `Simulation._rhs_for`. Before its own onset a hypothesis
*is* the healthy shadow, which is why the bank forks all hypotheses from one
shared healthy prefix (exactly, not approximately) instead of replaying each
from `t = 0`.

### 7b.1 The two rules that make it trustworthy

Both are measured in the notebook, not assumed:

- **Onset ladder, backwards only** (`LADDER_H = 0, −0.5, −1.0, −1.5, −2.0 h` —
  five rungs, out to two hours). The CUSUM onset is the weak link and its error
  is asymmetric in **both** directions, which is why the ladder exists at all: a
  raised arm means nothing on its own (on healthy data an arm leaves zero ~11×
  per batch on noise alone), so `onset_estimate`'s walk back can run through a
  pre-fault excursion and land early, while a slow fault whose early innovations
  sit below the slack `k` lands late. Only the *late* error is fatal — a late
  onset always favours whichever cause acts fastest — so the ladder only ever
  searches backwards.
- **Profile, do not marginalise.** The onset axis is collapsed with `max`, never
  `sum`. Summing rewards *slow* faults, which degrade gracefully as `t₀` moves
  and so accumulate mass across many onsets whatever is true; it produces
  confident wrong answers. Scoring also uses the **whole batch**, because once
  candidates carry different onsets a per-onset window would reward the latest
  candidate for having fewer samples to explain.

An **honest tie beats a confident error**, so the verdict reports every cause
within `TIE_MARGIN_NATS = 2.3` (10:1 odds) of the winner rather than the bare
argmax. Causes with numerically identical scores are reported as **provably
inseparable on this batch** — a degeneracy of the plant, not a shortfall of the
estimator.

### 7b.2 Hypothesis severity — now the same numbers the injector uses

A hypothesis has to commit to a severity: the bank asks *which mechanism*, not
*which mechanism, how hard*. `HYPOTHESIS_SEVERITY` is that commitment, and it
**now matches `FAULT_DEFAULTS`**. This used to be a deliberate disagreement — the
injection defaults sat below the detection rung (a 105 mL/h stuck pump is
recovered by the PI controller), so the bank had to override them. The registry
rebuild as four loops × two failure modes (§6) fixed that at the source: the
severity *is* the mechanism now, so the bank and the injector no longer need to
disagree.

**The bank is still wider than the injector — eight hypotheses to six catalogue
entries.** The injector collapses each continuously-parameterised loop into one
entry with a % field; the bank cannot, because every hypothesis must commit to a
number before it can be integrated. Dropping `heater_clog` / `agitation_aer_clog`
would not simplify the diagnosis, it would delete the answer *"the heater is weak
but not dead"* from the set of answers the panel can give. `HYPOTHESIS_LABELS`
spells the labels out locally for the same reason: `heater_clog` no longer has a
catalogue entry to borrow one from.

`candidate_causes()` gates on the **category** (an "Operation" bank exists; Design
and Ageing do not yet), never per fault id.

The cost of committing to a severity: a fault at a *different* severity is still
matched to the right mechanism, with a smaller margin — but a mechanism the bank
cannot express at any severity may lose to the null. When the null wins, the UI
says so explicitly — "the bank explains the alarm best with no fault at all" —
rather than naming a runner-up.

### 7b.3 Ethanol only vs. the probes — both, from one bank run

The scoring is the only thing that differs between the two answers, so both are
computed in a single pass and returned together (`scores.ethanol`,
`scores.probes`); the UI's mode tabs just re-render.

- **`ethanol`** — the notebook's question: what can the *soft sensor* tell you?
  At the severities where these faults are detectable at all, the answer is
  usually a **4-way tie**. That is the notebook's saturation result, reproduced:
  heater failure, vessel leak, a jammed acid pump and a dead sparger all end in
  "growth stalled", and ethanol cannot say which loop stalled it. The bank still
  rejects the null decisively and flags the causes that are *silent on ethanol*
  (never separating from "no fault", so unisolatable at any threshold).
- **`probes`** — the same bank with the rig's own T/pH/DO instruments in the
  likelihood (σ from `ekf.sigma_T_probe` / `sigma_pH_probe` / `sigma_DO_probe`).
  This breaks the tie, because that is exactly where the causes differ: a dead
  heater and a jammed acid pump produce the same ethanol curve at 22 °C vs
  28 °C. Verdicts become decisive on every implemented fault.

`ekf.mmae_probes` only chooses which tab opens first. Log-odds under the probes
run to 10⁴–10⁷ nats: the plant model is deterministic, so a matching mechanism
beats a non-matching one by more than any real instrument justifies — read the
*size*, not the digits, which is why the UI prints `≫ 10³` and per-row Δ rather
than raw totals.

### 7b.4 How much batch the verdict is standing on

The bank fires the instant the CUSUM latches, so by default it scores the
`[0, alarm]` window. The alarm fires as soon as the innovation is *improbable*,
which is well before the drift has a recognisable **shape** — so the ethanol
verdict at that moment has not settled, while the probes, which separate causes
by *which channel moved and in which direction*, need barely any extra evidence.
That asymmetry is what the two constants encode:

```python
ETHANOL_SETTLE_MIN = 120.0     # mmae.py
PROBE_SETTLE_MIN   =  30.0
```

#### Two questions, not one: the loop and the rung

The registry is four actuator **loops** × two severity **rungs** (§6). The bank
is therefore answering two questions at once, and they have very different
answers — so accuracy has to be reported against both:

- **loop** — *which control loop broke*. The actionable half: it says which
  cabinet to open.
- **rung** — *how dead that loop is*. `heater_off` vs `heater_clog`.

#### What the constants are actually worth

Measured over **30 trials** (6 alarming scenarios × 5 sensor-noise seeds, fault
at t = 2 h, robust organism) — `evaluation_accuracy_windows.csv`, reproduced by
`evaluation_accuracy_windows.py`. Accuracy is **pointwise**: the fraction of runs
whose verdict is correct if you read the panel at that offset.

| offset after alarm | ethanol · rung | ethanol · loop | probes · rung | probes · loop |
|---|---|---|---|---|
| 0 (alarm) | 40.0 % | 60.0 % | 83.3 % | **100 %** |
| 30 | 60.0 % | 76.7 % | 83.3 % | 100 % |
| 60 | 76.7 % | 76.7 % | 83.3 % | 100 % |
| 90 | 83.3 % | 83.3 % | 83.3 % | 100 % |
| **120** ← release | **83.3 %** | 83.3 % | 83.3 % | 100 % |
| 150 | 96.7 % | 96.7 % | 83.3 % | 100 % |
| 240 | 100 % | 100 % | 83.3 % | 100 % |

> ⚠️ **`ETHANOL_SETTLE_MIN = 120` releases the name at ~83 % rung accuracy, not
> at certainty.** One ethanol diagnosis in six is still wrong when the UI stops
> calling it provisional. 95 % is not reached until **150 min**, and 100 % not
> until 240. The constant is kept at 120 deliberately — the panel always shows
> the tie set alongside the name, so a wrong argmax inside a reported tie is a
> shortlist, not a false diagnosis — but it must not be described as a
> settled-verdict guarantee.
>
> An earlier single-seed measurement (`8_evaluation.ipynb`, mean 35 min / max
> 90) suggested 120 was ~3.4× conservative. That seed was lucky. Across 5 seeds
> the straggler is `base_pump_stuck`, stuck at 20 % correct from the alarm all
> the way through 120 min and only reaching 80 % at 150; `heater_off` is second
> at 80 % through 120.
>
> The curve is also **not monotone** — 83.3 % at 90, 80.0 % at 105, 83.3 % at
> 120. Verdicts flip near the margin, which is direct evidence that a fixed
> clock is a crude instrument for this. Gating on *verdict stability* (the argmax
> unchanged across the last N re-scores) would track the real thing; it is an
> open item, not implemented.

**The probe rung accuracy has a ceiling below 95 %, at every window, forever.**
The two oxygen rungs mutually confuse: `agitation_aer_stop` is named correctly in
2 of 5 seeds and `agitation_aer_clog` in 3 of 5, because a DO probe reads the same
floor under kLa = 0 % and kLa = 1 % and is reduced to guessing between them. So
for the probes there is no "time to 95 % on the rung" — only a ceiling. On the
**loop**, which is the actionable question, they are at 100 % from the alarm
itself. `PROBE_SETTLE_MIN = 30` is therefore conservative for the loop question
and irrelevant to the rung one.

`run_mmae` stamps each scoring with `provisional` and
`evidence_min_after_alarm`, and returns both window lengths in the payload.

**The UI does not name a cause while `provisional` is set.** It shows a progress
bar (`forming · 45 / 120 min`), the reason, and a *Show the provisional answer
anyway* button — nothing is hidden, only unshown. Naming a cause on screen and
quietly replacing it 40 minutes later reads as the method being unreliable, when
what it actually is, is honest about how much batch it has seen.

### 7b.5 When it runs (and when it is allowed to name a cause)

Diagnosis is automatic at **three** moments, in a worker thread
(`asyncio.to_thread`) so the run keeps streaming — a bank of 8 hypotheses × up to
5 onsets takes 2–5 s at the default batch:

1. **the alarm latches** — cheap, and the panel needs something to show;
2. **each settle window closes.** On latching, `main.py` schedules
   `redo_at = [alarm + PROBE_SETTLE_MIN, alarm + ETHANOL_SETTLE_MIN]` (dropping
   any past `t_final`) and re-kicks the bank as the sim clock passes each. This
   has to be automatic: since the UI withholds the name until the window matures,
   leaving the re-score to a click would ask the operator to know *when* the
   evidence matured — the one thing they cannot know.
3. **end of batch** — a run that alarmed late may never have had a diagnosis
   kicked off, and one kicked mid-run scored less evidence than now exists.

Plus on demand from the **Re-diagnose** button. Only one diagnosis runs at a
time; a request arriving while one is in flight is **dropped, not queued**
(the client gets `mmae_status: "busy"`), since newer evidence supersedes it. The
bank reads a snapshot of `Simulation.rows` / `meas_rows`; `_advance` touches no
mutable simulation state, so it is safe alongside a live run.

**Known limitation.** The per-sample likelihood carries the notebook's
`−½·log(2πS)` normalisation, and `S` differs slightly between hypotheses. At the
sub-nat margins of an ethanol-only tie, that term can reorder the top of the
ranking — which is precisely why the tie is reported as a tie.

---

## 8. HTTP & WebSocket API

### `GET /api/defaults`
Returns everything the frontend needs to render:
```json
{
  "params":   { "kinetics": {...}, "initial": {...}, "thermal": {...},
                "ph": {...}, "oxygen": {...}, "ekf": {...}, "sim": {...} },
  "faults":   [ { "category": "Operation", "enabled": true, "faults": [...] }, ... ],
  "cells":    [ { "id": "yeast", "label": "...", "enabled": true }, ... ],
  "reactors": [ { "id": "minifors2", "label": "...", "enabled": true }, ... ]
}
```

### `WS /ws/run`
JSON messages, both directions.

**Client → server**
| `action` | Payload | Meaning |
|---|---|---|
| `start` | `{ config: {...} }` | begin a run; `config` mirrors the `params` schema (partial configs are merged over defaults) |
| `inject_fault` | `{ fault: { fault: "<id>", t_fault, ...params } }` | schedule/inject a fault |
| `diagnose` | – | re-run the MMAE bank on all evidence so far |
| `stop` | – | end the run early |

**Server → client**
| `type` | Payload | Meaning |
|---|---|---|
| `started` | `{ t_final }` | run accepted |
| `config_error` | `{ note }` | **the plant refused to build** — a setpoint outside the organism's viable window, or a Rosso pole (§4.7). The run never starts; the UI raises an error toast and calls `endRun()` |
| `frame` | `{ t_h, true:{X,G,E,T,pH,DO,muG,muE}, healthy:{X,G,E,T,pH,DO}, ekf:{X,G,E,muG,muE,sX,sG,sE(,T,pH,DO)}\|null, cusum:{...}\|null, nrmse_X, faults:[...] }` | one tick (`DO` in % sat). `ekf` carries `T/pH/DO` and `cusum` is non-null only in Design A. |
| `fault_ack` | `{ ok, note }` | result of an `inject_fault` |
| `mmae_status` | `{ state: "running"\|"busy", t_evidence }` | the bank started (or was already running) |
| `mmae` | see below | a finished root-cause diagnosis |
| `complete` | `{ stopped, csv, png }` | run finished; `csv` is text, `png` is base64 |

**`cusum` payload** (Design A only):

```json
{ "nu": -0.42, "z": -1.28, "Sp": 0.0, "Sn": 1.83, "h": 1.5, "alarm_t": 3.50,
  "blind": false, "hE_peak": 2.29,
  "T_ss": 27.88, "fT_ss": 0.98, "fpH_set": 0.96 }
```

`blind` is `hE_peak < 3·σ_E_sensor` — the healthy reference makes less ethanol
than the sensor noise, so a quiet CUSUM means "cannot see". `fT_ss` is `f_T` at
the temperature the P-controller **actually reaches**, and `fpH_set` is the
**un-anchored** `absolute_modifier` at `pH_set` (`f_pH(pH_set)` is 1 by
construction and would say nothing). The UI names the smaller of the two as the
binding constraint — see [§4.7](#47-the-configuration-guards-and-the-blind-badge).

**`mmae` payload.** `state` is `ok`, `no_alarm` (nothing to isolate yet) or
`unavailable` (Design B — no healthy shadow). When `ok`:

```json
{
  "state": "ok", "alarm_t": 3.50, "t_onset": 3.17,
  "onsets": [3.15, 2.65, 2.15, 1.65, 1.15],   // the 5-rung ladder, snapped to frame bounds
  "t_evidence": 7.0, "n_samples": 84,
  "ethanol_settle_min": 120.0,           // how much batch after the alarm each
  "probe_settle_min": 30.0,              //   scoring needs before it settles
  "default_mode": "ethanol",             // which tab the UI opens on
  "silent": ["base_pump_stuck", "base_pump_off"],
  "scores": {
    "ethanol": { "verdict": "...", "verdict_label": "...", "tied": [...],
                 "margin_nats": 0.1, "degenerate": [["a","b"]], "order": [...],
                 "p": {id: 0.29}, "loglik": {id: 15.4},
                 "post": {id: [...]},    // running posterior, per sample
                 "onset": {id: 3.15},    // the profiled-best onset per cause
                 "provisional": true,    // → the UI withholds the cause NAME
                 "evidence_min_after_alarm": 45.0 },
    "probes":  { ... same shape ... }
  },
  "t_meas": [...], "z_meas": [...],      // the sample times and sensor values
  "meas": { "T": [...], "pH": [...], "DO": [...] },   // what the probes read
  "causes": [ { "id": "...", "label": "...", "sep": 2.2, "silent": false,
                "nu": [...], "E": [...], "T": [...], "pH": [...], "DO": [...] } ]
}
```

Per-cause `nu`/`E`/`T`/`pH`/`DO` are that hypothesis' trace at its
ethanol-profiled onset — the evidence the UI plots.

The server paces frames so one simulated hour takes `sim.sec_per_hour` wall-clock
seconds (set `sec_per_hour: 0` to run as fast as possible — used in tests).

**Minimal client example** (Python):
```python
import asyncio, json, websockets
async def main():
    async with websockets.connect("ws://localhost:8000/ws/run") as ws:
        await ws.send(json.dumps({"action": "start",
            "config": {"sim": {"t_final": 10, "sec_per_hour": 0}, "ekf": {"on": True}}}))
        async for raw in ws:
            m = json.loads(raw)
            if m["type"] == "complete":
                open("run.csv", "w").write(m["csv"]); break
asyncio.run(main())
```

---

## 9. Data exports

Delivered inside the `complete` WebSocket message (not separate endpoints), so
the client can offer instant downloads without extra round-trips.

- **CSV** — columns `time_h, T, pH, DO, G, E, X, mu_G, mu_E` (the notebook
  convention). One row per frame. `DO` is the live dissolved-oxygen state in
  **% saturation**. With the EKF on, `ekf_X, ekf_G, ekf_E` are appended; in
  Design A, `nu, Sp, Sn` follow.
- **PNG** — a matplotlib figure, `nrow × 3` (Biomass, Ethanol, Glucose,
  Temperature, pH, DO), with a third row in Design A holding the innovation and
  the CUSUM arms, and a fourth once a diagnosis exists: the posterior over causes under **both**
  scorings, the per-hypothesis innovation, and each hypothesis' predicted
  temperature against what the probe read.
  EKF estimates appear as dashed overlays on every panel that has one — all six
  in Design A, only X/E/G in Design B, where T/pH/DO are measured inputs.
  Scheduled faults are marked with vertical dotted lines and the alarm with a red
  one. Built with the `Agg` backend (headless-safe).

Builders live in `backend/main.py` (`build_csv`, `build_png`) and operate on
`Simulation.rows`, `Simulation.ekf_rows` and `Simulation.det_rows`.

---

## 10. Extending the app

### Add a parameter
Add the key to the relevant group in `DEFAULTS` (`core.py`). It appears in the UI
automatically. Optionally add a nicer label/unit to `LABELS` in `app.js`.

### Add an operational fault
1. Add defaults to `FAULT_DEFAULTS` and an `apply_faults()` branch in `faults.py`.
2. Add a catalogue entry (with `params`) under the Operation category; set
   `enabled: true`. The UI renders its inputs and inject button. Skip this step
   if the new rung is just another value of an existing loop's % field (§6.1).
3. Add it to `HYPOTHESIS_SEVERITY` **and** `HYPOTHESIS_LABELS` in `mmae.py`, or
   the bank cannot name it. The bank is keyed on those dicts, not on the
   catalogue, so a fault with no catalogue entry is still diagnosable — and one
   with no `HYPOTHESIS_SEVERITY` entry is injectable but invisible to the
   diagnosis.

### Add a measured channel (DO is the worked example — already done)
DO was added this way; use it as the template for a future channel (e.g. CO₂):
1. Promote the channel to a plant state in `plant.py` (add its ODE + a `f_*`
   modifier in `core.py`); add its params as a new `DEFAULTS` group.
2. Fold a *rate* modifier into `φ_g` at the **caller** — `plant.make_rhs`
   (plant growth) and `simulator.step` (the value handed to `ekf.predict`).
   Both EKFs use finite-difference Jacobians, so nothing needs rebuilding; a
   *measured-input* channel needs **no Jacobian rebuild**.
3. In `simulator.step`, emit the real value in the frame and CSV rows.
4. In `app.js`, add the plot id to the `ENV` map (and a `grp-*` config panel to
   `index.html`); in `main.build_png`, add its panel.

### Add a cell line / reactor
Append to `CELLS` / `REACTORS` in `main.py` and provide a matching parameter
preset. Set `enabled: true` when ready.

### Style
All colors are CSS variables at the top of `styles.css` (`--accent`, `--green`,
`--peach`, `--lilac`, …). Plot colors live in the `CONC`/`ENV` maps in `app.js`.

---

## 11. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| Run refused with an error toast about a setpoint | The anchor guard (§4.7): `pH_set`/`T_set` sits at or outside a cardinal bound, so the *healthy* reference would be dead and the detector would read faults backwards. Move the setpoint inside the window, or widen the window. |
| Run refused with "a Rosso denominator has changed sign" | `_assert_monotone` caught a pole hidden by the clip (§4.6). Check the cardinal triple — an optimum far below the window midpoint under n = 2. |
| Badge reads **BLIND** instead of quiet/alarm | The healthy counterfactual makes less ethanol than the sensor noise, so nothing can be detected on T or pH faults. Hover the badge: it names whichever of `f_T`/`f_pH` is binding and which knob to move. Oxygen faults still alarm. |
| Root cause says **forming · 45 / 120 min** | Working as intended (§7b.4): the evidence window has not closed, so the name is withheld. It appears on its own when the window closes; *Show the provisional answer anyway* forces it, or switch to the probes tab (30 min instead of 120). |
| The named ethanol cause turns out to be wrong | Expected at a measured ~17 % rate: 120 min is a release point, not a certainty guarantee (§7b.4). **Read the tie set, not the top row** — the right cause is usually inside it. Wait to 150 min for ~97 %, or use the probes tab, which is 100 % on the *loop*. |
| Root cause names the wrong oxygen **rung** (`stop` vs `clog`) | A ceiling, not a latency. A DO probe reads the same floor under kLa = 0 % and 1 %, so the probes guess between the two rungs — while still naming the oxygen **loop** correctly every time. Ethanol separates them, slowly. |
| A `RuntimeWarning` on the server about a crippled reference | The setpoint-quality check (§4.7 #3): the culture is below 50 % of its peak rate at the setpoint the controller holds. The run proceeds; the results are just not worth much. |
| Blank page, no charts | Plotly CDN blocked — check connection or vendor Plotly locally (see §3). |
| "WebSocket error" toast | Server not running, or a proxy is blocking WS. Confirm `uvicorn` is up on the same host/port. |
| Curves don't move | Run not started, or `sec_per_hour` set very high. Check the clock in the toolbar. |
| Fault rejected | `t_fault` is in the past relative to the sim clock — pick a future time. |
| DO panel sits near 100 %, barely moves | Expected in a healthy run — DO holds near saturation. Inject `agitation_aer_stop` (or raise `qO2_max`) to drive it down; the dip is deepest during peak growth. |
| pH dips then climbs back to 5 after a fault | Expected: the PI controller dips (proportional lag) then the integral drives pH back to setpoint. Tune `Kc_pump` (dip depth) / `Ki_pump` (recovery speed). |
| Acid-pump-stuck doesn't crash pH | Expected when `F_A_stuck < F_B` — pH dips then the PI controller recovers it to setpoint. Set stuck flow > 210 mL/h to force a crash. |
| Slow / stuttering plots | Very small `dt_min` with EKF on is compute-heavy; increase `dt_min` or lower `t_final`. |
| Root cause says "tie · 4 causes" | Expected on the **ethanol only** tab: those faults all end in "growth stalled" and ethanol cannot separate them (§7b.3). Switch to **+ T/pH/DO probes**. |
| Root cause names the wrong mechanism on the ethanol tab | Read the tie list, not the top row — at margins below 2.3 nats the ordering is not evidence. |
| "The bank explains the alarm best with no fault at all" | No hypothesis, at the severity it assumes, matches this batch — e.g. a heater derated to 40 % against a dead-heater hypothesis. Adjust `HYPOTHESIS_SEVERITY` in `mmae.py`. |
| A cause is flagged *silent on ethanol* | It never separates from "no fault" on the measured channel, so it cannot be isolated from ethanol at any threshold. It may still be wide open to the probes. |
| Re-diagnose button greyed out | It needs a live run (the bank replays on the server) and an alarm. |
| Port already in use | `uvicorn ... --port <other>`. |

**Backend sanity check without the browser:**
```bash
python -c "from backend.simulator import Simulation; s=Simulation({}); \
[s.step() for _ in range(600)]; print(s.rows[-1])"
```

---

## 12. Provenance & references

Models were ported cell-by-cell from the project notebooks; see the module
docstrings for the exact source cells.

- **Yousefi-Darani et al. (2020)** — YD biokinetic model, initial/operating conditions.
- **Rosso (1993); Rosso et al. (1995)** — cardinal-temperature (and, by extension, cardinal-pH) modifiers.
- **Salvadó et al. (2011)** — cardinal-T values for *S. cerevisiae*.
- **Cooney & Wang (1969); Roels (1980)** — metabolic heat yield `Y_Q/X`.
- **Doran (2013), Ch. 9** — stirred-tank energy balance (thermal ODE form).
- **Bar-Shalom, Li & Kirubarajan (2001)** — EKF / filter-consistency background.
- Minifors 2 hardware values (`V`, `Q_max`, pump `F_max`) — Infors HT manual.

For the full fault-injection design rationale (two-tier detection, organism ×
severity coupling), see the companion `fault_injection.md` in the project root.

---

*Runs locally only. No telemetry, no network calls except the Plotly CDN asset.*
