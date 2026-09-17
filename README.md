# A Digital Shadow of a Baker's-Yeast Batch Bioreactor — Extended Kalman Filter for State Estimation and Fault Identification

Master's thesis, University of Rostock.
Author: Rowshanak Hosseinzadehattar. 
Supervisor: `Dr. rer. nat. Sebastian Bader`.

This repository packages the modelling and estimation work behind the thesis into
8 notebooks and one reference README. The scientific core is the
**ethanol-only Extended Kalman Filter (EKF)** of Yousefi-Darani et al. (2020)
`[CITATION: Yousefi-Darani, Paquet-Durand & Hitzmann 2020]`, which we (i) reproduce,
(ii) test against hypotheses, (iii) test against a panel of
independent literature yeast models, (iv) extend into a fault-diagnosis digital shadow,
(v) use to answer a concrete instrumentation question — *are the temperature and
pH sensors worth buying?* — and finally (vi) invert the filter's own blindness into a
**detector** that says *that* a fault happened, and (vii) a **hypothesis bank** that says
*which* one.


---

## 1. How to read this repository

| # | Notebook | Question it answers | Source research notebook |
|---|----------|---------------------|--------------------------|
| 1 | `1_ekf_baseline_yd.ipynb` | Can we reproduce the YD 2020 EKF and its paper figures with our own ground-truth simulator? 
| 2 | `2_hypotheses_observability.ipynb` | What did the EKF's own hypotheses (noise/drift robustness; observability of CO₂) reveal about its limits? 
| 3 | `3_generalization.ipynb` | How far does an ethanol-only EKF transfer across different yeast cultivations, and when does a second sensor become necessary? 
| 4 | `4_fault_injection.ipynb` | How does the shadow behaive in presence of faults (heater, leak, pH-pump, aeration)? 
| 5 | `5_sensor_value_TpH.ipynb` | Are the temperature and pH sensors worth buying for state estimation? 
| 6 | `6_fault_detection.ipynb` | Can the EKF detect an actuator fault at all? 
| 7 | `7_mmae_fault_isolation.ipynb` | Given an alarm, *which* fault caused it? 
| 8 | `8_probes.ipynb` | Which *single* probe pays for itself — T, pH or DO — for detection and for isolation? 

- `data/` — all ground-truth and result CSVs the notebooks read.
- `figures/` — final figures reproduced by the notebooks.


**Dependencies:** Python ≥ 3.10, `numpy`, `scipy`, `sympy`, `matplotlib`, `pandas`.

---

## 2. The scientific core — the YD ethanol-only EKF

### 2.1 State vector and measurement

The estimator is a **5-state EKF** on

```
x = [ X,  G,  E,  µmax_G,  µmax_E ]
      g/L g/L g/L  1/h      1/h
```

where `X` = biomass, `G` = glucose, `E` = ethanol, and `µmax_G`, `µmax_E` are the two
maximum specific growth rates carried as **augmented (estimated) states** so the filter
can adapt kinetics online. The **only measurement is ethanol** (an off-gas proxy),
`H = [0, 0, 1, 0, 0]`, sampled at a 5-minute cadence.

### 2.2 Process model (diauxic Monod)

Baker's yeast grows diauxically: it consumes glucose first, then re-consumes the
ethanol it produced once glucose is depleted. The kinetics carry this switch explicitly:

```
µ_G  = µmax_G · G/(K_G + G)                       # glucose-limited growth
supp = max(1 − µ_G/µmax_G, 0)                      # diauxic suppression (≈0 while glucose high)
µ_E  = µmax_E · E/(K_E + E) · supp                 # ethanol consumption, switches on as G→0

dX/dt = (µ_G + µ_E) · X
dG/dt = −(µ_G / Y_XG) · X
dE/dt =  (µ_G · Y_EG / Y_XG − µ_E / Y_XE) · X      # +production (glucose) −reconsumption
dµmax_G/dt = dµmax_E/dt = 0                        # constant in truth; estimated in the EKF
```

**Parameter values** (YD MATLAB-code set) `[CITATION: Yousefi-Darani et al. 2020]`:
`K_G = K_E = 0.1 g/L`, `Y_XG (Ygx) = 0.15`, `Y_EG (Yge) = 0.34`, `Y_XE (Yex) = 0.43`.


## 3. Notebook 1 — Baseline reproduction

**Goal:** reproduce the YD 2020 EKF and its published figures with our own deterministic
ground-truth simulator, establishing the non-circular validation baseline (the plant is
integrated with a high-accuracy solver; the EKF sees only noisy ethanol samples).

- **Ground-truth simulator:** the diauxic-Monod ODEs of §2.2 integrated with
  `solve_ivp` (LSODA, `rtol≈1e-7…1e-8`, `atol≈1e-9`) to produce a reference batch;
  synthetic ethanol measurements are drawn at 5-min cadence with `σ_E`.
- **Reproduced figures:** the YD paper's Figure 4 (state tracking) and Figure 6
  (parameter convergence). Reference digitized paper data is in `data/data_paper.csv`;
  our simulated reference is `data/data.csv`.
- **Result to expect:** ethanol tracks to the noise floor; `µmax_E` converges essentially
  exactly (dense ethanol data in the E-phase pins it), while `µmax_G` retains a small
  low-side bias (~17 %) — a genuine observability limit, not a bug, because only one
  glucose-phase window informs it.

`[CITATION: Yousefi-Darani, Paquet-Durand & Hitzmann 2020 — original EKF, parameters, figures]`

---

## 4. Notebook 3 — Generalization across literature models

**Strategy (non-circular external validation):** instead of testing the EKF on the plant
it was designed for, we test it against **nine independently-fitted literature yeast
models**, each implemented as a faithful `solve_ivp` simulator that reproduces the source
paper's own figures and writes a ground-truth CSV (`data/sim1..sim9_groundtruth.csv`).
The plants were *not* built to flatter the filter.

**"Only initials change" rule:** across plants we keep the chosen filter fixed and change
only `G0` (initial sugar), `X0` (inoculum guess), run length, and noise scale. We
deliberately **do not retune yields per plant** — a yield mismatch that breaks recovery is
itself a result.

### 4.1 The nine models + control

| # | Tag | Paper | Organism / substrate | Stress on the EKF |
|---|-----|-------|----------------------|-------------------|
| 1 | sim1 | Damayanti (2022) | immobilized alginate beads | two glucose pools → carbon balance does not close |
| 2 | sim2 | Salakkam (2023) | sweet sorghum, high-gravity | substrate + product inhibition |
| 3 | sim3 | Pinheiro (2017) | cashew apple juice | Ghose–Tyagi + cell death (Xv/Xd) |
| 4 | sim4 | Wang (2004) | apple wine, 15 °C | logistic growth + 15 h ethanol lag |
| 5 | sim5 | Konopacka (2019) | low-gravity free-cell | Gompertz ethanol decoupled from sugar |
| 6 | sim6 | Ariyajaroenwong (2016) | immobilized on sorghum stalks | double Andrew–Levenspiel inhibition |
| 7 | sim7 | Zinnai (2014) | *S. bayanus*, glucose+fructose | two hexoses; **non-growing biocatalyst** (biomass fixed) |
| 8 | sim8 | Salazar (2023) | *Kluyveromyces marxianus* | **different genus**; ethanol rises then degrades |
| 9 | sim9 | González-Hernández (2022) | metabolic switching, aerobic | respiration: ethanol rises then is respired |
| 10 | — | extended plant | *S. cerevisiae*, near-nominal | **control**, biomass tracks because parameters near nominal |

`[CITATION: Damayanti 2022]` `[CITATION: Salakkam 2023]` `[CITATION: Pinheiro 2017]`
`[CITATION: Wang 2004]` `[CITATION: Konopacka 2019]` `[CITATION: Ariyajaroenwong 2016]`
`[CITATION: Zinnai 2014]` `[CITATION: Salazar 2023]` `[CITATION: González-Hernández 2022]`

Two EKF structures are used: a **generic simple EKF** `[X, S, E, µmax]` (single substrate,
monotonic ethanol) for the fermentative six; and the **YD-exact ethanol-switch EKF** of
§2 for sim7/sim8/sim9 (ethanol non-monotonic or respiration present). Simple-EKF yields
use the universal Gay-Lussac value `Yge_g = 0.46`, `Ygx_g = 0.10`.

---

## 5. Notebook 4 — Fault-injection digital twin

This notebook extends the YD plant into a fault-diagnosis twin by adding three
environmental dynamics — **temperature, pH, and dissolved oxygen** — each with a
first-principles ODE, an actuator/controller, and a fault. 


### 5.1 The environmental ODEs 

**Temperature — a stirred-tank energy balance** (textbook form, Doran 2013 Ch. 9)
`[CITATION: Doran 2013, Bioprocess Engineering Principles, Ch. 9]`:

```
Q_heat  = clip(K_p·(T_set − T), 0, Q_max)         # one-sided P-control (heater only, no active cooling)
Q_loss  = UA·(T − T_amb)                           # Newton's law of cooling
Q_metab = Y_QX·(µ_G + µ_E)·X·V                     # metabolic heat release
dT/dt   = (3600·(Q_heat − Q_loss) + Q_metab) / (V·ρCp)
```

*Why this form:* it is the standard first-principles energy balance; the heater is modelled
one-sided because the Minifors has a heated jacket. The metabolic
heat yield `Y_QX` links growth to self-heating.
Sources: `Y_QX` from Cooney & Wang (1969) / Roels (1980); `U` from Müller et al. (2013).
`[CITATION: Cooney & Wang 1969]` `[CITATION: Roels 1980]` `[CITATION: Müller et al. 2013]`

**pH — acidogenesis opposed by throttled PI acid/base pumps:**

```
err   = pH_set − pH                                # >0 → too acidic → add base
I    += err · dt                                    # error integral (with anti-windup)
u     = K_c·err + K_i·I                              # signed pump command [mL/h]
F_B   = clip( u, 0, F_B_max)                         # base (raises pH)
F_A   = clip(−u, 0, F_A_max)                         # acid (lowers pH)
yeast_acid = α_meta·(µ_G + µ_E)·X                    # cells acidify the broth
dpH/dt = −yeast_acid + α_pump·F_B − α_pump·F_A
```

*Why PI (decision and reason):* pure-P control leaves a steady-state offset and a bang-bang
relay limit-cycles under a constant disturbance. PI gives a smooth **dip-and-recover**
response — pH dips when a fault hits (proportional term still small), then the integral
restores the setpoint — which is exactly the readable fault signature we want. `K_c` sets
dip depth, `K_i` sets recovery speed, anti-windup freezes the integral while a pump
saturates. (The real Minifors pumps are digital ON/OFF; PI is a deliberate simplification
for a clean trace.) The acidogenesis coefficient `α_meta` was calibrated up to 0.4 so the
uncontrolled drift `Δph ≈ α_meta·ΔX ≈ 0.5 pH` is realistic and gives the base pump real work.
Sources: cardinal-pH from Furukawa et al. (1983); PI gains are controller tuning, not
measured constants .
`[CITATION: Furukawa et al. 1983]`

**Dissolved oxygen — a free-floating supply/demand balance (no DO controller):**

```
supply  = kLa_eff·(DO* − DO)                              # OTR  [%-sat/h]
demand  = qO2_max·(µ_G + µ_E)/(µmax_G + µmax_E)·X          # OUR  [mg/(L·h)]
dDO/dt  = supply − demand / DO_sat_mgL · 100              # %-sat/h (clamped ≥ 0)
```

*Why this form:* the YD setup holds agitation/aeration constant with no DO control, so DO
is a free balance of oxygen transfer against respiratory uptake. Crucially, **`kLa` is not a
constant**: `compute_kLa(N_rpm, Q_air, V)` derives it from agitation and aeration via the
**Van't Riet correlation** `kLa = C·(Pg/V)^a · v_s^b`, so an aeration/agitation fault
physically starves the supply term. The shaft carries **two 6-blade Rushton turbines**
(Minifors manual), modelled with **additive power** (`P0 = 2·Np·ρ·N³·D⁵`, `Np ≈ 5` each),
giving healthy kLa(500 rpm, 3.5 L/min) ≈ **44.8 1/h** and ≈ **10.5 1/h** at the 250-rpm
agitation fault. DO is expressed in **% saturation** to match a Minifors pO₂ probe.
`DO*` (saturation reference) = 100 %. OUR is a specific-uptake heuristic scaled
by growth activity.
Sources: Van't Riet / Hrnčiřík & Kohout (2024); qO2_max from Sonnleitner & Käppeli (1986).
`[CITATION: Van't Riet 1979]` `[CITATION: Hrnčiřík & Kohout 2024]` `[CITATION: Sonnleitner & Käppeli 1986]`


---


## 7. Notebook 5 — Fault detection with a healthy digital shadow




---

## 8. Notebook 6 — MMAE: which fault?

Notebook 5 says *that* something broke and roughly *when*. Notebook 6 says *what*, with a
**multiple-model adaptive estimation** bank: one EKF per candidate cause, each *believing* its own
fault started at the estimated onset `t̂₀`, all fed the **same real measurements**, scored by Bayes
over their innovations.

`p_j(k) ∝ p_j(k−1) · N(ν_j(k); 0, S_j(k))` — in logs, one `cumsum`.

The candidate faults are injected into the **filter's model**, not into a second simulator. There
is only ever one batch of data, and the bank competes to explain it.


---
## 9. Notebook 7 — Evaluation: how fast is the alarm, and how fast is the diagnosis?

The last notebook puts a number on the two questions the previous ones answered only
qualitatively: **how long after a fault does the CUSUM alarm**, and **how long after the
alarm does the MMAE verdict stop changing.

**Setup:** 8 scenarios = 4 actuator loops × {dead, running-but-wrong}, every fault injected
at `t = 2 h` into a 10 h batch, 5 noise seeds each. Two measurement sets are compared:
**ethanol only** (the soft sensor) and **ethanol + T/pH/DO probes at the same 5-min
cadence**, the probes used as monitors only — they are never fed back into the shadow, so
each residual stays a clean fault signal. Settling time is the smallest offset after the
alarm at which the verdict is correct at *every* later window.

| minutes after the fault | ethanol only | + T/pH/DO @ 5 min |
|---|---|---|
| alarm raised | 6/8 detected, mean 86, max 116 (heater 6 %) | 7/8 detected, mean 12, max 51 (base pump dead) |
| diagnosis settles (after alarm) | 6/6 settle, mean 35, max 90 | 5/6 settle, mean 0 |
| fault → settled diagnosis | mean 121, max 184 | mean 5 |

**What it shows.** On the six scenarios both configurations see, the probes cut the alarm
delay from 86 min to 5 min (~17×, i.e. the very first sample after the fault), and the
diagnosis is then already correct at the alarm. Two scenarios never alarm on ethanol, for
different reasons: `acid_pump_off` and `base_pump_off`. The pH probe recovers it regardless
(51 min). The one regression: `agitation_aer_stop` settles on the wrong cause with probes,
where ethanol alone gets it right at 60 min — the DO channel saturates and the two
aeration hypotheses become hard to separate.

**Outputs:** `figures/8_evaluation_latency.png` and
`evaluation_{latency,cusum_runs,mmae_windows}.csv`.

---
## 9b. `probes.ipynb` — one probe at a time

Notebook 7 compares ethanol alone against the whole T/pH/DO rack. `probes.ipynb` splits
that rack up, because the instrumentation budget is spent per probe: it runs the same rig
under **five measurement systems** — `E`, `E+T`, `E+pH`, `E+DO`, `E+all` — every channel at
the same 5-min cadence, with each system carried through **both** stages (its own CUSUM
bank for detection, the MMAE bank scored on its own channels for isolation). Part A is
copied verbatim from notebook 7, so the only thing that varies is the measurement set. The
one structural change is in the scorer: the probe likelihood is accumulated **per channel**
instead of lumped, so one set of replays scores all five systems.

**Result: each probe covers its own loop, and almost nothing else** (fault at t = 0.5 h,
one batch and one seed per scenario).

| fault | alarm E | alarm with the probe that helps | isolation |
|---|---|---|---|
| heater dead / underpowered | 40 / 45 min | **T**: 5 min | T carries ~7 000 nats; E's own margin is 0.2 |
| acid pump stuck ON | 40 min | **pH**: 5 min | pH carries ~2.9e5 nats |
| base pump dead | *never alarms* | **pH**: 40 min | invisible without pH |
| base pump stuck ON | 35 min | **pH**: 5 min | E settles at 335 min, E+pH at 5 |
| agitation/aeration stopped / restricted | 85 / 185 min | **DO**: 10 min | DO carries 2.7e4 / 2.0e5 nats |
| acid pump dead | never | never | state-identical to healthy — no instrument helps |

Isolation summary over the 8 scenarios: verdict correct 6 (E), 6 (E+T), 7 (E+pH), 6 (E+DO),
7 (E+all); *identified* at the 2.3-nat bar 2 / 5 / 4 / 3 / 7; mean settle time from the
fault 172 / 124 / 75 / 100 / **20** min. The single probes each fix their own loop and add
little elsewhere, and only the full rack is above the identification bar everywhere it
alarms.

`E` and `E+all` here **are** notebook 7's `E` and `E+P` — same channels, detector, bank,
batch and seeds — so §C.4 re-derives them and diffs them cell by cell against notebook 7's
published Part C (alarm and settle to the minute, margins to 0.5 %); it prints `ALL AGREE`.
That check is why the window ladder `OFFSETS_MIN` is carried over verbatim: settle time is
quantised to those rungs, so a trimmed ladder makes a system fall through to the next rung
it owns rather than measuring anything about the system.

**Outputs:** `figures/probes_{alarm_and_verdict,evidence_by_channel}.png` and
`data/probes_{alarm_delay,verdict_delay,margin_nats,evidence_by_channel}.csv`.

---
## 10. References (to be completed)

Fill in full bibliographic entries; the `[CITATION: …]` markers above indicate where each is used.

- `[CITATION: Yousefi-Darani, Paquet-Durand & Hitzmann (2020)]` — original ethanol-only EKF, YD biokinetic parameters, initial/operating conditions, reproduced figures.
- `[CITATION: Rosso (1993); Rosso et al. (1995)]` — cardinal-temperature (and cardinal-pH) growth modifier.
- `[CITATION: Salvadó et al. (2011)]` — cardinal-T values for *S. cerevisiae*.
- `[CITATION: Furukawa et al. (1983)]` — cardinal-pH values; DO Monod K_O.
- `[CITATION: Cooney & Wang (1969)]`, `[CITATION: Roels (1980)]` — metabolic heat yield Y_QX.
- `[CITATION: Doran (2013)]` — *Bioprocess Engineering Principles*, Ch. 9, stirred-tank energy balance.
- `[CITATION: Müller et al. (2013)]` — heat-transfer coefficient U in single-use stirred bioreactors.
- `[CITATION: Van't Riet (1979)]`, `[CITATION: Hrnčiřík & Kohout (2024)]` — kLa correlation for the DO ODE.
- `[CITATION: Sonnleitner & Käppeli (1986)]` — specific O₂ uptake rate qO2_max.
- `[CITATION: Isermann (2006)]` — *Fault-Diagnosis Systems*, fault taxonomy.
- `[CITATION: Page (1954)]` — CUSUM foundation for L-a detection.
- `[CITATION: Bar-Shalom, Li & Kirubarajan (2001)]` — EKF / filter-consistency (NIS, innovation whiteness).
- `[CITATION: Morris (1991)]` — elementary-effects sensitivity screening.
- Generalization panel: `[CITATION: Damayanti (2022)]`, `[CITATION: Salakkam (2023)]`, `[CITATION: Pinheiro (2017)]`, `[CITATION: Wang (2004)]`, `[CITATION: Konopacka (2019)]`, `[CITATION: Ariyajaroenwong (2016)]`, `[CITATION: Zinnai (2014)]`, `[CITATION: Salazar (2023)]`, `[CITATION: González-Hernández (2022)]`.

---

