"""
core.py — biokinetic model primitives and default parameters.

Ported from the notebooks:
  * carbon_fluxes (Crabtree)      — 7_mmae_fault_isolation.ipynb, A.3
    (replaced monod_rates, which folded f_DO into phi — wrong sign)
  * re-anchored Rosso f_T         — fault-injection.ipynb, Cell 1
  * re-anchored cardinal f_pH     — fault-injection.ipynb, Cell 9

Everything here is pure / parameter-driven so the same functions serve the
plant (ground truth) and the EKF (estimator).
"""
from __future__ import annotations

# ─────────────────────────────────────────────────────────────────────────
# Default parameters — grouped exactly the way the UI panels are grouped.
# Kinetic defaults use the "code" yield set (0.15 / 0.34 / 0.43).
# Thermal / pH defaults track the current fault-injection.ipynb values.
# ─────────────────────────────────────────────────────────────────────────
DEFAULTS: dict = {
    "kinetics": {
        "KG":  0.1,    # glucose half-saturation [g/L]
        "KE":  0.1,    # ethanol half-saturation [g/L]
        "Ygx": 0.15,   # biomass / glucose yield  (OVERALL; see the split below)
        "Yge": 0.34,   # ethanol / glucose yield  (OVERALL)
        "Yex": 0.43,   # biomass / ethanol yield
        # ── Respiro-fermentative split (Crabtree) ────────────────────────
        # Oxygen does NOT scale the growth rate. It sets the FATE of glucose
        # carbon: respiration has a hard capacity, and glucose arriving faster
        # than that capacity overflows into fermentation. Losing oxygen
        # therefore makes MORE ethanol, not less — the opposite sign to a T or
        # pH fault, and that sign is what tells them apart on ethanol alone.
        #
        # These are a DECOMPOSITION of the inherited YD yields, not a
        # replacement. With the fermentative branch at the theoretical
        # 0.51 g E/g G:
        #     phi_ferm = Yge / YGE_FERM = 0.34 / 0.51 = 2/3
        #     YGX_OX*(1 - phi_ferm) + YGX_RED*phi_ferm = Ygx = 0.15 -> YGX_OX = 0.25
        # i.e. the YD numbers already ARE a 1/3-respiratory, 2/3-fermentative
        # blend; this writes the blend down so the halves can move apart when
        # oxygen goes away.  (7_mmae_fault_isolation.ipynb, A.2)
        "YGE_FERM": 0.51,     # g ethanol / g glucose fermented (theoretical max)
        "YGX_OX":   0.25,     # g X / g glucose respired
        "YGX_RED":  0.10,     # g X / g glucose fermented
        # *** STALE CALIBRATION, flagged: fitted so the healthy batch peaked at
        # E = 1.968 g/L; it now peaks near 2.29 after kLa became geometry-derived
        # and the respiratory cap became phi_g-scaled. Left alone deliberately —
        # it should be anchored on a measured respiratory capacity, not on a
        # target ethanol peak. No result depends on the absolute ethanol level.
        "C_RESP":   0.2573,   # respiratory capacity [g glucose / g X / h]
        "YO_E_REL": 2.0,      # O2 per g ethanol respired, rel. to per g glucose
        # CAVEAT, inherited: kLa, qO2_max and DO_sat_mgL are not stoichio-
        # metrically consistent. Respiring C_RESP g glucose/gX/h needs ~274 mg
        # O2/gX/h, ~53x this model's peak OUR. K_O2 keeps the INHERITED demand
        # scale (which is what makes DO drain under a kLa fault) and only
        # redirects it onto the oxidative flux. Re-anchoring the oxygen balance
        # against a measured OUR is a separate job.
        "K_O2":     9.050,    # mg O2 per g oxidised substrate
    },
    "initial": {
        "X0":   2.5,    # biomass [g/L]
        "G0":   5.0,    # glucose [g/L]
        "E0":   0.0,    # ethanol [g/L]
        "muG0": 0.15,   # µmax_G [1/h]
        "muE0": 0.08,   # µmax_E [1/h]
    },
    "thermal": {
        "V":      1.35,     # working volume [L]
        "rhoCp":  4180.0,   # heat capacity per volume [J/L/K]
        "Y_QX":   12000.0,  # metabolic heat yield [J/g]
        "T_set":  30.0,     # setpoint [°C]
        "T_amb":  22.0,     # ambient [°C]
        "UA":     18.0,     # heat-loss coefficient [W/K]  (U·A = 242·0.07)
        "K_p":    50.0,     # heater P-gain [W/K]
        "Q_max":  630.0,    # heater saturation [W]  (Minifors manual)
        "T_min":  5.0,      # Rosso cardinal min [°C]
        "T_opt":  30.0,     # Rosso cardinal opt [°C]
        "T_max":  40.0,     # Rosso cardinal max [°C]
    },
    "ph": {
        "pH_set":     5.0,
        # Cardinal pH = the 'robust' organism profile of
        # 7_mmae_fault_isolation.ipynb (A.2). The previous triple here was
        # (3.0, 5.0, 12.0), which is not wrong under n = 1 but does not match
        # the notebook, so the two artefacts scored different biology.
        # The notebook's 'sensitive' profile is (4.6, 5.0, 6.0) with cardinal T
        # (20, 30, 40) — override both groups together to run that organism.
        # *** CONFIRM *** these cardinals: ph_modifier_handoff.md traces their
        # magnitude to Hinga (2002), a review of MARINE PHYTOPLANKTON. It
        # supports the shape and order of magnitude, not yeast values.
        "pH_min":     2.5,     # cardinal min
        "pH_opt":     4.7,     # cardinal opt (deliberately != pH_set)
        "pH_max":     8.0,     # cardinal max
        "alpha_meta": 0.4,     # acidogenesis [pH/h per g/L·h⁻¹] — uncontrolled
                               # drop ≈ alpha_meta·ΔX ≈ 0.5 pH over a batch
        "alpha_pump": 0.2563,  # buffer coefficient [pH per mL 0.1 M titrant]
        "F_A":        210.0,   # acid-pump max (saturation) flow [mL/h]
        "F_B":        210.0,   # base-pump max (saturation) flow [mL/h]
        "Kc_pump":    120.0,   # PI proportional gain [mL/h per pH] — sets the
                               # depth of the initial pH dip when a fault hits
        "Ki_pump":    250.0,   # PI integral gain [mL/h per pH·h] — drives pH
                               # slowly back to setpoint after the dip (no offset)
        "DT_CTRL_s":  2.0,     # controller/integration tick [s]
    },
    "oxygen": {
        "N_rpm":      500.0,   # agitation [rpm] — Van't Riet input & agitation-fault target
        "Q_air_Lmin": 3.5,     # aeration [L/min] — Van't Riet input & aeration-fault target
        "DO_star":    100.0,   # O₂ saturation reference [% sat] (probe = 100%)
        "DO0":        100.0,   # initial DO [% sat] — saturated at inoculation
        "qO2_max":    8.0,     # specific O₂ uptake [mg O₂/g/h] — *confirm* (Sonnleitner & Käppeli 1986)
        "K_O":        3.0,     # O₂ Monod half-saturation [% sat] (≈0.2 mg/L / 6.7)
        "DO_sat_mgL": 6.7,     # mg/L ↔ %-sat conversion constant [mg/L]
    },
    "ekf": {
        "on":          True,
        # Design A (default) is the 8-state HEALTHY SHADOW: T/pH/DO are
        # estimated from the healthy ODEs and the innovation is the fault
        # signal. Design B is the legacy 5-state filter that takes measured
        # T/pH/DO as inputs and therefore follows the fault — kept so the two
        # architectures can be compared side by side.
        "old_design":  False,
        "cadence_min": 5.0,    # ethanol measurement cadence [min]
        # Two DIFFERENT numbers on purpose. sigma_E is what the filter ASSUMES
        # (it sets R); sigma_E_sensor is what the sensor ACTUALLY does. Tying
        # them together would make a degrading sensor invisible — the filter
        # would simply trust it less. Raising sigma_E_sensor alone is how you
        # sweep sensor quality and watch the estimate degrade against truth.
        
        "sigma_E":        0.3162,  # ethanol σ the filter assumes → R [g/L]
        "sigma_E_sensor": 0.05,    # real ethanol sensor noise σ [g/L]
        "cusum_k":     0.15,   # CUSUM slack     [units of sqrt(S)]
        "cusum_h":     1.5,    # CUSUM threshold [units of sqrt(S)]
        "band":        False,  # show ±σ uncertainty band
        # ── MMAE root-cause bank (mmae.py) ──────────────────────────────
        # OFF is the notebook setting: the bank scores on ETHANOL ALONE, which
        # is the honest question "what can a soft sensor tell you". At the
        # severities where these faults are detectable at all, that answer is
        # usually a tie — every one of them ends in "growth stalled", and
        # ethanol cannot say which loop stalled it. Turning the probes ON adds
        # the T/pH/DO instruments the rig actually has to the likelihood, which
        # breaks the tie because that is exactly where the causes differ. It is
        # a different question, not a better tuning, so it is a switch.
        "mmae_probes":   False,
        "sigma_T_probe":  0.1,   # Pt100 [°C]
        "sigma_pH_probe": 0.02,  # glass electrode [-]
        "sigma_DO_probe": 1.0,   # pO₂ probe [% sat]
    },
    "sim": {
        "t_final":       10.0,  # batch length [h]
        "sec_per_hour":  5.0,   # wall-clock seconds per simulated hour
        "dt_min":        1.0,   # integration / frame step [min]
    },
}


# ─────────────────────────────────────────────────────────────────────────
# Respiro-fermentative carbon split (Crabtree).
#
# `phi_g = f_T * f_pH` scales the glucose uptake RATE; oxygen does not appear
# there. Oxygen caps the RESPIRATORY branch, so glucose arriving faster than
# `C_RESP * f_DO` overflows into fermentation. Ethanol is a respiratory
# substrate, so without oxygen it cannot be consumed back either.
#
# Ported from 7_mmae_fault_isolation.ipynb (A.3, `carbon_fluxes`), which is the
# canonical plant. The previous form here multiplied f_DO into phi alongside
# f_T and f_pH; that made a DO fault reduce ethanol, the opposite of the
# Crabtree sign, and it is the reason the UI and the notebook disagreed.
# ─────────────────────────────────────────────────────────────────────────
def carbon_fluxes(G, E, mu_max_G, mu_max_E, kin, phi_g=1.0, f_do=1.0):
    """Split glucose uptake between respiration and fermentation.

    Returns (q_ox, q_fm, mu_E, mu) — glucose respired, glucose fermented,
    specific ethanol uptake, and total specific growth rate, per g X per h.
    """
    KG, KE, Ygx = kin["KG"], kin["KE"], kin["Ygx"]
    if mu_max_G > 1e-9:
        q_G = (mu_max_G * G / (KG + G) * phi_g) / Ygx      # total glucose uptake
        suppression = 1.0 - (q_G * Ygx) / mu_max_G          # diauxic repression
    else:
        q_G, suppression = 0.0, 1.0
    # Respiratory capacity is scaled by phi_g TOO. Respiration is enzymatic, so
    # cold or acid slows it exactly as it slows uptake; leaving the cap at its
    # 30 degC value made a chilled culture MORE respiratory, which is backwards
    # and had a nasty consequence: slow the culture ~3.5x and q_G fell below a
    # fixed cap, fermentation stopped entirely, and ethanol - the only
    # measurement - went silent. A heater fault that costs 33 % of the biomass
    # then produced 0.6 sigma of ethanol signal and could not be detected at any
    # threshold, where the pre-Crabtree model saw it at 19 sigma.
    # Scaling the cap keeps the fermentative FRACTION invariant to T and pH,
    # while OXYGEN still sets carbon fate on its own through f_do - which is the
    # whole point of the Crabtree rework and is untouched by this.
    q_ox = min(q_G, kin["C_RESP"] * phi_g * f_do)
    q_fm = q_G - q_ox                        # overflow -> fermentation

    # ethanol is a RESPIRATORY substrate: no oxygen, no consumption
    mu_E = (mu_max_E * E / (KE + E) * suppression * phi_g * f_do
            ) if mu_max_E > 1e-9 else 0.0

    mu = q_ox * kin["YGX_OX"] + q_fm * kin["YGX_RED"] + mu_E
    return q_ox, q_fm, mu_E, mu


def growth_rates(X, G, E, mu_max_G, mu_max_E, kin, phi_g=1.0, f_do=1.0):
    """(mu_G, mu_E) for display/CSV: growth carried by glucose vs by ethanol."""
    q_ox, q_fm, mu_E, mu = carbon_fluxes(G, E, mu_max_G, mu_max_E, kin, phi_g, f_do)
    return mu - mu_E, mu_E


# ─────────────────────────────────────────────────────────────────────────
# Re-anchored Rosso cardinal modifier, reused for both T and pH.
#
# THE TWO FORMS ARE NOT INTERCHANGEABLE (fault_detection_CU.ipynb Cell 6).
# n = 2 (CTMI) is the cardinal *temperature* model. Its denominator has a root
# inside (v_min, v_max) whenever the optimum sits below the midpoint of the
# window, and the clip to [0, 1] hides that pole as a step function. That is
# what happened to f_pH on the triple (3.0, 5.0, 12.0): arms of 2 and 7, pole at
# ~3.84, so f_pH returned exactly 0 below it and exactly 1 above. It carried no
# gradient at all — it could say "dead" or "fine" and nothing in between, which
# made a stuck acid pump indistinguishable from anything else that stops growth
# (the MMAE bank read acid_pump_stuck as the oxygen fault). pH uses n = 1,
# which has no interior pole and handles asymmetric arms. `_assert_monotone` in
# plant.build_modifiers is the regression test; it fails on the n = 2 triple.
# ─────────────────────────────────────────────────────────────────────────
PHI_FLOOR = 1e-6   # never return a hard zero: that hands the EKF a zero
                   # Jacobian and the filter can never correct its way back.


def _rosso_n2(v, v_min, v_opt, v_max):
    """Rosso (1993) CTMI, n = 2 — the cardinal TEMPERATURE model."""
    if v <= v_min or v >= v_max:
        return 0.0
    num = (v - v_max) * (v - v_min) ** 2
    den = (v_opt - v_min) * (
        (v_opt - v_min) * (v - v_opt)
        - (v_opt - v_max) * (v_opt + v_min - 2 * v)
    )
    return 0.0 if abs(den) < 1e-12 else num / den


def _rosso_n1(v, v_min, v_opt, v_max):
    """Rosso n = 1 — the pH form. No interior pole, asymmetric arms handled."""
    if v <= v_min or v >= v_max:
        return 0.0
    num = (v - v_min) * (v - v_max)
    den = num - (v - v_opt) ** 2
    return 0.0 if abs(den) < 1e-12 else num / den


def absolute_modifier(v, v_min, v_opt, v_max, n=2):
    """The modifier at `v` relative to the organism's TRUE optimum, un-anchored.

    `make_cardinal_modifier` normalises at the SETPOINT, which makes f(setpoint)
    identically 1 — useful for reading "how far from where we hold it", useless
    for asking "is where we hold it any good". That second question needs the
    raw value scaled by the peak, and it is the one that catches a setpoint
    parked in a dead corner of the organism's window.
    """
    raw = _rosso_n1 if n == 1 else _rosso_n2
    peak = raw(v_opt, v_min, v_opt, v_max)
    return 0.0 if peak < 1e-12 else raw(v, v_min, v_opt, v_max) / peak


def make_cardinal_modifier(v_min, v_opt, v_max, anchor=None, n=2):
    """Build a re-anchored Rosso (1993) cardinal modifier on [v_min, v_max].

    `n` selects the form, and the choice is about the SHAPE OF THE WINDOW, not
    about which variable it models: 2 is valid only while the optimum sits at or
    above the window midpoint, 1 is pole-free for any arms. pH always takes 1;
    temperature takes 2 for the usual right-skewed window and falls back to 1 on
    a left-skewed one such as (27, 30, 40) — see `plant.build_modifiers`.

    Clipping happens AFTER re-anchoring, and only from below. Clipping to [0, 1]
    first — the previous behaviour — truncated the physically correct f > 1 that
    arises when the reactor sits nearer the optimum than the anchor does.
    """
    raw = _rosso_n1 if n == 1 else _rosso_n2

    # THE ANCHOR MUST LIE INSIDE THE VIABLE WINDOW, AND THIS USED TO BE A SILENT
    # FALLBACK. `a` is the raw modifier at the value the controller actually holds,
    # and the whole modifier is normalised by it. If it is ~0 the setpoint sits at
    # or outside a cardinal bound — the HEALTHY culture is dead at its own
    # setpoint — and everything downstream becomes nonsense in a way nothing
    # reports: the healthy reference never grows, so a fault that moves the
    # variable back INTO the window makes the culture grow BETTER than healthy,
    # and the detector, which compares against that dead reference, sees nothing
    # to alarm about. Observed with pH_set = pH_max = 5.0: healthy X_end 2.50
    # (i.e. no growth at all), base-pump-off X_end 4.00, no alarm.
  
    a = raw(v_opt if anchor is None else anchor, v_min, v_opt, v_max)
    if a < 1e-9:
        raise ValueError(
            f"anchor {anchor if anchor is not None else v_opt} lies outside the "
            f"viable window ({v_min}, {v_max}): the modifier is normalised at the "
            f"setpoint, so a setpoint the organism cannot grow at leaves the "
            f"HEALTHY reference dead. Move the setpoint inside the window, or "
            f"widen the window past it.")

    def f(v):
        return max(raw(v, v_min, v_opt, v_max) / a, PHI_FLOOR)

    return f


# ─────────────────────────────────────────────────────────────────────────
# Re-anchored oxygen limitation modifier f_DO(DO).
# Unlike f_T / f_pH (bell-shaped Rosso cardinals), oxygen limitation is a
# Monod saturation term DO/(K_O+DO). It is re-anchored so f_DO(DO_star) == 1,
# i.e. at saturation there is NO growth penalty; only the deficit below
# saturation damps growth. Returns a callable clipped to [0, 1].
# ─────────────────────────────────────────────────────────────────────────
def make_do_modifier(K_O, DO_star):
    """Build a re-anchored O₂ Monod modifier f_DO with f_DO(DO_star) == 1."""
    anchor = DO_star / (K_O + DO_star)
    if anchor < 1e-9:
        anchor = 1.0

    def f_DO(DO):
        DO = max(DO, 0.0)
        val = (DO / (K_O + DO)) / anchor
        return max(0.0, min(1.0, val))

    return f_DO
