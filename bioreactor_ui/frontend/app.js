"use strict";

// ── pretty labels / units for known parameters ──────────────────────────
const LABELS = {
  "kinetics.KG": ["K_G", "g/L"], "kinetics.KE": ["K_E", "g/L"],
  "kinetics.Ygx": ["Y_X/G", ""], "kinetics.Yge": ["Y_E/G", ""], "kinetics.Yex": ["Y_X/E", ""],
  "initial.X0": ["X₀", "g/L"], "initial.G0": ["G₀", "g/L"], "initial.E0": ["E₀", "g/L"],
  "initial.muG0": ["µmax_G₀", "1/h"], "initial.muE0": ["µmax_E₀", "1/h"],
  "thermal.V": ["V", "L"], "thermal.rhoCp": ["ρ·Cp", "J/L/K"], "thermal.Y_QX": ["Y_Q/X", "J/g"],
  "thermal.T_set": ["T_set", "°C"], "thermal.T_amb": ["T_amb", "°C"],
  "thermal.UA": ["U·A", "W/K"], "thermal.K_p": ["K_p", "W/K"], "thermal.Q_max": ["Q_max", "W"],
  "thermal.T_min": ["T_min", "°C"], "thermal.T_opt": ["T_opt", "°C"], "thermal.T_max": ["T_max", "°C"],
  "ph.pH_set": ["pH_set", ""], "ph.pH_min": ["pH_min", ""], "ph.pH_opt": ["pH_opt", ""],
  "ph.pH_max": ["pH_max", ""], "ph.alpha_meta": ["α_meta", ""], "ph.alpha_pump": ["α_pump", ""],
  "ph.F_A": ["acid flow (max)", "mL/h"], "ph.F_B": ["base flow (max)", "mL/h"],
  "ph.Kc_pump": ["gain K_c (P)", "mL/h/pH"], "ph.Ki_pump": ["gain K_i (I)", "mL/h/pH·h"],
  "ph.DT_CTRL_s": ["ctrl tick", "s"],
  "oxygen.N_rpm": ["agitation", "rpm"], "oxygen.Q_air_Lmin": ["aeration", "L/min"],
  "oxygen.DO_star": ["DO sat ref", "% sat"], "oxygen.DO0": ["DO₀", "% sat"],
  "oxygen.qO2_max": ["q_O₂ max", "mg/g/h"], "oxygen.K_O": ["K_O", "% sat"],
  "oxygen.DO_sat_mgL": ["DO sat", "mg/L"],
  "ekf.cadence_min": ["cadence", "min"],
  "ekf.sigma_E": ["σ assumed (R)", "g/L"], "ekf.sigma_E_sensor": ["σ sensor (true)", "g/L"],
  "ekf.cusum_k": ["CUSUM k (slack)", "√S"], "ekf.cusum_h": ["CUSUM h (threshold)", "√S"],
  "ekf.sigma_T_probe": ["σ T probe", "°C"], "ekf.sigma_pH_probe": ["σ pH probe", ""],
  "ekf.sigma_DO_probe": ["σ DO probe", "% sat"],
  "sim.t_final": ["batch length", "h"], "sim.sec_per_hour": ["wall s / sim-h", "s"],
  "sim.dt_min": ["step", "min"],
};
const BOOL_LABELS = {
  "ekf.on": "EKF on", "ekf.band": "±σ band",
  "ekf.old_design": "old design (B)",
  "ekf.mmae_probes": "RCA reads probes",
};

let DEFAULTS = null, CATALOGUE = null;
let ws = null, running = false;
let pendingFaults = [], faultLines = [];

// ── build config inputs ─────────────────────────────────────────────────
function field(group, key, val) {
  const id = `f-${group}-${key}`;
  if (typeof val === "boolean") {
    const lbl = BOOL_LABELS[`${group}.${key}`] || key;
    return `<div class="field"><label>${lbl}</label>
      <label class="switch"><input type="checkbox" id="${id}" ${val ? "checked" : ""}>
      <span class="slider"></span></label></div>`;
  }
  const [lbl, unit] = LABELS[`${group}.${key}`] || [key, ""];
  return `<div class="field"><label>${lbl} ${unit ? `<span class="unit">${unit}</span>` : ""}</label>
    <input type="number" step="any" id="${id}" value="${val}"></div>`;
}

function buildGroup(group) {
  const host = document.getElementById(`grp-${group}`);
  if (!host) return;
  let entries = Object.entries(DEFAULTS[group]);
  // EKF panel: the toggles, both σ's, and the detector tuning. The σ's are
  // separate on purpose — σ sensor is the real noise, σ assumed is what the
  // filter puts in R. Raise σ sensor alone to sweep sensor quality and watch
  // the NRMSE badge climb. Cadence stays at its default.
  if (group === "ekf")
    entries = entries.filter(([k]) =>
      ["on", "old_design", "band", "sigma_E", "sigma_E_sensor", "cusum_k", "cusum_h",
       "mmae_probes", "sigma_T_probe", "sigma_pH_probe", "sigma_DO_probe"].includes(k));
  host.innerHTML = entries.map(([k, v]) => field(group, k, v)).join("");
}

function getConfig() {
  const cfg = {};
  for (const group of Object.keys(DEFAULTS)) {
    cfg[group] = {};
    for (const [k, v] of Object.entries(DEFAULTS[group])) {
      const el = document.getElementById(`f-${group}-${k}`);
      if (!el) { cfg[group][k] = v; continue; }
      cfg[group][k] = (typeof v === "boolean") ? el.checked : parseFloat(el.value);
    }
  }
  return cfg;
}

// ── selectors + fault catalogue ─────────────────────────────────────────
function buildSelect(id, items) {
  document.getElementById(id).innerHTML = items.map(it =>
    `<option value="${it.id}" ${it.enabled ? "" : "disabled"}>${it.label}${it.enabled ? "" : " · coming soon"}</option>`
  ).join("");
}

function buildFaults() {
  const tabs = document.getElementById("cat-tabs");
  tabs.innerHTML = CATALOGUE.map((c, i) =>
    `<button data-cat="${i}" ${c.enabled ? "" : "disabled"}>${c.category}</button>`).join("");
  tabs.querySelectorAll("button").forEach(b => b.onclick = () => selectCategory(+b.dataset.cat));
  selectCategory(CATALOGUE.findIndex(c => c.enabled));
}

function selectCategory(idx) {
  document.querySelectorAll("#cat-tabs button").forEach(b =>
    b.classList.toggle("active", +b.dataset.cat === idx));
  const cat = CATALOGUE[idx];
  const host = document.getElementById("fault-list");
  if (!cat || !cat.enabled || cat.faults.length === 0) {
    host.innerHTML = `<div class="coming">No configurable faults in this category yet — coming soon.</div>`;
    return;
  }
  host.innerHTML = cat.faults.map(f => {
    if (!f.enabled)
      return `<div class="fault-row disabled"><div class="name">${f.label}</div></div>`;
    const params = (f.params || []).map(p =>
      `<span class="plabel">${p.label}</span>
       <input type="number" step="any" id="ff-${f.id}-${p.key}" value="${p.default}">`).join("");
    return `<div class="fault-row"><div class="name">${f.label}</div>
      <div class="fparams">${params}
      <button class="btn-add" data-fault="${f.id}">＋ inject</button></div></div>`;
  }).join("");
  host.querySelectorAll(".btn-add").forEach(b => b.onclick = () => addFault(b.dataset.fault));
}

function addFault(fid) {
  const cat = CATALOGUE.find(c => c.enabled);
  const spec = cat.faults.find(f => f.id === fid);
  const obj = { fault: fid };
  (spec.params || []).forEach(p => {
    obj[p.key] = parseFloat(document.getElementById(`ff-${fid}-${p.key}`).value);
  });
  if (running && ws && ws.readyState === 1) {
    ws.send(JSON.stringify({ action: "inject_fault", fault: obj }));
  } else {
    pendingFaults.push(obj);
    addChip(`${fid} @ ${obj.t_fault ?? "?"}h (queued)`);
    toast(`Queued ${fid} for next run`);
  }
}

// ── plots ────────────────────────────────────────────────────────────────
const CONC = { "plot-X": ["X", "Biomass X [g/L]", "#7c9ce0"],
               "plot-E": ["E", "Ethanol E [g/L]", "#d98cc0"],
               "plot-G": ["G", "Glucose G [g/L]", "#6fc59b"] };
const ENV  = { "plot-T": ["T", "Temperature [°C]", "#e8a598"],
               "plot-pH": ["pH", "pH", "#b58fe0"],
               "plot-DO": ["DO", "Dissolved O₂ [% sat]", "#5aa9c9"] };

const baseLayout = t => ({
  title: { text: t, font: { size: 13, color: "#3d4257" }, x: 0.02, xanchor: "left" },
  margin: { l: 44, r: 12, t: 30, b: 34 }, height: 230,
  xaxis: { title: { text: "time [h]", font: { size: 10 } }, gridcolor: "#eef0f7" },
  yaxis: { gridcolor: "#eef0f7" }, paper_bgcolor: "white", plot_bgcolor: "white",
  showlegend: false, shapes: [], font: { family: "-apple-system, sans-serif" },
});

function rgba(hex, a) {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}

const HEALTHY_LINE = { color: "#9a9fb2", width: 1.6, dash: "dot" };

// All six panels share the same 5-trace layout so any of them can carry an
// EKF dashed overlay + ±σ band. Only X/G/E have EKF traces (T/pH/DO are
// measured inputs, not estimated, so `frame.ekf` never carries their key and
// pushFrame's presence gate leaves those panels without an overlay).
function initPlots() {
  for (const [id, [, title, color]] of [...Object.entries(CONC), ...Object.entries(ENV)]) {
    Plotly.newPlot(id, [
      { x: [], y: [], mode: "lines", line: { width: 0 }, hoverinfo: "skip" },              // 0 band lo
      { x: [], y: [], mode: "lines", fill: "tonexty", fillcolor: rgba(color, 0.15),
        line: { width: 0 }, hoverinfo: "skip" },                                            // 1 band hi
      { x: [], y: [], mode: "lines", line: { color: color, width: 2.2 }, name: "true" },    // 2 true
      { x: [], y: [], mode: "lines", line: { color: "#7a7f95", width: 1.6, dash: "dash" }, name: "EKF" }, // 3 ekf
      { x: [], y: [], mode: "lines", line: HEALTHY_LINE, name: "no fault", visible: false },// 4 healthy
    ], baseLayout(title), { displayModeBar: false, responsive: true });
  }
  faultLines = [];
}

// ── detector panels (Design A only) ─────────────────────────────────────
// Left: the raw innovation ν = z_E − Ê, with the ±0.1 g/L band shaded — inside
// it the cultivation is running as expected, outside it is interrupted.
// Right: the two CUSUM arms accumulating evidence against the threshold h.
let nuShapes = [], cusumShapes = [], alarmDrawn = false;

function initDetector() {
  const hEl = document.getElementById("f-ekf-cusum_h");
  const h = hEl ? parseFloat(hEl.value) : 1.5;

  // ±0.1 g/L: inside this band the cultivation is running as expected.
  nuShapes = [{ type: "rect", xref: "paper", x0: 0, x1: 1, y0: -0.1, y1: 0.1,
    fillcolor: "rgba(122,127,149,0.10)", line: { width: 0 }, layer: "below" }];
  Plotly.newPlot("plot-nu", [
    { x: [], y: [], mode: "lines", line: { color: "#2a8a5a", width: 1.8 }, name: "ν" },
  ], Object.assign(baseLayout("Innovation ν = z_E − Ê  [g/L]"), { shapes: nuShapes }),
    { displayModeBar: false, responsive: true });

  cusumShapes = [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: h, y1: h,
    line: { color: "#7a7f95", width: 1.2, dash: "dash" } }];
  Plotly.newPlot("plot-cusum", [
    { x: [], y: [], mode: "lines", line: { color: "#b8862b", width: 1.8 }, name: "S⁺" },
    { x: [], y: [], mode: "lines", line: { color: "#9a4b30", width: 2.0 }, name: "S⁻" },
  ], Object.assign(baseLayout(`CUSUM arms  [√S] — h = ${h}`), { shapes: cusumShapes,
      showlegend: true,
      legend: { x: 0.02, y: 0.98, font: { size: 10 }, bgcolor: "rgba(255,255,255,0.7)" } }),
    { displayModeBar: false, responsive: true });

  alarmDrawn = false;
  setDetectorVisible(true);
  setAlarmBadge(null);
}

function setDetectorVisible(on) {
  document.getElementById("detector").style.display = on ? "" : "none";
}

function setAlarmBadge(cs) {
  const el = document.getElementById("cusum-state");
  if (!el) return;
  if (!cs) { el.textContent = "—"; el.className = "pill"; return; }
  if (cs.alarm_t != null) {
    el.textContent = `ALARM · t = ${cs.alarm_t.toFixed(2)} h`;
    el.className = "pill alarm";
  } else if (cs.blind) {
    // The healthy reference itself makes less ethanol than the sensor noise, so
    // a quiet CUSUM means "cannot see", not "nothing wrong". WHICH modifier is
    // responsible has to be read off the numbers — this used to blame the
    // cardinal T window unconditionally, which sends the user to the wrong knob
    // whenever pH is the one that is off.
    const fT = cs.fT_ss, fpH = cs.fpH_set;
    const culprit = (fpH != null && fT != null && fpH < fT) ? "pH" : "T";
    el.textContent = `BLIND · healthy peak E ${cs.hE_peak.toFixed(3)} g/L < 3σ`;
    el.className = "pill warn";
    el.title =
      `The HEALTHY reference itself makes almost no ethanol, so a quiet CUSUM `
    + `means "cannot see", not "nothing wrong".\n\n`
    + `Where it is being lost: f_T = ${fT != null ? fT.toFixed(2) : "?"} at the `
    + `${cs.T_ss != null ? cs.T_ss.toFixed(1) : "~27.9"} °C the vessel actually `
    + `holds, f_pH = ${fpH != null ? fpH.toFixed(2) : "?"} at the pH the `
    + `controller holds. Growth is the PRODUCT of the two, so the smaller one is `
    + `the binding constraint — here that is ${culprit}.\n\n`
    + (culprit === "T"
        ? `A proportional heater needs a standing error to make power, so the `
        + `vessel sits ~2 °C BELOW setpoint all batch — cardinals chosen against `
        + `T_set rather than against the temperature actually reached silently `
        + `cripple the culture. Fix: widen the cardinal T window until f_T there `
        + `is well above 0.5. BOTH arms matter — (27, 30, 32) gives 0.28 and `
        + `never alarms; (27, 30, 33) gives 0.50 and alarms at 1.17 h.`
        : `The pH modifier is normalised AT the setpoint, so it reads 1 there by `
        + `definition — the number above is the un-anchored one, relative to the `
        + `organism's optimum. A setpoint parked in a dead corner of the window `
        + `leaves the healthy reference crippled, and then a fault that moves pH `
        + `back TOWARDS pH_opt makes the culture grow BETTER than "healthy" and `
        + `the detector reads it backwards. Fix: move pH_set towards pH_opt, or `
        + `widen the cardinal pH window.`)
    + `\n\nA culture growing that slowly respires everything it takes up and `
    + `ferments nothing. This blinds T and pH faults only — oxygen faults still `
    + `alarm, because they redirect carbon into fermentation instead of slowing it.`;
  } else {
    el.textContent = `quiet · ${Math.max(cs.Sp, cs.Sn).toFixed(2)} / ${cs.h}`;
    el.className = "pill ok";
  }
}

function pushDetector(fr) {
  const cs = fr.cusum, t = fr.t_h;
  if (!cs) return;
  Plotly.extendTraces("plot-nu", { x: [[t]], y: [[cs.nu]] }, [0]);
  // one x array PER trace — Plotly throws if the lengths disagree
  Plotly.extendTraces("plot-cusum", { x: [[t], [t]], y: [[cs.Sp], [cs.Sn]] }, [0, 1]);
  setAlarmBadge(cs);
  if (cs.alarm_t != null && !alarmDrawn) {
    alarmDrawn = true;
    const line = { type: "line", x0: cs.alarm_t, x1: cs.alarm_t, yref: "paper",
      y0: 0, y1: 1, line: { color: "#9a4b30", width: 1.4 } };
    nuShapes.push(line); cusumShapes.push(line);
    Plotly.relayout("plot-nu", { shapes: nuShapes });
    Plotly.relayout("plot-cusum", { shapes: cusumShapes });
  }
}

// ── root cause: the MMAE bank ───────────────────────────────────────────
// The detector answers "is something wrong". This answers "which fault", and
// then shows its work: the posterior as it built up, the innovation each
// hypothesis had to explain, and what each one predicted on channels the
// ethanol-only bank never got to see.
const RCA_COLORS = ["#7c9ce0", "#d98cc0", "#6fc59b", "#e8a598", "#b58fe0",
                    "#5aa9c9", "#c9a227"];
let lastMMAE = null, rcaMode = "ethanol", rcaChannel = "T";
// A verdict taken at the alarm is not wrong so much as PREMATURE: on ethanol
// alone only half the causes are right at that moment, and the rest change
// their mind over the next 35 min on average (max 90). Naming one on screen
// and quietly replacing it later reads as the method being unreliable, when
// what it actually is, is honest about how much batch it has seen. So the
// name is withheld until the evidence window the bank itself declares has
// passed. The operator can still force it — nothing is hidden, only unshown.
let revealProvisional = false;

function rcaColor(res, id) {
  const i = res.causes.findIndex(c => c.id === id);
  return id === "healthy" ? "#9aa1b5" : RCA_COLORS[i % RCA_COLORS.length];
}

function _provisionalNote(sc) {
  // C.4: an ethanol-only verdict taken at the alarm has not settled yet.
  if (!sc || !sc.provisional) return "";
  const m = Math.round(sc.evidence_min_after_alarm || 0);
  const probes = rcaMode === "probes";
  return `<div class="note warn">PROVISIONAL — this verdict stands on `
       + `${m} min of batch after the alarm. `
       + (probes
          ? `With the probes the bank usually settles within 30 min; the two `
          + `heater rungs are the slow pair, because at the alarm neither has `
          + `finished cooling. Re-diagnose shortly.`
          : `On ethanol alone it can take up to 120 min more before it stops `
          + `changing, and the acid/base pump pair never separates. Turn on `
          + `<b>RCA reads probes</b>, or re-diagnose later.`)
       + `</div>`;
}

function setMmaeBadge(text, cls) {
  const el = document.getElementById("mmae-state");
  if (el) { el.textContent = text; el.className = "pill" + (cls ? " " + cls : ""); }
  const st = document.getElementById("rca-status");
  if (st) { st.textContent = text; st.className = "pill" + (cls ? " " + cls : ""); }
}

function showRca(on) { document.getElementById("rca").hidden = !on; }

// Log-odds run to 10^5 once the probes are scoring: the plant model is
// deterministic, so a mechanism that matches beats one that does not by more
// than any real instrument justifies. Report the size, not the digits.
function fmtMargin(m) {
  if (m == null) return "—";
  if (m > 1000) return "≫ 10³";
  return m.toFixed(1);
}

// Distance from the winner, which is the only part of a log-likelihood that
// means anything here — absolute values run to 10⁵ once the probes are scored.
function fmtDelta(d) {
  if (d === 0) return "best";
  const a = Math.abs(d);
  return "−" + (a >= 1e4 ? a.toExponential(1).replace("e+", "e") : a.toFixed(1));
}

function renderRca() {
  const res = lastMMAE;
  if (!res || res.state !== "ok") return;
  const sc = res.scores[rcaMode];
  const byId = Object.fromEntries(res.causes.map(c => [c.id, c]));
  const best = sc.loglik[sc.verdict];
  document.querySelectorAll("#rca-modes button").forEach(b =>
    b.classList.toggle("active", b.dataset.mode === rcaMode));

  // ── hold the name until the bank has the batch it says it needs ───────
  const need = rcaMode === "probes" ? res.probe_settle_min : res.ethanol_settle_min;
  const have = sc.evidence_min_after_alarm || 0;
  if (sc.provisional && !revealProvisional) {
    const pct = Math.min(100, 100 * have / Math.max(need, 1e-9));
    setMmaeBadge(`forming · ${Math.round(have)} / ${Math.round(need)} min`, "warn");
    document.getElementById("rca-verdict").innerHTML =
      `<span class="tag null">forming</span>
       <b>Not yet.</b> The bank has <b>${Math.round(have)} min</b> of batch after the
       alarm and wants <b>${Math.round(need)}</b> before it will put a name to this.
       <div class="rank" style="margin-top:10px">
         <div class="nm">evidence since the alarm</div>
         <div class="bar"><i style="width:${pct.toFixed(1)}%"></i></div>
         <div class="num">${Math.round(have)} / ${Math.round(need)} min</div>
       </div>
       <div class="note">` +
      (rcaMode === "probes"
        ? `With the T/pH/DO probes scoring, this is usually short — the causes differ
           in <i>which</i> channel moved and in which direction, which is one sample's
           worth of evidence. The two heater rungs are the slow pair: at the alarm
           neither has finished cooling, so they still look alike.`
        : `On ethanol alone the alarm fires as soon as the innovation is improbable,
           which is well before the drift has a recognisable <i>shape</i>. Measured on
           this configuration: half the causes are already right at the alarm, the rest
           need 60–90 min more, and the acid/base pump pair never separates at all.
           Switching on <b>RCA reads probes</b> answers it now instead.`) +
      ` It will name itself when the window closes — no need to re-diagnose.</div>
       <button class="btn btn-ghost" id="btn-reveal" style="margin-top:8px">
         Show the provisional answer anyway</button>`;
    document.getElementById("rca-ranks").innerHTML = "";
    document.getElementById("rca-note").innerHTML = "";
    const rb = document.getElementById("btn-reveal");
    if (rb) rb.onclick = () => { revealProvisional = true; renderRca(); };
    return;
  }

  setMmaeBadge((sc.tied.length > 1 ? `tie · ${sc.tied.length} causes`
                                  : sc.verdict_label)
               + (sc.provisional ? " · provisional" : ""), "alarm");

  // ── the verdict line ──────────────────────────────────────────────────
  const onset = sc.onset[sc.verdict];
  let tag, text;
  if (sc.verdict === "healthy") {
    tag = `<span class="tag null">no cause fits</span>`;
    text = `The bank explains the alarm best with <b>no fault at all</b> — no
      candidate mechanism, at the severity it assumes, matches this batch.`;
  } else if (sc.tied.length > 1) {
    tag = `<span class="tag tie">tie</span>`;
    text = `<b>${sc.tied.map(id => byId[id].label).join(" / ")}</b> —
      indistinguishable here (Δ log-lik = ${fmtMargin(sc.margin_nats)} nats,
      under the 2.3 that would be 10:1 odds). An honest tie, not a pick.`;
  } else {
    tag = `<span class="tag decisive">decisive</span>`;
    text = `<b>${byId[sc.verdict].label}</b>, starting near
      <b>${onset == null ? "—" : onset.toFixed(2) + " h"}</b> —
      ahead of the runner-up by ${fmtMargin(sc.margin_nats)} nats.`;
  }
  document.getElementById("rca-verdict").innerHTML =
    `${tag}${text}${_provisionalNote(sc)}`;

  // ── the ranking ───────────────────────────────────────────────────────
  document.getElementById("rca-ranks").innerHTML = sc.order.map(id => {
    const c = byId[id], p = sc.p[id];
    const flags = [];
    if (c.silent) flags.push(`<span class="flag" title="never separates from 'no fault' on ethanol (${c.sep.toFixed(1)} √S) — it cannot be isolated from the soft sensor at any threshold">silent on ethanol</span>`);
    if (sc.onset[id] != null) flags.push(
      `<span class="flag">onset ${sc.onset[id].toFixed(2)} h</span>`);
    return `<div class="rank ${id === sc.verdict ? "win" : ""}">
      <div class="nm">${c.label}${flags.join("")}</div>
      <div class="bar"><i style="width:${(100 * p).toFixed(1)}%"></i></div>
      <div class="num" title="log-likelihood ${sc.loglik[id].toFixed(1)}">p = ${p.toFixed(3)} · Δ ${fmtDelta(sc.loglik[id] - best)}</div>
    </div>`;
  }).join("");

  const notes = [];
  notes.push(`Onset walked back from the alarm at ${res.alarm_t.toFixed(2)} h to
    t̂₀ = ${res.t_onset.toFixed(2)} h; each cause also tried
    ${res.onsets.map(o => o.toFixed(2)).join(", ")} h and kept its best
    (profiled, not summed). Scored on ${res.n_samples} ethanol samples up to
    ${res.t_evidence.toFixed(2)} h.`);
  if (sc.degenerate.length)
    notes.push(`Provably inseparable on this batch: ` + sc.degenerate
      .map(([a, b]) => `${byId[a].label} ≡ ${byId[b].label}`).join(", ") +
      ` — identical scores, so no estimator can split them.`);
  if (rcaMode === "ethanol")
    notes.push(`Ethanol only — the soft-sensor question. Switch to the probes to
      score the T/pH/DO instruments the rig already has.`);
  else
    notes.push(`T/pH/DO probes included in the likelihood, at σ =
      ${DEFAULTS.ekf.sigma_T_probe} °C / ${DEFAULTS.ekf.sigma_pH_probe} /
      ${DEFAULTS.ekf.sigma_DO_probe} % sat. Each hypothesis still commits to one
      severity, so this identifies the mechanism, not how hard it hit.`);
  document.getElementById("rca-note").innerHTML = notes.join("<br>");

  rcaPlots(res, sc, byId);
}

function rcaPlots(res, sc, byId) {
  const tm = res.t_meas;
  const onsetLine = {
    type: "line", x0: res.t_onset, x1: res.t_onset, yref: "paper", y0: 0, y1: 1,
    line: { color: "#9a4b30", width: 1.2 },
  };
  const lay = (title, extra) => Object.assign(baseLayout(title), {
    shapes: [onsetLine], showlegend: true,
    legend: { x: 0.02, y: 0.98, font: { size: 9 }, bgcolor: "rgba(255,255,255,.72)" },
  }, extra || {});

  // (a) the posterior as the evidence came in — column k is the belief after
  // sample k, so a curve that climbs and stays is a cause the data kept backing.
  Plotly.react("plot-post", res.causes.map(c => ({
    x: tm, y: sc.post[c.id], mode: "lines", name: c.label,
    line: { color: rcaColor(res, c.id), width: c.id === sc.verdict ? 2.6 : 1.3 },
  })), lay("Posterior over causes", { yaxis: { range: [-0.03, 1.03], gridcolor: "#eef0f7" } }),
    { displayModeBar: false, responsive: true });

  // (b) what step 5 actually scores: small, unbiased ν = the right model.
  const top = sc.order.slice(0, 4);
  Plotly.react("plot-hnu", top.map(id => ({
    x: tm, y: byId[id].nu, mode: "lines", name: byId[id].label,
    line: { color: rcaColor(res, id), width: 1.4 },
  })), lay("ν per hypothesis  [g/L] — small = right model"),
    { displayModeBar: false, responsive: true });

  // (c) the one channel we measure, and what each hypothesis expected of it.
  Plotly.react("plot-hE", [
    { x: tm, y: res.z_meas, mode: "markers", name: "measured",
      marker: { size: 4, color: "#7a7f95" } },
    ...top.map(id => ({ x: tm, y: byId[id].E, mode: "lines", name: byId[id].label,
      line: { color: rcaColor(res, id), width: 1.4, dash: "dash" } })),
  ], lay("Ethanol — measured vs each hypothesis"),
    { displayModeBar: false, responsive: true });

  // (d) the channels ethanol cannot see. With the probes off this panel is the
  // visual check on a verdict reached without them; with the probes on it is
  // the evidence that broke the tie.
  const CH = { T: "Temperature [°C]", pH: "pH", DO: "Dissolved O₂ [% sat]" };
  Plotly.react("plot-hch", [
    { x: tm, y: res.meas[rcaChannel], mode: "lines", name: "probe",
      line: { color: "#3d4257", width: 2.2 } },
    ...top.map(id => ({ x: tm, y: byId[id][rcaChannel], mode: "lines",
      name: byId[id].label,
      line: { color: rcaColor(res, id), width: 1.4, dash: "dash" } })),
  ], lay(`${CH[rcaChannel]} — probe vs each hypothesis`),
    { displayModeBar: false, responsive: true });
}

function onMmae(m) {
  if (m.state !== "ok") {
    setMmaeBadge(m.state === "no_alarm" ? "no alarm" : "n/a");
    return;
  }
  lastMMAE = m;
  showRca(true);
  document.getElementById("btn-diagnose").disabled = !running;
  renderRca();
}

function resetRca() {
  lastMMAE = null;
  showRca(false);
  setMmaeBadge("—");
  document.getElementById("rca-verdict").innerHTML = "";
  document.getElementById("rca-ranks").innerHTML = "";
  document.getElementById("rca-note").innerHTML = "";
  document.getElementById("btn-diagnose").disabled = true;
}

function pushFrame(fr) {
  const t = fr.t_h, tr = fr.true, hh = fr.healthy;
  // running estimator error on the hidden state — the number to watch while
  // sweeping σ sensor to place a mild/moderate/severe line
  const out = document.getElementById("nrmse-X");
  if (out) out.textContent = (fr.nrmse_X == null) ? "—" : `${fr.nrmse_X.toFixed(1)} %`;
  for (const [id, [key]] of [...Object.entries(CONC), ...Object.entries(ENV)]) {
    Plotly.extendTraces(id, { x: [[t]], y: [[tr[key]]] }, [2]);                             // true
    if (hh) Plotly.extendTraces(id, { x: [[t]], y: [[hh[key]]] }, [4]);                     // healthy
    // EKF estimate: every panel in Design A (T/pH/DO are estimated states of
    // the healthy shadow), X/G/E only in Design B.
    if (fr.ekf && fr.ekf[key] !== undefined) {
      Plotly.extendTraces(id, { x: [[t]], y: [[fr.ekf[key]]] }, [3]);
      if (bandOn() && fr.ekf["s" + key] !== undefined) {
        const s = fr.ekf["s" + key];
        Plotly.extendTraces(id, { x: [[t], [t]], y: [[fr.ekf[key] - s], [fr.ekf[key] + s]] }, [0, 1]);
      }
    }
  }
  // Detector LAST, so a problem here can never blank the science panels.
  // Design B has no healthy shadow and so nothing worth detecting on: the
  // backend sends cusum: null and the row goes away.
  if (fr.cusum) pushDetector(fr); else setDetectorVisible(false);
}

function bandOn() {
  const el = document.getElementById("f-ekf-band");
  return el && el.checked;
}

function setCompareVisibility(on) {
  for (const id of [...Object.keys(CONC), ...Object.keys(ENV)])
    Plotly.restyle(id, { visible: on }, [4]);
}

function addFaultLine(t) {
  faultLines.push(t);
  const shapes = faultLines.map(x => ({ type: "line", x0: x, x1: x, yref: "paper", y0: 0, y1: 1,
    line: { color: "#9a4b30", width: 1, dash: "dot" } }));
  for (const id of [...Object.keys(CONC), ...Object.keys(ENV)])
    Plotly.relayout(id, { shapes });
}

// ── chips / toast ────────────────────────────────────────────────────────
function addChip(txt) {
  const c = document.createElement("span");
  c.className = "chip"; c.textContent = txt;
  document.getElementById("fault-chips").appendChild(c);
}
let toastTimer = null;
function toast(msg, err = false) {
  const el = document.getElementById("toast");
  el.textContent = msg; el.className = "show" + (err ? " err" : "");
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => (el.className = ""), 2600);
}

// ── run control ──────────────────────────────────────────────────────────
let lastCSV = null, lastPNG = null;

function run() {
  if (running) return;
  initPlots();
  initDetector();
  setCompareVisibility(document.getElementById("cmp-healthy").checked);
  document.getElementById("fault-chips").innerHTML = "";
  lastCSV = lastPNG = null;
  document.getElementById("btn-csv").disabled = true;
  document.getElementById("btn-png").disabled = true;
  resetRca();

  ws = new WebSocket(`ws://${location.host}/ws/run`);
  ws.onopen = () => ws.send(JSON.stringify({ action: "start", config: getConfig() }));
  ws.onmessage = ev => onMessage(JSON.parse(ev.data));
  ws.onerror = () => toast("WebSocket error", true);
  ws.onclose = () => { if (running) endRun(); };
}

function onMessage(m) {
  switch (m.type) {
    case "config_error":
      // The plant refused to build: a setpoint outside the organism's viable
      // window. Say so loudly — the whole point of the guard is that this used
      // to run silently and produce a dead healthy reference.
      toast(m.note, true);
      endRun();
      break;
    case "started":
      revealProvisional = false;
      running = true;
      document.getElementById("btn-run").disabled = true;
      document.getElementById("btn-stop").disabled = false;
      document.getElementById("tfinal").textContent = m.t_final.toFixed(1);
      // fire any pre-queued faults
      pendingFaults.forEach(f => ws.send(JSON.stringify({ action: "inject_fault", fault: f })));
      pendingFaults = [];
      document.getElementById("fault-chips").innerHTML = "";
      break;
    case "frame":
      pushFrame(m);
      document.getElementById("clock").textContent = m.t_h.toFixed(2);
      break;
    case "fault_ack":
      if (m.ok) {
        toast(m.note);
        const tf = parseFloat((m.note.match(/t = ([\d.]+)/) || [])[1]);
        if (!isNaN(tf)) { addFaultLine(tf); addChip(m.note.split(" scheduled")[0]); }
      } else {
        toast(m.note, true);
      }
      break;
    case "mmae_status":
      if (m.state === "running") {
        showRca(true);
        setMmaeBadge(`diagnosing… (${m.t_evidence.toFixed(2)} h)`);
        document.getElementById("btn-diagnose").disabled = true;
      } else if (m.state === "busy") {
        toast("A diagnosis is already running");
      }
      break;
    case "mmae":
      onMmae(m);
      break;
    case "complete":
      lastCSV = m.csv; lastPNG = m.png;
      document.getElementById("btn-csv").disabled = false;
      document.getElementById("btn-png").disabled = false;
      toast(m.stopped ? "Run stopped — data ready to download" : "Run complete ✓");
      endRun();
      break;
  }
}

function endRun() {
  running = false;
  document.getElementById("btn-run").disabled = false;
  document.getElementById("btn-stop").disabled = true;
  // the bank replays recorded samples on the server, so it needs a live socket
  document.getElementById("btn-diagnose").disabled = true;
  if (ws && ws.readyState === 1) ws.close();
}

function stop() {
  if (ws && ws.readyState === 1) ws.send(JSON.stringify({ action: "stop" }));
}

function reset() {
  if (running) stop();
  if (ws && ws.readyState === 1) ws.close();
  running = false;
  // restore every input to its default
  ["ekf", "initial", "kinetics", "thermal", "ph", "oxygen", "sim"].forEach(buildGroup);
  // clear plots, chips, queued faults, downloads, clock
  document.getElementById("cmp-healthy").checked = true;
  initPlots();
  initDetector();
  resetRca();
  setCompareVisibility(true);
  pendingFaults = [];
  document.getElementById("fault-chips").innerHTML = "";
  document.getElementById("clock").textContent = "0.00";
  document.getElementById("tfinal").textContent = DEFAULTS.sim.t_final.toFixed(1);
  lastCSV = lastPNG = null;
  document.getElementById("btn-run").disabled = false;
  document.getElementById("btn-stop").disabled = true;
  document.getElementById("btn-csv").disabled = true;
  document.getElementById("btn-png").disabled = true;
  toast("Reset to defaults");
}

// ── downloads ────────────────────────────────────────────────────────────
function download(name, href) {
  const a = document.createElement("a");
  a.href = href; a.download = name; a.click();
}

// ── boot ─────────────────────────────────────────────────────────────────
async function boot() {
  const data = await (await fetch("/api/defaults")).json();
  DEFAULTS = data.params; CATALOGUE = data.faults;
  ["ekf", "initial", "kinetics", "thermal", "ph", "oxygen", "sim"].forEach(buildGroup);
  buildSelect("sel-cell", data.cells);
  buildSelect("sel-reactor", data.reactors);
  buildFaults();
  initPlots();
  initDetector();
  setCompareVisibility(document.getElementById("cmp-healthy").checked);
  document.getElementById("tfinal").textContent = DEFAULTS.sim.t_final.toFixed(1);

  // accordion toggles
  document.querySelectorAll(".acc-head").forEach(h =>
    h.onclick = () => h.parentElement.classList.toggle("collapsed"));

  document.getElementById("btn-run").onclick = run;
  document.getElementById("btn-stop").onclick = stop;
  document.getElementById("btn-reset").onclick = reset;

  // Root-cause controls. The mode tabs re-render from the SAME result — both
  // scorings come back from one bank run, so switching costs nothing and the
  // two answers are always over identical evidence.
  rcaMode = DEFAULTS.ekf.mmae_probes ? "probes" : "ethanol";
  document.querySelectorAll("#rca-modes button").forEach(b =>
    b.onclick = () => { rcaMode = b.dataset.mode; renderRca(); });
  document.getElementById("rca-channel").onchange = e => {
    rcaChannel = e.target.value; renderRca();
  };
  document.getElementById("btn-diagnose").onclick = () => {
    if (ws && ws.readyState === 1) ws.send(JSON.stringify({ action: "diagnose" }));
  };
  resetRca();
  document.getElementById("cmp-healthy").onchange = e => setCompareVisibility(e.target.checked);
  document.getElementById("btn-csv").onclick = () =>
    lastCSV && download("bioreactor_run.csv", "data:text/csv;charset=utf-8," + encodeURIComponent(lastCSV));
  document.getElementById("btn-png").onclick = () =>
    lastPNG && download("bioreactor_run.png", "data:image/png;base64," + lastPNG);
}
boot();
