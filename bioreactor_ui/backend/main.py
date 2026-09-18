"""
main.py — FastAPI app: static UI, defaults, WebSocket run loop, CSV/PNG export.

Run:  uvicorn backend.main:app  (from the bioreactor_ui/ directory)
Then open http://localhost:8000
"""
from __future__ import annotations

import asyncio
import base64
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .core import DEFAULTS
from .faults import FAULT_CATALOGUE
from .mmae import ETHANOL_SETTLE_MIN, PROBE_SETTLE_MIN, run_mmae
from .simulator import Simulation

FRONTEND = Path(__file__).resolve().parent.parent / "frontend"

app = FastAPI(title="Bioreactor Digital-Shadow Simulator")

# Cell / bioreactor selectors (only one enabled each for now).
CELLS = [
    {"id": "yeast", "label": "Yeast (S. cerevisiae)", "enabled": True},
    {"id": "cho",   "label": "CHO (mammalian)",       "enabled": False},
    {"id": "ecoli", "label": "E. coli",                "enabled": False},
]
REACTORS = [
    {"id": "minifors2", "label": "Minifors 2 (1.35 L)", "enabled": True},
    {"id": "biostat",   "label": "BIOSTAT (coming soon)", "enabled": False},
]


@app.get("/api/defaults")
def api_defaults():
    return {
        "params": DEFAULTS,
        "faults": FAULT_CATALOGUE,
        "cells": CELLS,
        "reactors": REACTORS,
    }


# ── CSV / PNG builders ────────────────────────────────────────────────────
CSV_COLS = ["time_h", "T", "pH", "DO", "G", "E", "X", "mu_G", "mu_E"]
# Design-A extras: the estimate and the detector trace, keyed by time so a row
# without them (EKF off) still lines up.
EKF_COLS = ["X", "G", "E"]
DET_COLS = ["nu", "Sp", "Sn"]


def build_csv(sim: Simulation) -> str:
    ek = {r["time_h"]: r for r in sim.ekf_rows}
    det = {r["time_h"]: r for r in sim.det_rows}
    cols = list(CSV_COLS)
    if ek:
        cols += [f"ekf_{c}" for c in EKF_COLS]
    if det:
        cols += DET_COLS
    lines = [",".join(cols)]
    for r in sim.rows:
        vals = [f"{r[c]:.6g}" for c in CSV_COLS]
        t = r["time_h"]
        if ek:
            e = ek.get(t)
            vals += [f"{e[c]:.6g}" if e else "" for c in EKF_COLS]
        if det:
            d = det.get(t)
            vals += [f"{d[c]:.6g}" if d else "" for c in DET_COLS]
        lines.append(",".join(vals))
    return "\n".join(lines) + "\n"


def _rca_panels(row, rca: dict):
    """The root-cause row: the two verdicts, and the evidence behind them.

    Left is the answer, middle and right are the proof — the innovation each
    hypothesis had to explain (what the ethanol-only bank scores) and the
    temperature each predicted against what the probe read (what breaks the
    tie when the ethanol scores are indistinguishable).
    """
    ax_p, ax_nu, ax_T = row
    causes = [c for c in rca["causes"]]
    labels = [c["label"] for c in causes]
    y = range(len(causes))
    pe = [rca["scores"]["ethanol"]["p"][c["id"]] for c in causes]
    pp = [rca["scores"]["probes"]["p"][c["id"]] for c in causes]
    ax_p.barh([v + 0.2 for v in y], pe, height=0.38, color="#7c9ce0",
              label="ethanol only")
    ax_p.barh([v - 0.2 for v in y], pp, height=0.38, color="#e0a3c8",
              label="+ T/pH/DO probes")
    ax_p.set_yticks(list(y)); ax_p.set_yticklabels(labels, fontsize=8)
    ax_p.set_xlim(0, 1); ax_p.legend(fontsize=7)
    ax_p.set_title("Posterior over causes  (onset {:.2f} h)".format(rca["t_onset"]))

    tm = rca["t_meas"]
    top = rca["scores"]["ethanol"]["order"][:4]
    for cid in top:
        c = next(x for x in causes if x["id"] == cid)
        ax_nu.plot(tm, c["nu"], lw=1.2, label=c["label"])
        ax_T.plot(tm, c["T"], "--", lw=1.1, label=c["label"])
    ax_nu.axhline(0, color="k", lw=0.6); ax_nu.legend(fontsize=7)
    ax_nu.set_title("ν per hypothesis — small = right model")
    ax_T.plot(tm, rca["meas"]["T"], "k-", lw=1.6, label="probe")
    ax_T.legend(fontsize=7)
    ax_T.set_title("Predicted T vs the probe [°C]")
    for ax in (ax_nu, ax_T):
        ax.axvline(rca["t_onset"], color="#9a4b30", lw=1.2)
        ax.set_xlabel("time [h]"); ax.grid(alpha=0.3)


def build_png(sim: Simulation) -> str:
    rows = sim.rows
    t = [r["time_h"] for r in rows]
    panels = [
        ("Biomass X [g/L]", "X", "#7c9ce0"),
        ("Ethanol E [g/L]", "E", "#e0a3c8"),
        ("Glucose G [g/L]", "G", "#8fd6b4"),
        ("Temperature [°C]", "T", "#e8a598"),
        ("pH", "pH", "#c9a3e0"),
        ("Dissolved O₂ [% sat]", "DO", "#5aa9c9"),
    ]
    ek = {r["time_h"]: r for r in sim.ekf_rows}
    ek_t = [r["time_h"] for r in sim.ekf_rows]
    # Overlay the no-fault baseline only when a fault was actually injected.
    show_healthy = bool(sim.active_faults)
    det = sim.det_rows                      # non-empty only in Design A

    rca = sim.mmae if (sim.mmae or {}).get("state") == "ok" else None
    nrow = (3 if det else 2) + (1 if rca else 0)
    fig, axes = plt.subplots(nrow, 3, figsize=(15, 4 * nrow))
    for ax, (title, key, color) in zip(axes.flat, panels):
        legend = False
        if show_healthy:
            ax.plot(t, [r[key] for r in sim.healthy_rows], ":", color="#9a9fb2",
                    lw=1.5, label="no fault")
            legend = True
        ax.plot(t, [r[key] for r in rows], color=color, lw=1.8, label="true")
        # EKF dashed overlay. In Design A every panel has one (T/pH/DO are
        # estimated states of the healthy shadow); in Design B only X/E/G do,
        # because there T/pH/DO are measured inputs and carry no estimate.
        if sim.ekf_on and ek_t and key in ek[ek_t[0]]:
            ax.plot(ek_t, [ek[tt][key] for tt in ek_t], "--", color="#555",
                    lw=1.2, label="EKF")
            legend = True
        if legend:
            ax.legend(fontsize=8)
        for f in sim.active_faults:
            ax.axvline(f["t_fault"], color="k", ls=":", alpha=0.5)
        ax.set_title(title); ax.set_xlabel("time [h]"); ax.grid(alpha=0.3)

    if det:
        dt_ = [r["time_h"] for r in det]
        alarm = det[-1]["alarm_t"]
        ax = axes[2][0]
        ax.axhspan(-0.1, 0.1, color="#7a7f95", alpha=0.12)
        ax.plot(dt_, [r["nu"] for r in det], color="#2a8a5a", lw=1.4)
        ax.set_title("Innovation ν = z_E − Ê [g/L]")
        ax = axes[2][1]
        ax.plot(dt_, [r["Sp"] for r in det], color="#b8862b", lw=1.4, label="S⁺")
        ax.plot(dt_, [r["Sn"] for r in det], color="#9a4b30", lw=1.6, label="S⁻")
        ax.axhline(det[-1]["h"], color="#7a7f95", ls="--", lw=1.1,
                   label=f"h = {det[-1]['h']}")
        ax.legend(fontsize=8)
        ax.set_title("CUSUM arms [√S]")
        for ax in (axes[2][0], axes[2][1]):
            if alarm is not None:
                ax.axvline(alarm, color="#9a4b30", lw=1.3)
            for f in sim.active_faults:
                ax.axvline(f["t_fault"], color="k", ls=":", alpha=0.5)
            ax.set_xlabel("time [h]"); ax.grid(alpha=0.3)
        axes[2][2].axis("off")
        axes[2][2].text(0.02, 0.6,
                        "alarm: " + (f"t = {alarm:.2f} h" if alarm is not None
                                     else "none"),
                        fontsize=12,
                        color="#9a4b30" if alarm is not None else "#2a7a52")

    if rca:
        _rca_panels(axes[nrow - 1], rca)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# ── WebSocket run loop ────────────────────────────────────────────────────
@app.websocket("/ws/run")
async def ws_run(ws: WebSocket):
    await ws.accept()
    inbox: asyncio.Queue = asyncio.Queue()

    async def receiver():
        try:
            while True:
                await inbox.put(await ws.receive_json())
        except WebSocketDisconnect:
            await inbox.put({"action": "__disconnect__"})

    recv_task = asyncio.create_task(receiver())
    sim: Simulation | None = None
    try:
        while True:
            msg = await inbox.get()
            if msg.get("action") == "__disconnect__":
                break
            if msg.get("action") != "start":
                continue

            # A config whose setpoint sits outside the organism's viable window is
            # rejected by plant.build_modifiers rather than run silently (see the
            # anchor guard in core.make_cardinal_modifier). Catch it here: an
            # unhandled raise would take the socket down with no explanation, which
            # is the same silence the guard exists to end.
            try:
                sim = Simulation(msg.get("config", {}))
            except ValueError as exc:
                await ws.send_json({"type": "config_error", "note": str(exc)})
                sim = None
                continue
            stop = False
            frame_wall = sim.params["sim"]["sec_per_hour"] * sim.dt_h  # sec/frame
            await ws.send_json({"type": "started", "t_final": sim.t_final})

            # Root-cause isolation runs OFF the event loop: a bank of ~19 filter
            # replays takes seconds, and the run must keep streaming while it
            # thinks. One diagnosis at a time; a request arriving while one is
            # in flight is dropped rather than queued, since the newer evidence
            # will be picked up by the next explicit re-diagnosis anyway.
            diag: asyncio.Task | None = None
            # Sim-times [h] at which to re-run the bank without being asked.
            #
            # The bank fires once at the alarm because it is cheap and the panel
            # needs something to show, but on ethanol alone that verdict has NOT
            # settled: measured mean 35 min after the alarm, max 90
            # (8_evaluation.ipynb). The UI withholds the cause NAME until the
            # window the bank itself declares has passed — so the run has to
            # produce a fresh verdict at that moment on its own. Leaving it to a
            # "re-diagnose" click asks the operator to know when the evidence
            # matured, which is the one thing they cannot know.
            redo_at: list[float] = []

            async def diagnose():
                await ws.send_json({"type": "mmae_status", "state": "running",
                                    "t_evidence": sim.t})
                res = await asyncio.to_thread(run_mmae, sim)
                sim.mmae = res
                await ws.send_json({"type": "mmae", **res})

            def kick():
                nonlocal diag
                if diag is None or diag.done():
                    diag = asyncio.create_task(diagnose())
                    return True
                return False

            while not sim.done and not stop:
                t0 = asyncio.get_event_loop().time()
                frame = sim.step()
                await ws.send_json({"type": "frame", **frame})

                # The alarm is the trigger: the moment the CUSUM latches, ask
                # *which* fault, without waiting to be told to.
                if sim.alarm_t is not None and sim.mmae is None and diag is None:
                    kick()
                    # ...and again when each scoring's evidence window closes.
                    redo_at = sorted(
                        t for t in (sim.alarm_t + PROBE_SETTLE_MIN / 60.0,
                                    sim.alarm_t + ETHANOL_SETTLE_MIN / 60.0)
                        if t <= sim.t_final)

                # A matured window: re-score so the panel can name the cause.
                if redo_at and sim.t >= redo_at[0]:
                    if kick():
                        redo_at.pop(0)

                # drain any pending control messages without blocking
                while not inbox.empty():
                    m = inbox.get_nowait()
                    a = m.get("action")
                    if a == "inject_fault":
                        ok, note = sim.inject_fault(m.get("fault", {}))
                        await ws.send_json({"type": "fault_ack", "ok": ok, "note": note})
                    elif a == "diagnose":
                        if not kick():
                            await ws.send_json({"type": "mmae_status",
                                                "state": "busy"})
                    elif a == "stop":
                        stop = True
                    elif a == "__disconnect__":
                        stop = True

                elapsed = asyncio.get_event_loop().time() - t0
                await asyncio.sleep(max(0.0, frame_wall - elapsed))

            if diag is not None and not diag.done():
                await diag                      # let an in-flight bank finish
            # A run that alarmed late may never have had a diagnosis kicked off,
            # and one kicked off mid-run scored less evidence than now exists.
            if sim.alarm_t is not None:
                await diagnose()

            await ws.send_json({
                "type": "complete",
                "stopped": stop,
                "csv": build_csv(sim),
                "png": build_png(sim),
            })
    except WebSocketDisconnect:
        pass
    finally:
        recv_task.cancel()


# ── Static frontend (mounted last so /api and /ws take priority) ──────────
# NO-CACHE on every asset. This is a dev tool that is edited while it runs, and
# a browser holding a stale app.js looks exactly like a broken backend: new
# fields silently missing, new frame keys ignored. Correctness of what you see
# beats a few kB of caching here.
_NO_CACHE = {"Cache-Control": "no-store, must-revalidate", "Pragma": "no-cache"}


class NoCacheStatic(StaticFiles):
    def file_response(self, *args, **kwargs):
        resp = super().file_response(*args, **kwargs)
        resp.headers.update(_NO_CACHE)
        return resp


@app.get("/")
def index():
    return FileResponse(FRONTEND / "index.html", headers=_NO_CACHE)


app.mount("/", NoCacheStatic(directory=FRONTEND), name="static")
