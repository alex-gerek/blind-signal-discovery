"""BSD — Local Measurement Workstation (Minimal Streamlit GUI).

Thin presentation layer over `bsd.io_loader` and `bsd.workstation_runner`.
Zero detector mathematics is reimplemented here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal as sp_signal
import streamlit as st

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bsd.io_loader import (
    SUPPORTED_EXTENSIONS,
    InputValidationError,
    O3HoldoutBlockedError,
    load_timeseries_file,
)
from bsd.upload_storage import save_uploaded_bytes
from bsd.workstation_runner import (
    CAL_STATUS_INSUFFICIENT,
    CAL_STATUS_PROVISIONAL,
    CAL_STATUS_VALIDATED,
    WorkstationConfig,
    WorkstationRunResult,
    export_workstation_run,
    build_workstation_report,
    build_workstation_csv,
    run_synthetic_self_test,
    run_workstation_analysis,
)

st.set_page_config(page_title="BSD Local Measurement Workstation", layout="wide")


def _save_upload_to_local_tmp(uploaded_file: Any) -> Path:
    return save_uploaded_bytes(
        uploaded_file.name, uploaded_file.getvalue(), PROJECT_ROOT / ".tmp" / "uploads"
    )


def _plot_window_waveform_and_spec(result: WorkstationRunResult, win_idx: int) -> plt.Figure:
    obs = result.target_observations[win_idx]
    win = result.windows[win_idx]
    t = np.arange(obs.signal.size, dtype=np.float64) / obs.sample_rate
    units = result.run_metadata.get("target_preprocessing", {}).get("units", "arbitrary")
    win_dur = float(result.run_metadata.get("config", {}).get("window_duration_s", obs.duration))
    t_abs_0 = win_idx * win_dur
    t_abs_1 = (win_idx + 1) * win_dur

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.2, 5.4), sharex=True)
    ax1.plot(t, obs.signal, color="#1f4e79", lw=0.8)
    ax1.set_ylabel(f"Amplitude ({units})")
    ax1.set_title(
        f"Window {win_idx} ({obs.sample_id}) | Absolute Time: [{t_abs_0:.1f}, {t_abs_1:.1f}] s | "
        f"Regime: {win.credibility_status} | "
        f"Warnings: {', '.join(win.warning_causes) if win.warning_causes else 'None'}"
    )
    ax1.grid(alpha=0.25)

    freqs, times, zxx = sp_signal.stft(
        obs.signal - np.median(obs.signal),
        fs=obs.sample_rate,
        window="hann",
        nperseg=128,
        noverlap=96,
        boundary=None,
        padded=False,
    )
    pwr_db = 10.0 * np.log10(np.maximum(np.abs(zxx) ** 2, 1e-24))
    pcm = ax2.pcolormesh(times, freqs, pwr_db, shading="auto", cmap="viridis")
    fig.colorbar(pcm, ax=ax2, label="Relative STFT Power (10 log10 |Z|^2, dB re 1 unit^2; Hann nperseg=128)")
    for c in win.raw_stft_candidates:
        rect = plt.Rectangle(
            (c.t_start, c.f_min),
            max(0.03, c.t_end - c.t_start),
            max(4.0, c.f_max - c.f_min),
            fill=False,
            edgecolor="#ff4b4b",
            linewidth=1.5,
        )
        ax2.add_patch(rect)
    ax2.set_ylabel("Frequency (Hz)")
    ax2.set_xlabel("Time Within Window (s)")
    ax2.set_ylim(0, obs.sample_rate / 2.0)
    fig.tight_layout()
    return fig


def _plot_timeline(result: WorkstationRunResult) -> plt.Figure:
    xs = [w.index for w in result.windows]
    mode = result.run_metadata["config"]["search_mode"]
    if mode == "transient":
        scores = [w.stft_max_score for w in result.windows]
        thresh = [w.stft_threshold for w in result.windows]
        ylabel = "Standardized Global Anomaly Score (dimensionless)"
    else:
        scores = [w.broadband_excess_score for w in result.windows]
        thresh = [w.broadband_threshold for w in result.windows]
        ylabel = "Frozen-Ref Broadband Excess Score (dimensionless)"

    thresh_plot = [float(t) if (t is not None and np.isfinite(t)) else np.nan for t in thresh]
    finite_th = [t for t in thresh_plot if np.isfinite(t)]
    n_unavail = len(thresh_plot) - len(finite_th)
    if finite_th:
        th_min, th_max = min(finite_th), max(finite_th)
        base_lbl = (
            f"Active Calibrated Threshold = {th_min:.4f} (dimensionless)"
            if abs(th_max - th_min) < 1e-6
            else f"Active Calibrated Threshold = {th_min:.4f}..{th_max:.4f} (dimensionless)"
        )
        th_lbl = f"{base_lbl}; {n_unavail} window(s) UNAVAILABLE" if n_unavail > 0 else base_lbl
    else:
        th_lbl = "Active Calibrated Threshold = UNAVAILABLE (null)"

    fig, ax = plt.subplots(figsize=(10.2, 3.2))
    ax.plot(xs, scores, color="#1f4e79", lw=1.2, label=f"Stream A Score ({mode}, dimensionless)")
    ax.plot(xs, thresh_plot, color="#d62728", ls="--", lw=1.3, label=th_lbl)
    for w in result.windows:
        if "CONTINUUM_DRIFT" in w.warning_causes:
            ax.axvspan(w.index - 0.5, w.index + 0.5, color="#fbd38d", alpha=0.25)
        if "GLITCH_STORM" in w.warning_causes:
            ax.axvspan(w.index - 0.5, w.index + 0.5, color="#d6bcfa", alpha=0.25)
        if not w.detection_available or "CALIBRATION_UNAVAILABLE" in w.warning_causes:
            ax.axvspan(w.index - 0.5, w.index + 0.5, color="#cbd5e0", alpha=0.35)
    ax.set_xlabel("Window Index")
    ax.set_ylabel(ylabel)
    ax.set_title(
        "Full-Timeline Anomaly Score vs. Threshold "
        "(Orange = Continuum Drift, Purple = Glitch Storm, Gray = Cal Unavailable)"
    )
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig


def _self_test_display_values(gt: dict[str, Any]) -> dict[str, str]:
    """Display current availability-conditioned fields; never invent legacy exposure."""
    far = gt.get("negative_control_far_available")
    n_neg = gt.get("negative_control_available_windows")
    n_pos = gt.get("positive_injections_available")
    efficiency = gt.get("detection_efficiency_available")
    ci = gt.get("detection_efficiency_available_ci90")
    return {
        "far": "UNAVAILABLE" if far is None else f"{far * 100:.2f}% ({gt['negative_control_false_alarms']}/{n_neg})",
        "efficiency": "UNAVAILABLE" if efficiency is None or ci is None else f"{efficiency * 100:.1f}% [{ci[0] * 100:.0f}%, {ci[1] * 100:.0f}%]",
        "exposure": "UNAVAILABLE" if n_neg is None else f"{n_neg} windows / {gt['negative_control_available_exposure_s']:g} s",
        "detected": f"{gt['positive_injections_detected']} / {n_pos}" if n_pos is not None else "UNAVAILABLE (legacy audit)",
        "missed": f"{gt['positive_injections_missed_available']} / {n_pos}" if n_pos is not None else "UNAVAILABLE (legacy audit)",
        "abstained": str(gt['positive_injections_abstained']) if n_pos is not None else "UNAVAILABLE (legacy audit)",
    }


def main() -> None:
    st.title("Blind Signal Discovery (BSD) — Local Measurement Workstation")
    st.info(
        "**Research Prototype Notice (v0.3 Phase 2b):** BSD distinguishes *measured features*, "
        "*statistical anomalies*, and *detection candidates*, and **never** classifies a single-channel "
        "anomaly as an independently confirmed physical event. Reserved O3 holdout data is hard-blocked."
    )

    with st.sidebar:
        st.header("1. Measurement Source")
        run_kind = st.radio(
            "Workflow",
            ("Synthetic Injection Self-Test", "Local File Analysis"),
            index=0,
        )

        st.header("2. Search & Statistical Settings")
        search_mode = st.selectbox("Search Mode", ("transient", "persistent"), index=0)
        proc_mode = st.selectbox("Processing Mode", ("causal_trailing", "offline_symmetric"), index=0)
        if proc_mode == "offline_symmetric":
            st.warning("`offline_symmetric` uses future windows (`future_data_used=True`).")
        if search_mode == "persistent":
            st.warning(
                "Single-channel `persistent` broadband excess is physically degenerate with background drift "
                "(`SINGLE_CHANNEL_BROADBAND_AMBIGUOUS`)."
            )
        target_fap = st.number_input("Target False-Alarm Probability (FAP)", 0.002, 0.10, 0.01, step=0.005, format="%.3f")
        min_nth_req = int(np.ceil((1.0 - float(target_fap)) / float(target_fap)))
        min_bg_req = max(60, int(np.ceil(min_nth_req / 0.70)))
        st.caption(f"Minimum background windows for attainable `target_fap={target_fap:.3f}`: `{min_bg_req}` (`N_thresh >= {min_nth_req}`).")
        cal_seed = int(st.number_input("Calibration Seed", 0, 999999, 0, step=1))

        with st.expander("Optional Heuristic Vetoes (Off by Default for Arbitrary Sensors)"):
            en_glitch = st.checkbox("Enable Glitch Veto (reversible)", value=False)
            en_line = st.checkbox("Enable Narrowband Line Veto (reversible)", value=False)

        cfg = WorkstationConfig(
            search_mode=search_mode,
            processing_mode=proc_mode,
            window_duration_s=4.0,
            target_fap=float(target_fap),
            enable_glitch_veto=en_glitch,
            enable_line_veto=en_line,
            calibration_seed=cal_seed,
        )

        run_btn = False
        if run_kind == "Synthetic Injection Self-Test":
            st.header("3. Self-Test Parameters")
            n_bg = int(st.slider("Calibration Background Windows (4 s each)", 60, 260, max(160, min_bg_req), step=20))
            n_neg = int(st.slider("Negative Control Windows (Pure Noise)", 20, 80, 40, step=10))
            n_pos = int(st.slider("Positive Injections", 4, 20, 10, step=2))
            noise_model = st.selectbox("Noise Model", ("pink_colored", "gaussian_white", "shifted_ar2"))
            sig_family = st.selectbox(
                "Injected Signal Family", ("linear_chirp", "sinusoid", "broadband_pulse", "sinusoidal_fm")
            )
            snr = float(
                st.slider(
                    "Injected Root-Energy SNR (||s||_2 / sigma, linear dimensionless — NOT dB)",
                    8.0,
                    30.0,
                    18.0,
                    step=1.0,
                )
            )
            st_seed = int(st.number_input("Self-Test Seed", 0, 999999, 2026, step=1))
            run_btn = st.button("Run Synthetic Self-Test", type="primary")
        else:
            st.header("3. Local Files (Separate Target & Background)")
            st.caption(f"Supported formats: {', '.join(SUPPORTED_EXTENSIONS)} (sample_rate >= 1024 Hz)")
            tgt_up = st.file_uploader("Target Dataset", type=[e.lstrip(".") for e in SUPPORTED_EXTENSIONS])
            bg_up = st.file_uploader("Independent Background / Calibration Dataset (>= 60 windows = 240 s)", type=[e.lstrip(".") for e in SUPPORTED_EXTENSIONS])
            allow_resamp = st.checkbox("Allow explicit anti-aliased resampling to 1024 Hz (Nyquist 512 Hz)", value=False)
            man_fs = st.number_input("Manual Sample Rate Hz (for .npy / headerless .csv, 0 = auto)", 0.0, 1000000.0, 0.0, step=1.0)
            ch_idx_raw = st.number_input("Channel Index (for multi-channel .wav/.npy, -1 = auto/mono)", -1, 64, -1, step=1)
            csv_value_col = int(st.number_input("CSV value column (0-based)", 0, 1000, 0, step=1))
            csv_time_col_raw = int(st.number_input("CSV time column in seconds (-1 = none)", -1, 1000, -1, step=1))
            st.caption("CSV columns apply to both files. Time columns must be finite, increasing and uniformly sampled; supply a common time reference.")
            units_lbl = st.text_input("Units Label", value="arbitrary")
            run_btn = st.button("Run Local File Measurement", type="primary")

    if run_btn:
        # A rejected new input must not leave the previous run available as its result.
        st.session_state.pop("bsd_last_result", None)
        try:
            with st.spinner("Running BSD calibration and detection..."):
                if run_kind == "Synthetic Injection Self-Test":
                    result = run_synthetic_self_test(
                        cfg,
                        n_bg_windows=n_bg,
                        n_neg_control_windows=n_neg,
                        n_pos_injections=n_pos,
                        noise_model=noise_model,
                        signal_family=sig_family,
                        target_snr=snr,
                        seed=st_seed,
                    )
                else:
                    if tgt_up is None or bg_up is None:
                        raise InputValidationError("Both a Target file and a separate Background file are required.")
                    tgt_path = _save_upload_to_local_tmp(tgt_up)
                    bg_path = _save_upload_to_local_tmp(bg_up)
                    fs_opt = float(man_fs) if man_fs > 0 else None
                    ch_opt = int(ch_idx_raw) if ch_idx_raw >= 0 else None
                    tgt_ts = load_timeseries_file(
                        tgt_path, sample_rate_hz=fs_opt, channel_index=ch_opt, allow_resample=allow_resamp, units=units_lbl,
                        csv_value_column=csv_value_col, csv_time_column=csv_time_col_raw if csv_time_col_raw >= 0 else None
                    )
                    bg_ts = load_timeseries_file(
                        bg_path, sample_rate_hz=fs_opt, channel_index=ch_opt, allow_resample=allow_resamp, units=units_lbl,
                        csv_value_column=csv_value_col, csv_time_column=csv_time_col_raw if csv_time_col_raw >= 0 else None
                    )
                    result = run_workstation_analysis(tgt_ts, bg_ts, cfg)
            st.session_state["bsd_last_result"] = result
        except (InputValidationError, O3HoldoutBlockedError, FileNotFoundError) as exc:
            st.error(f"Validation / Safety Error: {exc}")
        except Exception as exc:
            st.error(f"Unexpected Error: {exc}")

    result: WorkstationRunResult | None = st.session_state.get("bsd_last_result")
    if result is None:
        st.caption("Configure a Synthetic Self-Test or upload Target + Independent Background files in the sidebar to start.")
        return

    tab1, tab2, tab3, tab4 = st.tabs(
        [
            "1. Calibration Status & Warnings",
            "2. Waveform & Spectrogram",
            "3. Candidates & Self-Test Audit",
            "4. Export & Run Metadata",
        ]
    )

    with tab1:
        ca = result.calibration_audit
        if ca.status == CAL_STATUS_VALIDATED:
            st.success(
                f"Calibration Status: **{ca.status}** "
                "(Internal finite-sample & linear lag-1..10 screening met; does NOT guarantee real-sensor stationarity or real-data FAR control)"
            )
        elif ca.status == CAL_STATUS_PROVISIONAL:
            st.warning(f"Calibration Status: **{ca.status}** — Unmet validation criteria: {'; '.join(ca.unmet_validation_conditions)}")
        else:
            st.error(f"Calibration Status: **{ca.status}** — {'; '.join(ca.unmet_validation_conditions)}")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Background Windows", f"{ca.n_bg_windows}")
        c2.metric("Effective Samples (N_eff / N_th)", f"{ca.n_eff:.1f} / {ca.n_thresh}")
        c3.metric("Ljung-Box Linear-Lag p", f"{ca.ljung_box_p:.3f}")
        c4.metric("90% FAR CI (N_eff)", f"[{ca.far_ci90_neff[0]*100:.2f}%, {ca.far_ci90_neff[1]*100:.2f}%]")

        sm = result.summary_metrics
        if sm.get("n_calibration_unavailable_windows", 0) > 0:
            st.warning(
                f"Detection was **UNAVAILABLE** on `{sm['n_calibration_unavailable_windows']} / {sm['n_target_windows']}` "
                "window(s) because calibration did not yet have enough clean background windows to attain the requested "
                "`target_fap`. Unavailable windows are excluded from detection-available and clean-calibrated exposure."
            )
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Clean / Det.-Avail. Duty Cycle", f"{sm['clean_duty_cycle']*100:.1f}% / {sm.get('detection_available_duty_cycle', 1.0)*100:.1f}%")
        m2.metric("Active-Mode Anomaly Windows", f"{sm['n_active_mode_anomaly_windows']} / {sm.get('n_detection_available_windows', sm['n_target_windows'])}")
        m3.metric("Continuum Drift Warnings", f"{sm['n_continuum_drift_warnings']}")
        m4.metric("Glitch-Storm Warnings", f"{sm['n_glitch_storm_warnings']}")

        st.subheader("Background Independence & Preprocessing Provenance")
        if result.independence_audit.warnings:
            for w_msg in result.independence_audit.warnings:
                st.warning(w_msg)
        else:
            st.success(f"Independence Check: `{result.independence_audit.status}`")

        col_t, col_b = st.columns(2)
        with col_t:
            st.markdown("**Target Preprocessing**")
            st.json(result.run_metadata["target_preprocessing"])
        with col_b:
            st.markdown("**Background Preprocessing**")
            st.json(result.run_metadata["background_preprocessing"])

        st.subheader("Scientific Output Semantics")
        st.json(result.scientific_semantics)

    with tab2:
        st.pyplot(_plot_timeline(result))
        win_idx = st.slider("Inspect Target Window Index", 0, len(result.windows) - 1, 0)
        st.pyplot(_plot_window_waveform_and_spec(result, win_idx))

    with tab3:
        st.subheader("Preserved Raw Candidates & Statistical Anomalies (Stream A)")
        th_now = result.windows[0].stft_threshold if result.windows else None
        th_str = f"`{th_now:.4f}`" if th_now is not None else "`UNAVAILABLE (null)`"
        st.caption(
            f"**Decision Rule:** Candidates are nominated when **`standardized_score > active_threshold`** "
            f"(dimensionless; initial threshold = {th_str}). `raw_anomaly_score` is the unstandardized "
            f"cluster excess-power sum and is **not** compared directly against the threshold. "
            f"Minimum resolvable empirical p-value is `p_min = 1 / (N_thresh + 1) = {result.calibration_audit.p_min:.4f}`."
        )
        if result.candidates_rows:
            st.dataframe(pd.DataFrame(result.candidates_rows), use_container_width=True)
        else:
            st.info("No candidates or statistical anomalies exceeded threshold in this run.")

        if result.ground_truth_audit is not None:
            gt = result.ground_truth_audit
            st.subheader("Synthetic Self-Test Ground-Truth Audit")
            display = _self_test_display_values(gt)
            g1, g2, g3 = st.columns(3)
            g1.metric("Negative Control FAR (available)", display["far"])
            g2.metric("Conditional Efficiency (90% CI)", display["efficiency"])
            g3.metric("Available Noise Exposure", display["exposure"])
            g4, g5, g6 = st.columns(3)
            g4.metric("Positive Detected / Available", display["detected"])
            g5.metric("Positive Missed / Available", display["missed"])
            g6.metric("Positive Abstained", display["abstained"])
            st.caption("Conditional rates exclude unavailable decisions. Legacy all-window counters remain in the export with explicit semantics; Wilson intervals assume independent evaluation windows.")

            st.markdown("**Detected Injections**")
            if gt["detected_records"]:
                st.dataframe(pd.DataFrame(gt["detected_records"]), use_container_width=True)
            else:
                st.caption("None detected.")
            st.markdown("**Missed Injections While Detection Was Available**")
            if gt.get("missed_available_records"):
                st.dataframe(pd.DataFrame(gt["missed_available_records"]), use_container_width=True)
            else:
                st.caption("No recorded misses among available decisions. Check abstentions and available-injection count before interpreting sensitivity.")
            st.markdown("**Abstained Injections (No Available Decision)**")
            if gt.get("abstained_records"):
                st.dataframe(pd.DataFrame(gt["abstained_records"]), use_container_width=True)
            else:
                st.caption("No abstentions recorded in the current availability audit; legacy audits may lack this field.")

    with tab4:
        st.subheader("Reproducible Export (CSV + JSON)")
        if st.button("Save Reproducible Bundle Locally"):
            out_dir = PROJECT_ROOT / "results" / "workstation_runs" / result.run_metadata["run_id"]
            paths = export_workstation_run(result, out_dir)
            st.success(f"Exported to `{out_dir.relative_to(PROJECT_ROOT)}` ({len(paths)} files).")

        csv_text = build_workstation_csv(result)
        d1, d2, d3 = st.columns(3)
        d1.download_button("Download candidates.csv", csv_text["candidates_csv"].encode("utf-8"), "candidates.csv", "text/csv")
        d2.download_button("Download timeline.csv", csv_text["timeline_csv"].encode("utf-8"), "timeline.csv", "text/csv")
        report_json = json.dumps(build_workstation_report(result), indent=2, allow_nan=False)
        d3.download_button("Download run_report.json", report_json.encode("utf-8"), "run_report.json", "application/json")

        st.subheader("Complete Run Metadata")
        st.json(result.run_metadata)


if __name__ == "__main__":
    main()
