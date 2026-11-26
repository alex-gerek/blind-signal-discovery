"""Availability accounting and additive precision contracts, without O3 input."""
from __future__ import annotations

import csv
import importlib.util
import io
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from bsd import workstation_runner as wr
from bsd.adaptive_regime_v2b import CalibrationUncertaintySummary, DualOutputWindowV2b
from bsd.detector_v2 import MultiResolutionCandidate
from bsd.io_loader import InputValidationError


def window(index=0, available=True, score=2.88893, threshold=2.888922, nominated=None):
    if not available:
        threshold = None
    candidate = MultiResolutionCandidate(1, "long_500ms", .123456789, 1.23456789, 100.123456789, 201.23456789,
                                         .7654321, 150.123456, 123.456789, score, 9., 4, 1/155, True)
    flag = (score > threshold if threshold is not None else True) if nominated is None else nominated
    u = CalibrationUncertaintySummary(154, 140., 0., .9, 1/155, 1/155, (.001, .03), (.001, .03), False)
    return DualOutputWindowV2b(index, f"w{index}", (candidate,), (candidate,), score, threshold, flag, flag,
        5.123456789, 4.123456789 if available else None, False, False, None, 0, "CALIBRATED_STATIONARY",
        "CALIBRATED_STATIONARY", "NO_GLITCH_STORM", (), "EMPIRICAL_PROVISIONAL", u, "causal_trailing", "transient", False,
        0, 0., 0., diagnostics={"union_stft_threshold": threshold}, detection_available=available)


def result(windows, mode="transient"):
    candidates, timeline = wr._build_rows(windows, 4., mode)
    serializable = SimpleNamespace(to_dict=lambda: {})
    return wr.WorkstationRunResult({}, {}, serializable, serializable, {}, windows, [], candidates, timeline)


@pytest.mark.parametrize("mode", ["transient", "persistent"])
def test_all_abstained_has_no_available_rate_or_efficiency(mode):
    negatives = [window(0, False), window(1, False)]
    positives = [{"detected": False, "detection_available": False, "evaluation_status": "ABSTAINED"}]
    audit = wr._self_test_availability_metrics(negatives, positives, 4., mode)
    assert audit["negative_control_available_exposure_s"] == 0
    assert audit["negative_control_abstained_windows"] == 2
    assert audit["negative_control_far_available"] is None
    assert audit["negative_control_far_available_ci90"] is None
    assert audit["negative_control_far_status"] == "UNAVAILABLE_NO_EXPOSURE"
    assert audit["positive_injections_abstained"] == 1
    assert audit["positive_injections_missed_available"] == 0
    assert audit["detection_efficiency_available"] is None
    assert audit["detection_efficiency_available_ci90"] is None
    json.dumps(audit, allow_nan=False)


def test_mixed_availability_is_separate_from_legacy_counts(monkeypatch):
    windows = [window(i, False) for i in range(10)]
    windows[1] = window(1, True, nominated=False)
    windows[2] = window(2, True, nominated=True)
    windows[4] = window(4, True, nominated=True)
    windows[6] = window(6, True, nominated=False)
    # Three injections occupy indices 4, 6, 8: detected, missed, abstained.
    monkeypatch.setattr(wr, "run_workstation_analysis", lambda *args: result(windows))
    output = wr.run_synthetic_self_test(n_bg_windows=60, n_neg_control_windows=3, n_pos_injections=3, noise_model="gaussian_white")
    audit = output.ground_truth_audit
    assert audit["negative_control_far"] == .3333  # legacy denominator retained
    assert audit["positive_injections_missed"] == 2  # legacy non-detection count retained
    assert audit["detection_efficiency"] == .3333
    assert audit["negative_control_available_windows"] == 2
    assert audit["negative_control_available_exposure_s"] == 8.
    assert audit["negative_control_total_exposure_s"] == 12.
    assert audit["negative_control_far_available"] == .5
    assert audit["positive_injections_available"] == 2
    assert audit["positive_injections_missed_available"] == 1
    assert audit["positive_injections_abstained"] == 1
    assert audit["detection_efficiency_available"] == .5
    assert [r["evaluation_status"] for r in audit["detected_records"]] == ["DETECTED"]
    assert {r["evaluation_status"] for r in audit["missed_records"]} == {"MISSED", "ABSTAINED"}
    assert len(audit["abstained_records"]) == len(audit["missed_available_records"]) == 1


def test_no_evaluation_windows_has_null_conditional_metrics():
    audit = wr._self_test_availability_metrics([], [], 4., "transient")
    assert audit["negative_control_far_available"] is None
    assert audit["detection_efficiency_available"] is None
    assert audit["negative_control_available_exposure_s"] == 0


@pytest.mark.parametrize("duration", [0., -1., np.nan, np.inf])
def test_exposure_must_be_finite_and_positive(duration):
    with pytest.raises(InputValidationError):
        wr._self_test_availability_metrics([], [], duration, "transient")


@pytest.mark.parametrize("mode", ["transient", "persistent"])
def test_exact_crossing_round_trips_without_changing_legacy_numbers(tmp_path, mode):
    threshold = 2.888922
    score = np.nextafter(threshold, np.inf)
    w = window(score=score, threshold=threshold)
    if mode == "persistent":
        w = replace(w, broadband_excess_score=score, broadband_threshold=threshold, broadband_nominated=True)
    output = result([w], mode)
    row = output.timeline_rows[0]
    prefix = "stft" if mode == "transient" else "broadband"
    legacy_score = "stft_max_score" if mode == "transient" else "broadband_excess_score"
    assert row[legacy_score] == row[prefix + "_threshold"] == 2.8889
    assert row["active_mode_anomaly_flag"] is True
    assert row["active_score_full_precision"] == score
    assert row["active_threshold_full_precision"] == threshold
    paths = wr.export_workstation_run(output, tmp_path / mode)
    parsed = json.loads(paths["run_report_json"].read_text())
    saved = parsed["window_availability_timeline"][0]
    assert saved["active_score_full_precision"] > saved["active_threshold_full_precision"]
    assert saved["active_score_full_precision"] == score
    assert saved["evaluation_status"] == "DETECTED"
    assert saved["precision_status"] == "FULL_PRECISION"
    csv_row = next(csv.DictReader(io.StringIO(paths["timeline_csv"].read_text())))
    assert float(csv_row["active_score_full_precision"]) == score
    assert float(csv_row["active_threshold_full_precision"]) == threshold
    assert float(csv_row["active_score_full_precision"]) > float(csv_row["active_threshold_full_precision"])
    assert parsed == wr.build_workstation_report(output)
    for key, content in wr.build_workstation_csv(output).items():
        assert paths[key].read_bytes() == content.encode("utf-8")
    candidates = list(csv.DictReader(io.StringIO(paths["candidates_csv"].read_text())))
    assert float(candidates[0]["empirical_p_value_full_precision"]) == 1/155
    assert float(candidates[0]["t_abs_start_s_full_precision"]) == .123456789
    assert parsed["candidate_records"][0]["standardized_score_full_precision"] == score


def test_unavailable_record_keeps_diagnostics_and_null_threshold(tmp_path):
    output = result([window(available=False, score=20.)])
    row = output.timeline_rows[0]
    assert row["active_mode_diagnostic_nomination"] is True
    assert row["active_mode_anomaly_flag"] is False
    assert row["evaluation_status"] == "ABSTAINED"
    assert row["active_threshold_full_precision"] is None
    assert output.candidates_rows  # diagnostic candidate preserved
    paths = wr.export_workstation_run(output, tmp_path / "abstained")
    data = json.loads(paths["run_report_json"].read_text(), parse_constant=lambda x: pytest.fail(x))
    assert data["window_availability_timeline"][0]["active_threshold_full_precision"] is None
    assert data["candidate_records"][0]["channel_threshold_full_precision"] is None
    csv_row = next(csv.DictReader(io.StringIO(paths["timeline_csv"].read_text())))
    assert csv_row["active_threshold_full_precision"] == ""


def test_legacy_precision_is_not_reconstructed_from_rounded_values():
    output = result([window()])
    output.timeline_rows[0] = {k:v for k,v in output.timeline_rows[0].items() if not k.endswith("_full_precision")}
    saved = wr.build_workstation_report(output)["window_availability_timeline"][0]
    assert saved["precision_status"] == "UNAVAILABLE_LEGACY_RECORD"
    assert saved["active_score_full_precision"] is None


def test_empty_csv_keeps_schema_and_nonfinite_report_writes_nothing(tmp_path):
    output = result([])
    texts = wr.build_workstation_csv(output)
    assert "active_score_full_precision" in next(csv.reader(io.StringIO(texts["timeline_csv"])))
    assert "channel_threshold_full_precision" in next(csv.reader(io.StringIO(texts["candidates_csv"])))
    output.summary_metrics["invalid"] = float("nan")
    with pytest.raises(ValueError):
        wr.export_workstation_run(output, tmp_path / "invalid")
    assert not (tmp_path / "invalid").exists()
    with pytest.raises(InputValidationError):
        result([window(score=float("inf"))])


def test_ui_uses_unavailable_and_shared_download_builders():
    app = Path(__file__).parents[1] / "app/bsd_workstation.py"
    spec = importlib.util.spec_from_file_location("bsd_metrics_ui_test", app)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    gt = {"negative_control_far_available": None, "negative_control_available_windows": 0,
          "negative_control_available_exposure_s": 0, "negative_control_false_alarms": 0,
          "positive_injections_available": 0, "positive_injections_detected": 0,
          "positive_injections_missed_available": 0, "positive_injections_abstained": 2,
          "detection_efficiency_available": None}
    display = module._self_test_display_values(gt)
    assert display["far"] == display["efficiency"] == "UNAVAILABLE"
    assert display["abstained"] == "2"
    assert module.build_workstation_report is wr.build_workstation_report
    assert module.build_workstation_csv is wr.build_workstation_csv


def test_actual_streamlit_download_bytes_match_saved_contract(monkeypatch):
    from contextlib import nullcontext
    output = wr.run_synthetic_self_test(wr.WorkstationConfig(target_fap=.02), n_bg_windows=80,
                                       n_neg_control_windows=2, n_pos_injections=1, seed=41)
    app = Path(__file__).parents[1] / "app/bsd_workstation.py"
    spec = importlib.util.spec_from_file_location("bsd_download_contract_test", app)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    class FakeStreamlit:
        def __init__(self):
            self.session_state = {"bsd_last_result": output}
            self.sidebar = nullcontext()
            self.downloads = {}
        def __enter__(self): return self
        def __exit__(self, *args): return False
        def radio(self, *args, **kwargs): return "Local File Analysis"
        def selectbox(self, label, options, **kwargs): return options[0]
        def number_input(self, label, *args, **kwargs): return args[2]
        def checkbox(self, *args, **kwargs): return False
        def text_input(self, *args, **kwargs): return "arbitrary"
        def button(self, *args, **kwargs): return False
        def file_uploader(self, *args, **kwargs): return None
        def expander(self, *args, **kwargs): return nullcontext()
        def tabs(self, labels): return [nullcontext() for _ in labels]
        def columns(self, n): return [self for _ in range(n)]
        def slider(self, label, *args, **kwargs): return args[2]
        def download_button(self, label, data, filename, mime): self.downloads[filename] = data
        def __getattr__(self, name): return lambda *args, **kwargs: None

    fake = FakeStreamlit()
    monkeypatch.setattr(module, "st", fake)
    monkeypatch.setattr(module, "_plot_timeline", lambda result: None)
    monkeypatch.setattr(module, "_plot_window_waveform_and_spec", lambda result, index: None)
    module.main()
    assert fake.downloads["timeline.csv"] == wr.build_workstation_csv(output)["timeline_csv"].encode("utf-8")
    assert fake.downloads["candidates.csv"] == wr.build_workstation_csv(output)["candidates_csv"].encode("utf-8")
    assert fake.downloads["run_report.json"] == json.dumps(wr.build_workstation_report(output), indent=2, allow_nan=False).encode("utf-8")
