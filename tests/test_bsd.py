"""Automated pytest suite for Blind Signal Discovery (BSD) v0.1."""

from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path

import numpy as np
import pytest

from bsd.calibration import CalibrationResult, calibrate_detector, clopper_pearson_ci, wilson_ci
from bsd.detector import DetectorConfig, TemplateFreeDetector
from bsd.evaluation import evaluate_split, summarize_records
from bsd.generator import (
    ALL_SIGNAL_FAMILIES,
    GroundTruth,
    Observation,
    generate_noise_dataset,
    generate_sample,
    generate_signal_dataset,
    make_split_rng,
)
from bsd.noise import ALL_NOISE_MODELS, generate_noise


def test_synthetic_data_reproducibility() -> None:
    """Identical seeds must produce bit-identical datasets; distinct splits must differ."""
    ds1 = generate_signal_dataset(
        split_name="test_signal",
        n_samples_count=8,
        sample_rate=1024.0,
        duration=2.0,
        signal_families=("sinusoid", "linear_chirp"),
        snr_levels=(10.0, 15.0),
        noise_models=("gaussian_white", "pink_colored"),
        master_seed=42,
    )
    ds2 = generate_signal_dataset(
        split_name="test_signal",
        n_samples_count=8,
        sample_rate=1024.0,
        duration=2.0,
        signal_families=("sinusoid", "linear_chirp"),
        snr_levels=(10.0, 15.0),
        noise_models=("gaussian_white", "pink_colored"),
        master_seed=42,
    )
    ds_other = generate_signal_dataset(
        split_name="test_withheld",
        n_samples_count=8,
        sample_rate=1024.0,
        duration=2.0,
        signal_families=("sinusoid", "linear_chirp"),
        snr_levels=(10.0, 15.0),
        noise_models=("gaussian_white", "pink_colored"),
        master_seed=42,
    )

    for o1, o2, oo in zip(ds1.observations, ds2.observations, ds_other.observations):
        assert np.array_equal(o1.signal, o2.signal)
        assert not np.allclose(o1.signal, oo.signal)


@pytest.mark.parametrize("noise_model", ALL_NOISE_MODELS)
def test_noise_only_generation(noise_model: str) -> None:
    """All noise models produce finite, correctly sized noise-only samples."""
    rng = np.random.default_rng(123)
    noise, meta = generate_noise(noise_model, n_samples=2048, sample_rate=1024.0, rng=rng, sigma=1.0)
    assert noise.shape == (2048,)
    assert np.all(np.isfinite(noise))
    assert meta["model"] == noise_model
    if noise_model != "impulsive":
        assert 0.85 <= float(np.std(noise)) <= 1.15
    else:
        # Impulsive noise has elevated kurtosis / crest factor due to spikes
        assert float(np.max(np.abs(noise))) >= 5.0


@pytest.mark.parametrize("family", ALL_SIGNAL_FAMILIES)
def test_signal_injection_correctness(family: str) -> None:
    """All signal families (including withheld) inject valid energy SNR and TF bounds."""
    rng = np.random.default_rng(77)
    obs, gt = generate_sample(
        sample_id=f"sig_{family}",
        sample_rate=1024.0,
        duration=4.0,
        noise_model="gaussian_white",
        rng=rng,
        signal_family=family,
        target_snr=15.0,
        noise_sigma=1.0,
    )
    assert obs.signal.shape == (4096,)
    assert np.all(np.isfinite(obs.signal))
    assert gt.is_signal is True
    assert gt.signal_family == family
    assert gt.target_snr == pytest.approx(15.0)
    # Empirical energy SNR should closely match target SNR (within finite-sample noise RMS variation)
    assert gt.empirical_energy_snr is not None
    assert 13.0 <= gt.empirical_energy_snr <= 17.5
    assert gt.t_start is not None and gt.t_end is not None and gt.t_start < gt.t_end
    assert gt.f_min is not None and gt.f_max is not None and gt.f_min < gt.f_max


def test_detector_input_isolation() -> None:
    """Observation must contain no ground-truth fields and detector must reject non-Observation inputs."""
    obs_field_names = {f.name for f in fields(Observation)}
    forbidden = {"is_signal", "signal_family", "noise_model", "target_snr", "t_start", "f_min", "params"}
    assert obs_field_names.isdisjoint(forbidden)

    detector = TemplateFreeDetector()
    rng = make_split_rng(99, "development")
    obs, gt = generate_sample("iso_0", 1024.0, 2.0, "gaussian_white", rng, "sinusoid", 15.0)

    with pytest.raises(TypeError):
        detector.analyze(gt)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        detector.analyze((obs, gt))  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        detector.analyze(obs.signal)  # type: ignore[arg-type]


def test_threshold_calibration_and_p_values() -> None:
    """Calibration produces a valid threshold, empirical p-value bounds, and CIs."""
    detector = TemplateFreeDetector(DetectorConfig())
    cal_split = generate_noise_dataset(
        split_name="calibration",
        n_samples_count=60,
        sample_rate=1024.0,
        duration=2.0,
        noise_models=("gaussian_white", "pink_colored"),
        master_seed=55,
    )
    cal = calibrate_detector(
        detector=detector,
        calibration_observations=cal_split.observations,
        target_fap=0.05,
        calibration_noise_models=("gaussian_white", "pink_colored"),
        calibration_seed=55,
    )
    assert isinstance(cal, CalibrationResult)
    assert cal.n_calibration_samples == 60
    assert cal.threshold > 0.0
    assert cal.empirical_calibration_fap <= 0.05 + 1e-9
    assert cal.bonferroni_per_bin_exp_threshold > cal.naive_per_bin_exp_threshold

    # Score above maximum calibration score yields exact upper-bound empirical p-value
    p_hi, is_ub = cal.empirical_p_value(cal.sorted_null_scores[-1] + 100.0)
    assert is_ub is True
    assert p_hi == pytest.approx(1.0 / 61.0)

    lo_w, hi_w = wilson_ci(5, 500)
    lo_cp, hi_cp = clopper_pearson_ci(5, 500)
    assert 0.0 < lo_w < 0.01 < hi_w < 0.05
    assert 0.0 < lo_cp < 0.01 < hi_cp < 0.05


@pytest.mark.parametrize("family", ALL_SIGNAL_FAMILIES)
def test_candidate_extraction_strong_signals(family: str) -> None:
    """Strong injected signals (SNR=25) across all families must be detected and localized."""
    detector = TemplateFreeDetector(DetectorConfig())
    cal_split = generate_noise_dataset(
        split_name="calibration",
        n_samples_count=50,
        sample_rate=1024.0,
        duration=4.0,
        noise_models=("gaussian_white", "pink_colored", "nonstationary"),
        master_seed=11,
    )
    cal = calibrate_detector(detector, cal_split.observations, target_fap=0.02, calibration_seed=11)

    rng = make_split_rng(11, f"test_strong_{family}")
    obs, gt = generate_sample(
        sample_id=f"strong_{family}",
        sample_rate=1024.0,
        duration=4.0,
        noise_model="pink_colored",
        rng=rng,
        signal_family=family,
        target_snr=25.0,
    )
    res = detector.analyze(obs, calibration=cal, return_maps=True)
    assert res.detected is True
    assert len(res.candidates) >= 1
    assert res.tf_maps is not None
    top = res.candidates[0]
    assert top.anomaly_score > cal.threshold
    assert top.empirical_p_value is not None and top.empirical_p_value <= 0.05
    # Verify candidate overlaps true frequency band
    assert gt.f_min is not None and gt.f_max is not None
    assert top.f_max >= gt.f_min and top.f_min <= gt.f_max


def test_result_serialization_and_end_to_end(tmp_path: Path) -> None:
    """End-to-end execution produces valid JSON, CSV, figures, and Markdown report."""
    import sys

    exp_dir = Path(__file__).resolve().parent.parent / "experiments"
    if str(exp_dir) not in sys.path:
        sys.path.insert(0, str(exp_dir))
    import run_experiment as runner

    cfg_text = f"""
[experiment]
name = "pytest_e2e"
seeds = [7]
output_dir = "{str(tmp_path / 'out').replace(chr(92), '/')}"

[signal]
sample_rate = 1024.0
duration = 2.0
snr_levels = [8.0, 20.0]
known_families = ["sinusoid", "broadband_pulse", "linear_chirp", "intermittent_narrowband"]
withheld_families = ["sinusoidal_fm"]

[noise]
calibration_models = ["gaussian_white", "pink_colored", "nonstationary"]
evaluation_models = ["gaussian_white", "pink_colored", "nonstationary"]
shifted_models = ["shifted_ar2"]
impulsive_models = ["impulsive"]

[detector]
nperseg = 128
noverlap = 96
window = "hann"
freq_median_kernel = 11
time_median_kernel = 15
bin_seed_threshold = 4.0
min_cluster_bins = 2

[calibration]
target_fap = 0.02
n_calibration_samples = 40

[evaluation]
n_noise_test = 30
n_signal_test = 32
n_withheld_test = 12
n_shifted_noise = 15
n_shifted_signal = 16
n_robustness_per_condition = 10
"""
    cfg_path = tmp_path / "test_cfg.toml"
    cfg_path.write_text(cfg_text, encoding="utf-8")

    res = runner.run_experiment(cfg_path)
    out_dir = tmp_path / "out"
    assert (out_dir / "results.json").is_file()
    assert (out_dir / "per_sample_results.csv").is_file()
    assert (out_dir / "summary_by_family_snr.csv").is_file()
    assert (out_dir / "report.md").is_file()

    loaded = json.loads((out_dir / "results.json").read_text(encoding="utf-8"))
    assert loaded["metadata"]["experiment_name"] == "pytest_e2e"
    assert len(loaded["figures"]) == 9
    for fig_name in loaded["figures"]:
        assert (out_dir / "figures" / fig_name).is_file()
    assert res["held_out_noise_only"]["n_noise"] == 30
