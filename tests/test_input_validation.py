"""Input-integrity regressions; all data are synthetic and O3 metadata is stubbed."""
from __future__ import annotations

import importlib.util
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from bsd import io_loader
from bsd.io_loader import InputValidationError, check_background_independence, load_timeseries_file
from bsd.upload_storage import save_uploaded_bytes


@pytest.fixture(autouse=True)
def synthetic_holdout_metadata(monkeypatch):
    # Exercise filename guard mechanics without reading the real reserved manifest.
    monkeypatch.setattr(io_loader, "_load_o3_event_ids", lambda: ("TEST_RESERVED_EVENT",))


def test_uploads_keep_distinct_contents_and_basename(tmp_path):
    paths = [save_uploaded_bytes("recording.csv", payload, tmp_path) for payload in (b"target", b"background")]
    assert paths[0] != paths[1]
    assert [p.name for p in paths] == ["recording.csv", "recording.csv"]
    assert [p.read_bytes() for p in paths] == [b"target", b"background"]
    assert all(p.resolve().is_relative_to(tmp_path.resolve()) for p in paths)
    with ThreadPoolExecutor(max_workers=4) as pool:
        concurrent = list(pool.map(lambda i: save_uploaded_bytes("recording.csv", str(i).encode(), tmp_path), range(12)))
    assert len(set(concurrent)) == 12
    assert [p.read_bytes() for p in concurrent] == [str(i).encode() for i in range(12)]


@pytest.mark.parametrize("name", ["../recording.csv", r"C:\fakepath\recording.csv"])
def test_upload_path_components_cannot_escape(tmp_path, name):
    path = save_uploaded_bytes(name, b"valid bytes", tmp_path)
    assert path.name == "recording.csv"
    assert path.resolve().is_relative_to(tmp_path.resolve())


@pytest.mark.parametrize("name,payload", [("", b"x"), ("..", b"x"), ("file.exe", b"x"), ("file.csv", b""), ("file.csv:stream", b"x"), ("CON.csv", b"x")])
def test_invalid_uploads_are_rejected_without_files(tmp_path, name, payload):
    with pytest.raises(InputValidationError):
        save_uploaded_bytes(name, payload, tmp_path)
    assert not list(tmp_path.iterdir())


def test_streamlit_adapter_uses_unique_storage(tmp_path, monkeypatch):
    app_path = Path(__file__).parents[1] / "app/bsd_workstation.py"
    spec = importlib.util.spec_from_file_location("bsd_upload_adapter_test", app_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "PROJECT_ROOT", tmp_path)
    target = module._save_upload_to_local_tmp(SimpleNamespace(name="same.csv", getvalue=lambda: b"TARGET"))
    background = module._save_upload_to_local_tmp(SimpleNamespace(name="same.csv", getvalue=lambda: b"BACKGROUND"))
    assert target != background
    assert target.read_bytes() == b"TARGET"
    assert background.read_bytes() == b"BACKGROUND"


@pytest.mark.parametrize("text,match", [
    ("", "empty"), ("\n1\n2\n", "blank"), ('""\n1\n2\n', "blank"),
    ("value\n1\n\n2\n", "row 3"), ("value\n1\nnope\n", "row 3"),
    ("value\n1\nNaN\n", "non-finite"), ("value\n1\nInfinity\n", "non-finite"),
    ("time,value\n0,1\n1\n", "row 3"),
])
def test_csv_bad_rows_are_not_dropped(tmp_path, text, match):
    path = tmp_path / "bad.csv"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(InputValidationError, match=match):
        load_timeseries_file(path, sample_rate_hz=1024.)


@pytest.mark.parametrize("column", [-1, True, 4, "missing"])
def test_csv_invalid_column_is_explicit(tmp_path, column):
    path = tmp_path / "columns.csv"
    path.write_text("value\n1\n2\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="column|Column"):
        load_timeseries_file(path, sample_rate_hz=1024., csv_value_column=column)


def test_csv_time_metadata_and_rate_validation(tmp_path):
    path = tmp_path / "timed.csv"
    times = np.arange(4096) / 1024.
    path.write_text("time,value\n" + "\n".join(f"{t:.12f},{np.sin(i * .1):.12f}" for i, t in enumerate(times)), encoding="utf-8")
    ts = load_timeseries_file(path, csv_time_column="time", csv_value_column="value", sample_rate_hz=1024.)
    assert ts.n_samples == 4096
    assert ts.sample_rate_hz == pytest.approx(1024.)
    assert ts.signal[2] == pytest.approx(np.sin(.2))
    with pytest.raises(InputValidationError, match="conflicts"):
        load_timeseries_file(path, csv_time_column="time", csv_value_column="value", sample_rate_hz=2048.)
    with pytest.raises(InputValidationError, match="different"):
        load_timeseries_file(path, csv_time_column=0, csv_value_column=0)
    path.write_text("time,value\n0,1\n.001,2\n.003,3\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="non-uniform"):
        load_timeseries_file(path, csv_time_column=0, csv_value_column=1)


def test_pickle_npy_and_invalid_encoding_are_validation_errors(tmp_path):
    path = tmp_path / "object.npy"
    np.save(path, np.array([{}], dtype=object))
    with pytest.raises(InputValidationError, match="safe numeric NPY"):
        load_timeseries_file(path, sample_rate_hz=1024.)
    path = tmp_path / "encoding.csv"
    path.write_bytes(b"\xff\xfe\x00")
    with pytest.raises(InputValidationError, match="Failed to read CSV"):
        load_timeseries_file(path, sample_rate_hz=1024.)


def test_matching_rates_and_disjoint_intervals_are_required(tmp_path):
    path = tmp_path / "data.npy"
    np.save(path, np.sin(np.arange(4096) * .1))
    target = load_timeseries_file(path, sample_rate_hz=1024.)
    background = replace(target, sha256="different", time_interval=(0., 4.), time_reference="GPS", provenance={"detector": "L1"})
    target = replace(target, time_interval=(2., 6.), time_reference="GPS", provenance={"detector": "H1"})
    with pytest.raises(InputValidationError, match="overlaps"):
        check_background_independence(target, background)
    background = replace(background, time_interval=(10., 14.))
    assert check_background_independence(target, background).status == "VERIFIED_DISJOINT_TIMESTAMPS"
    for invalid_rate in [512., float("nan"), float("inf")]:
        with pytest.raises(InputValidationError, match="sample rates"):
            check_background_independence(target, replace(background, sample_rate_hz=invalid_rate))
    assert check_background_independence(replace(target, time_interval=None), background).status == "UNVERIFIED_EXTERNAL_CHECK_REQUIRED"


def test_preserved_upload_basename_still_triggers_guard(tmp_path):
    path = save_uploaded_bytes("TEST_RESERVED_EVENT.csv", b"1\n2\n", tmp_path)
    with pytest.raises(io_loader.O3HoldoutBlockedError):
        load_timeseries_file(path, sample_rate_hz=1024.)


@pytest.mark.skipif(shutil.which("node") is None, reason="Node.js is required for browser input regression checks")
def test_browser_input_regressions():
    result = subprocess.run(["node", "--test", str(Path(__file__).with_name("test_browser_input_validation.js"))], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("array", [np.array([1+2j, 3+4j]), np.array(["1", "2"]), np.zeros(2, dtype=[("value", "f8")])])
def test_npy_rejects_non_real_samples_without_coercion(tmp_path, array):
    path = tmp_path / "unsupported.npy"
    np.save(path, array)
    with pytest.raises(InputValidationError, match="real numeric samples"):
        load_timeseries_file(path, sample_rate_hz=1024.)


def test_streamlit_rejected_input_clears_previous_run(tmp_path, monkeypatch):
    from contextlib import nullcontext

    app_path = Path(__file__).parents[1] / "app/bsd_workstation.py"
    spec = importlib.util.spec_from_file_location("bsd_rejected_input_test", app_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "PROJECT_ROOT", tmp_path)

    class FakeStreamlit:
        def __init__(self):
            self.session_state = {"bsd_last_result": object()}
            self.sidebar = nullcontext()
            self.errors = []
        def radio(self, *args, **kwargs): return "Local File Analysis"
        def selectbox(self, label, options, **kwargs): return options[0]
        def number_input(self, label, *args, **kwargs): return args[2]
        def checkbox(self, *args, **kwargs): return False
        def text_input(self, *args, **kwargs): return "arbitrary"
        def button(self, *args, **kwargs): return True
        def expander(self, *args, **kwargs): return nullcontext()
        def spinner(self, *args, **kwargs): return nullcontext()
        def file_uploader(self, label, **kwargs):
            payload = b"value\n1\nbad\n" if label == "Target Dataset" else b"value\n1\n2\n"
            return SimpleNamespace(name="same.csv", getvalue=lambda: payload)
        def error(self, message): self.errors.append(message)
        def tabs(self, *args, **kwargs): raise AssertionError("Stale result must not be rendered")
        def __getattr__(self, name): return lambda *args, **kwargs: None

    fake = FakeStreamlit()
    monkeypatch.setattr(module, "st", fake)
    module.main()
    assert fake.session_state.get("bsd_last_result") is None
    assert len(fake.errors) == 1 and "row 3" in fake.errors[0]
    stored = list((tmp_path / ".tmp/uploads").glob("*/same.csv"))
    assert len(stored) == 2
    assert {p.read_bytes() for p in stored} == {b"value\n1\nbad\n", b"value\n1\n2\n"}
