"""Independent SciPy reference for the browser's fixed-parameter WAV FIR."""
from __future__ import annotations

import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from scipy.signal import firwin, resample_poly


@pytest.mark.skipif(shutil.which("node") is None, reason="Node.js is required for browser resampler reference checks")
@pytest.mark.parametrize("source_rate", [1024, 44100, 48000])
def test_browser_wav_matches_independent_scipy_fir(source_rate):
    # Float32 is the Web Audio decoder's sample representation. No real files or O3 data.
    rng=np.random.default_rng(78)
    samples=rng.normal(size=source_rate//4).astype(np.float32).astype(np.float64)
    script=Path(__file__).with_name("run_browser_wav_resampler.js")
    process=subprocess.run(["node",str(script)],input=json.dumps({"source_rate":source_rate,"signal":samples.tolist()}),
                           text=True,capture_output=True,timeout=60)
    assert process.returncode==0,process.stderr
    decoded=json.loads(process.stdout)
    observed=np.array(decoded["signal"],dtype=np.float64)
    common=math.gcd(source_rate,512);up,down=512//common,source_rate//common
    half=math.ceil(24*source_rate/512)
    # FIR designed independently by SciPy on the upsampled grid. Its global DC
    # normalization differs slightly from the browser's per-phase normalization.
    coefficients=firwin(2*half*up+1,224.,fs=source_rate*up,window=("kaiser",8.6))
    expected=resample_poly(samples,up,down,window=coefficients,padtype="constant")[:len(observed)]
    assert decoded["fs"]==512
    assert len(observed)==len(samples)*512//source_rate
    assert np.all(np.isfinite(observed))
    assert np.max(np.abs(observed-expected))<1e-5
    assert decoded["preprocessing"]["padding"]=="zero"
