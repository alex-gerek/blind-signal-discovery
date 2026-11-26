def _gen_obs_list(rng: np.random.Generator, n_win: int, noise_model: str, sigma_fn: Any, prefix: str) -> list[Observation]:
    out: list[Observation] = []
    for i in range(n_win):
        x, _ = generate_noise(noise_model, 4096, 1024.0, rng, sigma=float(sigma_fn(i)))
        out.append(Observation(f'{prefix}_{i:03d}', x, 1024.0, 4.0))
    return out

def _make_chirp_chunks(n_windows: int, f0: float, f1: float, snr_per_win: float) -> list[np.ndarray]:
    dur = n_windows * 4.0
    t = np.arange(n_windows * 4096, dtype=np.float64) / 1024.0
    slope = (f1 - f0) / dur
    w = np.sin(2.0 * np.pi * (f0 * t + 0.5 * slope * t ** 2))
    return [_scale_to_energy_snr(w[k * 4096:(k + 1) * 4096], 1.0, snr_per_win) for k in range(n_windows)]

def _wilson_ci90(k: int, n: int) -> tuple[float, float]:
    if n <= 0:
        return (0.0, 1.0)
    z = 1.6448536
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * np.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / denom
    return (max(0.0, float(center - half)), min(1.0, float(center + half)))

def _run_suite_a_corrected() -> dict[str, Any]:
    """Identical seed (3026) and scenarios to 8ed1813, with Corrections 1 & 2."""
    rng = make_split_rng(3026, 'v03_phase2b_suite_a')
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    print('[v0.3-P2b-Corr] Fitting stationary reference R_0 (pink_colored, N=320)...')
    ref_obs = _gen_obs_list(rng, 320, 'pink_colored', lambda _: 1.0, 'ref')
    scal = calibrate_separate_search_modes(det, bb, ref_obs, target_fap_transient=0.01, target_fap_persistent=0.01, target_fwer_union=0.01, norm_fraction=0.3, calibration_seed=3026)
    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode='causal_trailing', policy='transient'))
    cref = mgr.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)
    u0 = scal.transient_stft_uncertainty
    rows: list[dict[str, Any]] = []
    print('[v0.3-P2b-Corr]   Suite A1: stationary_null_200win...')
    stat_obs = _gen_obs_list(rng, 200, 'pink_colored', lambda _: 1.0, 'stat')
    out_stat = mgr.process_sequence(stat_obs, det, bb, scal, cref, gref)
    tr_fp = sum((1 for o in out_stat if o.stft_nominated_post_veto))
    ps_fp = sum((1 for o in out_stat if o.broadband_nominated))
    un_fp = sum((1 for o in out_stat if o.joint_nominated_post_veto))
    ci_tr = _wilson_ci90(tr_fp, len(out_stat))
    ci_ps = _wilson_ci90(ps_fp, len(out_stat))
    ci_un = _wilson_ci90(un_fp, len(out_stat))
    rows.append({'scenario': 'stationary_null_200win', 'n_windows': len(out_stat), 'duty_cycle_clean': float(np.mean([o.is_clean_calibrated for o in out_stat])), 'transient_mode_far': tr_fp / len(out_stat), 'transient_far_ci90': f'[{ci_tr[0] * 100:.2f}%, {ci_tr[1] * 100:.2f}%]', 'persistent_mode_far': ps_fp / len(out_stat), 'persistent_far_ci90': f'[{ci_ps[0] * 100:.2f}%, {ci_ps[1] * 100:.2f}%]', 'union_sidak_far': un_fp / len(out_stat), 'union_far_ci90': f'[{ci_un[0] * 100:.2f}%, {ci_un[1] * 100:.2f}%]', 'cal_n_thresh': u0.n_thresh, 'cal_n_eff': u0.n_eff, 'cal_ljung_box_p': u0.ljung_box_p, 'cal_beta_ci90_neff_hi': u0.far_ci90_neff[1], 't_provisional_win': None, 't_mature_win': None, 'mature_surv_far': None, 'transient_eff_full': None, 'transient_eff_clean_cal': None, 'glitch_storm_flag_rate': float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_stat])), 'n_regimes': max((o.regime_id for o in out_stat)) + 1, 'glitch_ref_refreshes': out_stat[-1].diagnostics['glitch_ref_refresh_count']})
    print('[v0.3-P2b-Corr]   Suite A2: abrupt_step_500win...')
    ab_obs = _gen_obs_list(rng, 100, 'pink_colored', lambda _: 1.0, 'ab0') + _gen_obs_list(rng, 400, 'shifted_ar2', lambda _: 1.5, 'ab1')
    trans_idx = [108, 125, 145]
    post_idx = list(range(330, 430, 10))
    is_sig_ab = [False] * 500
    for idx in trans_idx + post_idx:
        clean, *_ = generate_clean_signal('linear_chirp', 4096, 1024.0, rng)
        s = _scale_to_energy_snr(clean, 1.5, 18.0)
        ab_obs[idx] = Observation(ab_obs[idx].sample_id, ab_obs[idx].signal + s, 1024.0, 4.0)
        is_sig_ab[idx] = True
    out_ab = mgr.process_sequence(ab_obs, det, bb, scal, cref, gref)
    t_prov = next((o.index - 100 for o in out_ab[100:] if o.regime_id >= 1 and o.credibility_status == CRED_CALIBRATED), None)
    t_mat = next((o.index - 100 for o in out_ab[100:] if o.regime_id >= 1 and o.diagnostics.get('is_mature')), None)
    mature_bg = [o for o, s in zip(out_ab, is_sig_ab) if not s and o.regime_id >= 1 and o.diagnostics.get('is_mature') and o.is_clean_calibrated]
    mature_fp = sum((1 for o in mature_bg if o.stft_nominated_post_veto))
    mature_far = mature_fp / len(mature_bg) if mature_bg else None
    post_eff_full = float(np.mean([out_ab[i].stft_nominated_post_veto for i in post_idx]))
    post_eff_clean = float(np.mean([out_ab[i].stft_nominated_post_veto and out_ab[i].is_clean_calibrated for i in post_idx]))
    trans_pres = float(np.mean([out_ab[i].stft_nominated_raw and out_ab[i].significance_status == SIG_UNCERTAIN for i in trans_idx]))
    u_mat = out_ab[-1].calibration_uncertainty
    post_rec_glitch_rate = float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_ab[170:]]))
    rows.append({'scenario': 'abrupt_step_500win', 'n_windows': len(out_ab), 'duty_cycle_clean': float(np.mean([o.is_clean_calibrated for o in out_ab])), 'transient_mode_far': float(np.mean([o.stft_nominated_post_veto for o, s in zip(out_ab, is_sig_ab) if not s])), 'transient_far_ci90': '', 'persistent_mode_far': None, 'persistent_far_ci90': '', 'union_sidak_far': None, 'union_far_ci90': '', 'cal_n_thresh': u_mat.n_thresh, 'cal_n_eff': u_mat.n_eff, 'cal_ljung_box_p': u_mat.ljung_box_p, 'cal_beta_ci90_neff_hi': u_mat.far_ci90_neff[1], 't_provisional_win': t_prov, 't_mature_win': t_mat, 'mature_surv_far': mature_far, 'transient_eff_full': post_eff_full, 'transient_eff_clean_cal': post_eff_clean, 'glitch_storm_flag_rate': float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_ab])), 'post_recovery_glitch_flag_rate': post_rec_glitch_rate, 'n_regimes': max((o.regime_id for o in out_ab)) + 1, 'glitch_ref_refreshes': out_ab[-1].diagnostics['glitch_ref_refresh_count'], 'transition_raw_preservation': trans_pres})
    print('[v0.3-P2b-Corr]   Suite A3: transient_chirps_stationary...')
    tr_obs = _gen_obs_list(rng, 120, 'pink_colored', lambda _: 1.0, 'tr')
    tr_idx = list(range(15, 115, 5))
    tr_sig = [False] * 120
    for idx in tr_idx:
        clean, *_ = generate_clean_signal('linear_chirp', 4096, 1024.0, rng)
        s = _scale_to_energy_snr(clean, 1.0, 18.0)
        tr_obs[idx] = Observation(tr_obs[idx].sample_id, tr_obs[idx].signal + s, 1024.0, 4.0)
        tr_sig[idx] = True
    out_tr = mgr.process_sequence(tr_obs, det, bb, scal, cref, gref)
    rows.append({'scenario': 'transient_chirps_stationary', 'n_windows': 120, 'duty_cycle_clean': float(np.mean([o.is_clean_calibrated for o in out_tr])), 'transient_mode_far': float(np.mean([o.stft_nominated_post_veto for o, s in zip(out_tr, tr_sig) if not s])), 'persistent_mode_far': float(np.mean([o.broadband_nominated for o, s in zip(out_tr, tr_sig) if not s])), 'union_sidak_far': float(np.mean([o.joint_nominated_post_veto for o, s in zip(out_tr, tr_sig) if not s])), 'cal_n_thresh': u0.n_thresh, 'cal_n_eff': u0.n_eff, 'cal_ljung_box_p': u0.ljung_box_p, 'cal_beta_ci90_neff_hi': u0.far_ci90_neff[1], 'transient_eff_full': float(np.mean([out_tr[i].stft_nominated_post_veto for i in tr_idx])), 'transient_eff_clean_cal': float(np.mean([out_tr[i].stft_nominated_post_veto and out_tr[i].is_clean_calibrated for i in tr_idx])), 'glitch_storm_flag_rate': float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_tr])), 'n_regimes': max((o.regime_id for o in out_tr)) + 1, 'glitch_ref_refreshes': out_tr[-1].diagnostics['glitch_ref_refresh_count']})
    print('[v0.3-P2b-Corr]   Suite A4: glitch_storm_120win...')
    gl_obs = _gen_obs_list(rng, 60, 'pink_colored', lambda _: 1.0, 'gl0') + _gen_obs_list(rng, 60, 'impulsive', lambda _: 1.0, 'gl1')
    out_gl = mgr.process_sequence(gl_obs, det, bb, scal, cref, gref)
    gl_post_flag = float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_gl[60:]]))
    rows.append({'scenario': 'glitch_storm_120win', 'n_windows': 120, 'duty_cycle_clean': float(np.mean([o.is_clean_calibrated for o in out_gl])), 'transient_mode_far': float(np.mean([o.stft_nominated_post_veto for o in out_gl])), 'persistent_mode_far': float(np.mean([o.broadband_nominated for o in out_gl])), 'union_sidak_far': float(np.mean([o.joint_nominated_post_veto for o in out_gl])), 'cal_n_thresh': u0.n_thresh, 'cal_n_eff': u0.n_eff, 'cal_ljung_box_p': u0.ljung_box_p, 'cal_beta_ci90_neff_hi': u0.far_ci90_neff[1], 'mature_surv_far': float(np.mean([o.stft_nominated_post_veto for o in out_gl if o.is_clean_calibrated])), 'glitch_storm_flag_rate': gl_post_flag, 'n_regimes': max((o.regime_id for o in out_gl)) + 1, 'glitch_ref_refreshes': out_gl[-1].diagnostics['glitch_ref_refresh_count']})
    ab_tl = {'z_drift': [o.max_feature_z_drift for o in out_ab], 'glitch_flag': [int(o.glitch_storm_status == GLITCH_ELEVATED) for o in out_ab], 'stage': [0 if o.credibility_status == CRED_CALIBRATED and o.regime_id == 0 else 1 if o.credibility_status in (CRED_TRANSITION, CRED_RECALIBRATING) else 2 if not o.diagnostics.get('is_mature') else 3 for o in out_ab], 'post_idx': post_idx, 'post_hits': [int(out_ab[i].stft_nominated_post_veto and out_ab[i].is_clean_calibrated) for i in post_idx]}
    return {'rows': rows, 'ab_tl': ab_tl, 'shared_ref': (det, bb, scal, cref, gref)}

def _run_suite_b_corrected(shared_ref: tuple[Any, ...]) -> dict[str, Any]:
    """Identical persistent matrix (seed 3027) under dedicated Persistent Search mode."""
    det, bb, scal, cref, gref = shared_ref
    rng = make_split_rng(3027, 'v03_phase2b_suite_b')
    mgr_pers = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode='causal_trailing', policy='persistent'))
    t_win = np.arange(4096, dtype=np.float64) / 1024.0
    sos_bb = sp_signal.butter(4, [80.0, 300.0], btype='bandpass', fs=1024.0, output='sos')
    rows: list[dict[str, Any]] = []
    for bg_mode, sigma_fn in (('stationary', lambda _: 1.0), ('drifting_1.0_to_1.4x', lambda i: 1.0 + 0.4 * (i / 119.0))):
        for morph in ('80s_sinusoid', '80s_multitone_comb', '40s_slow_chirp', '80s_broadband_stochastic'):
            print(f'[v0.3-P2b-Corr]   Suite B (persistent): {morph} | bg={bg_mode}...')
            obs = _gen_obs_list(rng, 120, 'pink_colored', sigma_fn, f'p_{morph[:3]}')
            is_sig = [False] * 120
            if morph == '80s_sinusoid':
                for idx in range(50, 70):
                    s = _scale_to_energy_snr(np.sin(2.0 * np.pi * 180.0 * (t_win + idx * 4.0)), 1.0, 20.0)
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + s, 1024.0, 4.0)
                    is_sig[idx] = True
            elif morph == '80s_multitone_comb':
                for idx in range(50, 70):
                    tt = t_win + idx * 4.0
                    comb = sum((np.sin(2.0 * np.pi * f * tt) for f in (90.0, 140.0, 190.0, 240.0, 290.0)))
                    s = _scale_to_energy_snr(comb, 1.0, 22.0)
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + s, 1024.0, 4.0)
                    is_sig[idx] = True
            elif morph == '40s_slow_chirp':
                for k, ch in enumerate(_make_chirp_chunks(10, 60.0, 260.0, 16.0)):
                    idx = 50 + k
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + ch, 1024.0, 4.0)
                    is_sig[idx] = True
            else:
                for idx in range(50, 70):
                    raw_bb = sp_signal.sosfiltfilt(sos_bb, rng.normal(0.0, 1.0, size=4096))
                    s = _scale_to_energy_snr(raw_bb, 1.0, 26.0)
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + s, 1024.0, 4.0)
                    is_sig[idx] = True
            out = mgr_pers.process_sequence(obs, det, bb, scal, cref, gref)
            sig_wins = [o for o, s in zip(out, is_sig) if s]
            bg_wins = [o for o, s in zip(out, is_sig) if not s]
            rows.append({'morphology': morph, 'background': bg_mode, 'transient_stft_eff_full': float(np.mean([o.stft_nominated_post_veto for o in sig_wins])), 'persistent_bb_eff_full': float(np.mean([o.broadband_nominated for o in sig_wins])), 'union_sidak_eff_full': float(np.mean([o.joint_nominated_post_veto for o in sig_wins])), 'clean_cal_stft_eff': float(np.mean([o.stft_nominated_post_veto and o.is_clean_calibrated for o in sig_wins])), 'bg_union_far_full': float(np.mean([o.joint_nominated_post_veto for o in bg_wins])), 'n_regimes': max((o.regime_id for o in out)) + 1})
    return {'rows': rows}