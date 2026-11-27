"""Descriptive warning accounting, separate from monitor decisions."""
def episodes(flags,offset=284):
    result=[];start=None
    for i,flag in enumerate(list(flags)+[False]):
        if flag and start is None:start=i
        if not flag and start is not None:
            result.append({'start_index':start+offset,'end_index_exclusive':i+offset,'windows':i-start,'duration_s':4.*(i-start)});start=None
    return result

def summarize(records,assessments,search,stable_flags=None):
    flags=[bool(a['regime_shifted']) for a in assessments];ep=episodes(flags)
    eligible=[a['credibility_status']!='WARMUP_INSUFFICIENT_HISTORY' for a in assessments]
    signals=[i for i,r in enumerate(records) if r['signal_kind'] is not None]
    first_post=next((i+284 for i,v in enumerate(flags) if i+284>=360 and v),None)
    first_new=next((e['start_index'] for e in ep if e['start_index']>=360),None)
    summary={'target_windows':len(records),'monitor_eligible_windows':sum(eligible),'warmup_windows':sum(not x for x in eligible),'warning_windows':sum(flags),'warning_episodes':ep,'any_warning':any(flags),
      'prechange_warning_windows':sum(flags[:76]),'warning_active_before_change':flags[75],
      'first_post_change_warning_index':first_post,'first_post_change_warning_delay_s':4.*(first_post+1-360) if first_post is not None else None,
      'first_new_episode_after_change':first_new,'first_new_episode_delay_s':4.*(first_new+1-360) if first_new is not None else None,
      'trigger_windows':{name:sum(bool(a.get('diagnostics',{}).get(name,False)) for a in assessments) for name in ['ks_triggered','z_triggered','cusum_triggered']},
      'signal_windows':len(signals),'signal_warning_windows':sum(flags[i] for i in signals),'signal_warmup_windows':sum(not eligible[i] for i in signals),
      'associated_search_signal_windows':sum(search[i]['associated'] for i in signals),'associated_recoveries_that_warning_gate_would_remove':sum(search[i]['associated'] and flags[i] for i in signals),
      'search_noise_alarms':sum(s['detected'] for s,r in zip(search,records) if r['signal_kind'] is None),'search_noise_windows':sum(r['signal_kind'] is None for r in records),
      'detector_available_windows':sum(s['available'] for s in search),'background_updated':False,'candidates_suppressed':False}
    if stable_flags is not None:
        summary['extra_warning_windows_vs_stable']=sum(a and not b for a,b in zip(flags,stable_flags))
        summary['suppressed_warning_windows_vs_stable']=sum(b and not a for a,b in zip(flags,stable_flags))
        # Paired changes after a transient can persist through the running CUSUM.
        last=max(signals) if signals else None
        summary['extra_warning_windows_after_last_signal']=sum(flags[i] and not stable_flags[i] for i in range(last+1,len(flags))) if last is not None else 0
    return summary
