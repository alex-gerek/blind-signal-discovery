"""Experiment-only timing/state eligibility. It does not infer background purity."""
import math,json
POLICIES=("frozen_naive","age_only","age_and_state")
def eligibility(policy,start,end,source_start,source_end,source_state,target_state,max_age=2400.,clock="SIM",source_clock="SIM",ready=True):
    if policy not in POLICIES:raise ValueError("Unknown policy")
    if not ready:return False,"CALIBRATION_NOT_READY"
    if clock is None or source_clock is None or clock!=source_clock:return False,"CLOCK_UNVERIFIED"
    if any(x is None or not math.isfinite(x) for x in [start,end,source_start,source_end,max_age]):return False,"TIMING_UNVERIFIED"
    if end<=start or source_end<=source_start or max_age<=0:return False,"INVALID_INTERVAL_OR_POLICY"
    if source_end>start:return False,"SOURCE_TARGET_OVERLAP"
    if policy=="frozen_naive":return True,"NAIVE_REUSE_UNQUALIFIED"
    if policy=="age_and_state":
        if source_state is None or target_state is None:return False,"STATE_UNVERIFIED"
        if source_state!=target_state:return False,"STATE_MISMATCH"
    if end-source_start>max_age:return False,"AGE_EXPIRED"
    return True,"ELIGIBLE_UNVALIDATED_POLICY"

def interval_exposure(records):
    intervals=[];seen={}
    for ident,a,b in records:
        if not(math.isfinite(a) and math.isfinite(b)) or b<a:raise ValueError("Invalid interval")
        if ident in seen and seen[ident]!=(a,b):raise ValueError("Conflicting duplicate ID")
        seen[ident]=(a,b)
        if b>a:intervals.append((a,b))
    total=0.;last=None
    for a,b in sorted(intervals):
        if last is None or a>=last:total+=b-a;last=b
        elif b>last:total+=b-last;last=b
    return total

def focused_checks():
    import numpy as np
    from bsd.calibration import conformal_order_rank
    checks=[]
    def check(name,value):assert value,name;checks.append(name)
    def e(**kw):
        args=dict(policy="age_and_state",start=2396.,end=2400.,source_start=0.,source_end=1136.,source_state="A",target_state="A")
        args.update(kw);return eligibility(**args)
    check("age_boundary_equal",e()[0])
    check("age_nextafter_above",e(end=float(np.nextafter(2400.,np.inf)))[1]=="AGE_EXPIRED")
    check("age_nextafter_below",e(end=float(np.nextafter(2400.,0.)))[0])
    check("state_changed",e(target_state="B")[1]=="STATE_MISMATCH")
    check("state_unknown",e(target_state=None)[1]=="STATE_UNVERIFIED")
    check("clock_unknown",e(clock=None)[1]=="CLOCK_UNVERIFIED")
    check("clock_mismatch",e(clock="OTHER")[1]=="CLOCK_UNVERIFIED")
    check("time_unknown",e(end=None)[1]=="TIMING_UNVERIFIED")
    check("nonfinite",e(end=float('inf'))[1]=="TIMING_UNVERIFIED")
    check("empty_source",e(source_end=0.)[1]=="INVALID_INTERVAL_OR_POLICY")
    check("source_overlap",e(source_end=2397.)[1]=="SOURCE_TARGET_OVERLAP")
    check("source_touches_target",e(source_end=2396.)[0])
    check("not_ready",e(ready=False)[1]=="CALIBRATION_NOT_READY")
    check("age_only_ignores_state",e(policy="age_only",target_state="B")[0])
    check("frozen_ignores_age_and_state",e(policy="frozen_naive",end=2600.,target_state="B")[0])
    check("empty_exposure",interval_exposure([])==0.)
    check("zero_duplicate_overlap_gap",interval_exposure([('a',0.,0.),('b',0.,4.),('b',0.,4.),('c',2.,6.),('d',8.,10.)])==8.)
    for name,rows in [('reversed',[('x',2.,1.)]),('duplicate_conflict',[('x',0.,1.),('x',0.,2.)])]:
        try:interval_exposure(rows)
        except ValueError:checks.append(name)
        else:raise AssertionError(name)
    check("fap005_199attainable",conformal_order_rank(199,.005)[1])
    check("fap005_198insufficient",not conformal_order_rank(198,.005)[1])
    check("union005_199insufficient",not conformal_order_rank(199,float(-np.expm1(np.log1p(-.005)/2)))[1])
    fixture={'available':False,'threshold':None,'p':None,'diagnostic_candidates':[{'score':20.}]}
    check("strict_json_diagnostic_preservation",json.loads(json.dumps(fixture,allow_nan=False))==fixture)
    return {'status':'PASS','checks':checks,'count':len(checks),'scope':'Experimental controller/ledger fixtures, not production UI tests'}
