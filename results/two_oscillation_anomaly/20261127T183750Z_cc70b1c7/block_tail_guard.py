"""Experimental whole-block guard; never edits a core calibration or deletes samples."""
import numpy as np

def block_guard(scores):
    z=np.asarray(scores,dtype=float)
    if z.ndim!=1 or z.size<40 or not np.all(np.isfinite(z)):
        return {"admitted":False,"reason":"INVALID_OR_INSUFFICIENT_NULL","median":None,"spread":None,"limit":None,"maximum":None}
    median=float(np.median(z));spread=max(float(np.quantile(z,.95))-median,1.)
    limit=median+10.*spread;maximum=float(np.max(z))
    if not np.isfinite(limit):
        return {"admitted":False,"reason":"NONFINITE_GUARD_LIMIT","median":median,"spread":spread,"limit":None,"maximum":maximum}
    admitted=maximum<=limit
    return {"admitted":bool(admitted),"reason":"ADMITTED" if admitted else "WHOLE_BLOCK_TAIL_QUARANTINE","median":median,"spread":spread,"limit":limit,"maximum":maximum}

def focused_checks():
    checks=[]
    def check(name,ok):
        assert ok,name
        checks.append(name)
    check("constant_ties",block_guard(np.full(40,3.))["admitted"])
    a=np.r_[np.zeros(99),10.]
    check("strict_boundary_equal",block_guard(a)["admitted"])
    check("nextafter_above",not block_guard(np.r_[np.zeros(99),np.nextafter(10.,np.inf)])["admitted"])
    check("nextafter_below",block_guard(np.r_[np.zeros(99),np.nextafter(10.,0.)])["admitted"])
    for name,z in [("empty",[]),("small",np.zeros(39)),("nan",np.r_[np.zeros(39),np.nan]),("infinity",np.r_[np.zeros(39),np.inf]),("matrix",np.zeros((40,2)))]:
        check(name,block_guard(z)["reason"]=="INVALID_OR_INSUFFICIENT_NULL")
    check("tail_spike",not block_guard(np.r_[np.zeros(99),100.])["admitted"])
    from bsd.calibration import conformal_order_rank
    check("fap005_attainable",conformal_order_rank(199,.005)[1])
    union=-np.expm1(np.log1p(-.005)/2.)
    check("union005_unattainable",not conformal_order_rank(199,float(union))[1])
    return {"status":"PASS","checks":checks,"count":len(checks)}
