"""
Differentiate redconc, whiteconc, ptvalue

Simulate the rate of change of redconc, whiteconc and ptvalues
Requires the autograd package

"""
import sys
from autograd import grad
from microdialysis_equations import *
import numpy as np
T0 = 80.0
L0_CONC = 50.0
REDVOL = 100.0
WHITEVOL = 300.0
PCVALUE = 1.0

KDS_to_simulate = [1., 100., 200., 300., 400., 500.]


g_lred_wrt_kd = grad(lambda kd: qud_lred(T0,L0_CONC,kd,REDVOL,WHITEVOL,1.0))
g_lwhite_wrt_kd = grad(lambda kd: qud_lwhite(T0,L0_CONC,kd,REDVOL,WHITEVOL,1.0))
g_ptvalue_wrt_kd = grad(lambda kd: qud_ptvalue(T0,L0_CONC,kd,REDVOL,WHITEVOL,1.0))

print(f"{'KD':>10},{'lred':>10},{'dlreddKD':>10},{'lwhite':>10},{'dlwhitedKD':>10},{'ptvalue':>10},{'dlptdKD':>10}")
for kd in KDS_to_simulate:
    lwhite = qud_lwhite(T0, L0_CONC, kd, REDVOL, WHITEVOL, 1.0)
    lred = qud_lred(T0, L0_CONC, kd, REDVOL, WHITEVOL, 1.0)
    pt = qud_ptvalue(T0, L0_CONC, kd, REDVOL, WHITEVOL, 1.0)
    print(f"{kd:>10.0f},{lred:>10.4f},{g_lred_wrt_kd(kd):>10.4f},{lwhite:>10.4f},{g_lwhite_wrt_kd(kd):>10.4f},{pt:>10.4f},{g_ptvalue_wrt_kd(kd):>10.4f}")
print()


g_kd_wrt_lred = grad(lambda lred: qud_Kd_from_lred(lred,T0,L0_CONC,REDVOL,WHITEVOL,1.0))
g_kd_wrt_lwhite = grad(lambda lwhite: qud_Kd_from_lwhite(lwhite,T0,L0_CONC,REDVOL,WHITEVOL,1.0))
g_kd_wrt_ptval = grad(lambda pt: qud_Kd_from_ptvalue(pt,T0,L0_CONC,REDVOL,WHITEVOL,1.0))
print(f"{'KD':>10},{'lred':>10},{'dKDdlred':>10},{'lwhite':>10},{'dlKDdlwhite':>10},{'ptvalue':>10},{'dKDdlpt':>10}")
for kd in KDS_to_simulate:
    lred = qud_lred(T0, L0_CONC, kd, REDVOL, WHITEVOL, 1.0)
    lwhite = qud_lwhite(T0, L0_CONC, kd, REDVOL, WHITEVOL, 1.0)
    pt = qud_ptvalue(T0, L0_CONC, kd, REDVOL, WHITEVOL, 1.0)
    print(f"{kd:>10.0f},{lred:>10.4f},{g_kd_wrt_lred(lred):>10.4f},{lwhite:>10.4f},{g_kd_wrt_lwhite(lwhite):>10.4f},{pt:>10.4f},{g_kd_wrt_ptval(pt):>10.4f}")
