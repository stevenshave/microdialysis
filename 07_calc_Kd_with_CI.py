"""
Perform microdialysis KD determination with confidence intervals

Use LMFit to calculate 95% CI for KD assignment
"""
635
import numpy as np
PROTEIN_CONC=80
from microdialysis_equations import *
from itertools import permutations


L0_CONC=25
REDVOL=100
WHITEVOL=300


pc=1.02
pc_std=0.01

pt=1.21
pt_std=0.01



dist_pt=[pt-pt_std,pt,pt+pt_std]
dist_pc=[pc-pc_std,pc,pc+pc_std]

        
pvalues=[(pc_val,pt_val) for pc_val in dist_pc for pt_val in dist_pt]
kds=[qud_Kd_from_ptvalue(x[1],PROTEIN_CONC,L0_CONC,REDVOL,WHITEVOL,x[0]) for x in pvalues]
print(qud_Kd_from_ptvalue(pt,PROTEIN_CONC,L0_CONC,REDVOL,WHITEVOL,pc))
print(list(kds))
print(np.mean(kds), np.std(kds))