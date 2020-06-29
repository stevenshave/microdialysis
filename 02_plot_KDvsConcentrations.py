"""
Simulate microdialysis experiment behaviour (KD vs Conc)

Simulate KD vs compound concentrations in the red and white chambers
of a microdialysis experiment
"""

from matplotlib import pyplot as plt
import numpy as np
t0=80
from microdialysis_equations import *

l0=50
redvol=100
whitevol=300
KD_beginning = 0
KD_end = 500
pc = 1.0
NUM_POINTS_ON_XAXIS = 1000 # Publication used 1000 pts along X
x_axis = np.linspace(KD_beginning,KD_end, NUM_POINTS_ON_XAXIS)

y=np.full((2,NUM_POINTS_ON_XAXIS), np.nan)
y[0]=qud_lred(t0, l0, x_axis, redvol, whitevol, pc)
y[1]=qud_lwhite(t0, l0, x_axis, redvol, whitevol, pc)

fig, ax = plt.subplots(1,1, figsize=(7.204724, 5.09424929292), sharex=True)
#fig, ax = plt.subplots(1,1, figsize=(8, 6), sharex=True)
fig.suptitle("qµD simulation",y=0.95, fontsize=14)#, K$_\mathrm{D}$ vs red and white chamber compound concentration."+
#f"\nRed chamber = {redvol} µl, White chamber = {whitevol} µl,\n[t0]={t0} µM, [l0]={l0} µM", fontsize=12)
ax.plot(x_axis,y[0],'r-', label='Red chamber')
ax.plot(x_axis,y[1],'k-', label='White chamber')
ax.legend()
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)

ax.set_xlim(KD_beginning,KD_end)
ax.grid()
ax.set_ylabel(r"[Compound] (µM)", fontsize=12)
ax.set_xlabel(r"Compound K$_\mathrm{D}$ (µM)", fontsize=12)

plt.tight_layout(rect=(0,0,1,1))
plt.show()
