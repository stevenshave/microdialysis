"""
Simulate microdialysis KD determination with error levels

Simulate KD vs Pt in a microdialysis experiment, using the 5% error level

"""
from matplotlib import pyplot as plt
import numpy as np
t0=80
from microdialysis_equations import *

lwhite_concs_to_get_kds_from=[35,40,45]

l0=50
redvol=100
whitevol=300
XAXIS_BEGINNING = 0
XAXIS_END = 500
NUM_POINTS_ON_XAXIS = 1000 # Publication used 1000 pts along X
x_axis = np.linspace(XAXIS_BEGINNING,XAXIS_END, NUM_POINTS_ON_XAXIS)

y=np.full((NUM_POINTS_ON_XAXIS), np.nan)
y=qud_lwhite(t0, l0, x_axis, redvol, whitevol, 1.0)

fig, ax = plt.subplots(1,1, figsize=(7.204724, 5.09424929292), sharex=True)
fig.suptitle("qµD simulation with 5% error", y=0.95, fontsize=14)
ax.plot(x_axis,y,'k', label='White chamber')

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)

for conc_i, conc in enumerate(lwhite_concs_to_get_kds_from):
    ax.axhline(conc*0.975)
    ax.axhline(conc*1.025)
    lowKd=qud_Kd_from_lwhite(conc*0.975,t0,l0,redvol,whitevol,1.0)
    highKd=qud_Kd_from_lwhite(conc*1.025,t0,l0,redvol,whitevol,1.0)
    ax.axvline(lowKd)
    ax.axvline(highKd)
    print(conc, lowKd, highKd)

ax.set_xlim(XAXIS_BEGINNING,XAXIS_END)
ax.set_ylim(25,51)

ax.grid()
ax.set_ylabel(r"[$lwhite$] (µM)", fontsize=12)
ax.set_xlabel(r"Compound K$_\mathrm{D}$ (µM)", fontsize=12)

plt.tight_layout(rect=(0,0,1,1))
plt.show()
