"""
Simulate microdialysis experiment behaviour (KD vs Pt)

Simulate KD vs Pt in a microdialysis experiment, using both the full direct
analytical solution, and the simplified Weidemann equation

"""


from matplotlib import pyplot as plt
import numpy as np

from microdialysis_equations import *

t0=80
l0=50
redvol=100
whitevol=300
XAXIS_BEGINNING = 0
XAXIS_END = 500
pc=1.0
NUM_POINTS_ON_XAXIS = 1000 # Publication used 1000 pts along X
x_axis = np.linspace(XAXIS_BEGINNING,XAXIS_END, NUM_POINTS_ON_XAXIS)

y=np.full((2,NUM_POINTS_ON_XAXIS), np.nan)
y[0]=qud_lred(t0, l0, x_axis, redvol, whitevol, pc)
y[1]=qud_lwhite(t0, l0, x_axis, redvol, whitevol, pc)

fig, ax = plt.subplots(1,1, figsize=(7.204724, 5.09424929292), sharex=True)
fig.suptitle("qµD simulation", y=0.95, fontsize=14)

# Uncomment the line bellow to display the system as simulated by the
# simplified Weidemann equation, as detailed in PMID: 20681515
#ax.plot((t0/x_axis)+1, 'k--', label='Simplified Weidemann equation')

ax.plot(x_axis,y[0]/y[1], 'k', label='System')

ax.set_xlim(XAXIS_BEGINNING,XAXIS_END)
ax.grid()
ax.set_ylim(0,np.max(y[0]/y[1]))
ax.hlines(1,XAXIS_BEGINNING,XAXIS_END, linestyle='dotted', label="Equlibrium ($p_t$ = 1)")
ax.legend()
ax.set_ylabel("$p_t$",  fontsize=12)
ax.set_xlabel(r"Compound K$_\mathrm{D}$ (µM)", fontsize=12)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)

plt.tight_layout(rect=(0,0,1,1))
plt.show()
