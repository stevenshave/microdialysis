from matplotlib import pyplot as plt
import numpy as np
PROTEIN_CONC=40

from microdialysis_equations import *

L0_CONC=2
REDVOL=100
WHITEVOL=300
XAXIS_BEGINNING = 0  # pKD of 3 is µM
XAXIS_END = 500  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 1000 # Publication used 2000 pts along X
x_axis = np.linspace(XAXIS_BEGINNING,XAXIS_END, NUM_POINTS_ON_XAXIS)

y=np.full((2,NUM_POINTS_ON_XAXIS), np.nan)
y[0]=ud_red_c_conc(PROTEIN_CONC, L0_CONC, x_axis, REDVOL, WHITEVOL)
y[1]=ud_white_c_conc(PROTEIN_CONC, L0_CONC, x_axis, REDVOL, WHITEVOL)

fig, ax = plt.subplots(1,1, figsize=(8, 6), sharex=True)
fig.suptitle("qµD simulation, K$_\mathrm{D}$ vs P$_t$."+
f"\nRed chamber = {REDVOL} µl, White chamber = {WHITEVOL} µl,\n[Total protein]={PROTEIN_CONC} µM, [Total Compound]={L0_CONC} µM")
ax.plot((PROTEIN_CONC/x_axis)+1, 'k--', label='Simplified Weidemann equation')
ax.plot(x_axis,y[0]/y[1], 'k', label='System')

ax.set_xlim(XAXIS_BEGINNING,XAXIS_END)
ax.grid()
ax.set_ylim(0,10)
ax.hlines(1,XAXIS_BEGINNING,XAXIS_END, linestyle='dotted', label="Equlibrium (P$_t$ = 1)")
ax.legend()
ax.set_ylabel("P$_t$")
ax.set_xlabel(r"Compound K$_\mathrm{D}$ (µM)")
plt.tight_layout(h_pad=0.04, rect=(0,0,1,0.9))
plt.show()
