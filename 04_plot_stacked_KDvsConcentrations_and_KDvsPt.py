"""
Simulate microdialysis experiments, (KD vs compound concs & KD vs Pt,)

Produce a stacked plots showing KD vs compound concentrations in the red and
white chambers, and KD vs Pt in a microdialysis experiments
"""




from matplotlib import pyplot as plt
import numpy as np
PROTEIN_CONC=80
from microdialysis_equations import *

L0_CONC=25
REDVOL=100
WHITEVOL=300
XAXIS_BEGINNING = 0  # pKD of 3 is µM
XAXIS_END = 200  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 1000 # Publication used 2000 pts along X
x_axis = np.linspace(XAXIS_BEGINNING,XAXIS_END, NUM_POINTS_ON_XAXIS)

y=np.full((2,NUM_POINTS_ON_XAXIS), np.nan)
y[0]=ud_red_l_conc(PROTEIN_CONC, L0_CONC, x_axis, REDVOL, WHITEVOL)
y[1]=ud_white_l_conc(PROTEIN_CONC, L0_CONC, x_axis, REDVOL, WHITEVOL)

fig, ax = plt.subplots(2,1, figsize=(8, 6), sharex=True)
fig.suptitle("qµD simulation, K$_\mathrm{D}$ vs red and white chamber compound concentration and p-values"+
f"\nRed chamber = {REDVOL} µl, White chamber = {WHITEVOL} µl, [Protein]={PROTEIN_CONC} µM, [Compound]={L0_CONC} µM")
ax[0].plot(x_axis,y[0],'k-', label='Red chamber')
ax[0].plot(x_axis,y[1],'k--', label='White chamber')
ax[1].plot(x_axis,y[0]/y[1], 'k', label='P$_t$')
ax[0].legend()
ax[0].set_xlim(XAXIS_BEGINNING,XAXIS_END)
ax[0].grid()
ax[1].set_ylim(0,10)
ax[1].hlines(1,XAXIS_BEGINNING,XAXIS_END, linestyle='dotted', label="Equlibrium (P$_t$ = 1)")
ax[1].legend()
ax[1].grid()
ax[0].set_ylabel(r"[Compound] (µM)")
ax[1].set_ylabel(r"P$_t$")
ax[1].set_xlabel(r"Compound K$_\mathrm{D}$ (µM)")
plt.tight_layout(h_pad=0.04, rect=(0,0,1,0.9))
plt.show()
