"""
Simulate microdialysis experiment behaviour (KD vs Conc)

Simulate KD vs compound concentrations in the red and white chambers
of a microdialysis experiment
"""

from pathlib import Path
from matplotlib import pyplot as plt
plt.rcParams['animation.ffmpeg_path'] = Path("C:\\Users\\steve\\Downloads\\ffmpeg-20200716-d11cc74-win64-static\\bin\\ffmpeg.exe")
from matplotlib import animation
import numpy as np
t0s=list(range(1,101))+[100]*20
from microdialysis_equations import *

l0=50
redvol=100
whitevol=300
KD_beginning = 0
KD_end = 500
pc = 1.0
NUM_POINTS_ON_XAXIS = 1000 # Publication used 1000 pts along X
x_axis = np.linspace(KD_beginning,KD_end, NUM_POINTS_ON_XAXIS)


y=np.full((len(t0s),2,NUM_POINTS_ON_XAXIS), np.nan)
for t0_i,t0 in enumerate(t0s):
    y[t0_i,0]=qud_lred(t0, l0, x_axis, redvol, whitevol, pc)
    y[t0_i,1]=qud_lwhite(t0, l0, x_axis, redvol, whitevol, pc)

fig, ax = plt.subplots(1,1, figsize=(7.204724, 5.09424929292), sharex=True)
#fig, ax = plt.subplots(1,1, figsize=(8, 6), sharex=True)
fig.suptitle("qµD simulation",y=0.95, fontsize=14)#, K$_\mathrm{D}$ vs red and white chamber compound concentration."+
#f"\nRed chamber = {redvol} µl, White chamber = {whitevol} µl,\n[t0]={t0} µM, [l0]={l0} µM", fontsize=12)
red_line,=ax.plot(x_axis,y[0,0],'r-', label='Red chamber')
white_line,=ax.plot(x_axis,y[0,1],'k-', label='White chamber')

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=7200)


def update(num, x, y, line1, line2, ax, fig,t0s):
    line1.set_data(x, y[num,0])
    line2.set_data(x, y[num,1])
    fig.suptitle("T0 = "+str(t0s[num])+ " µM", fontsize=16,y=0.95)
    return [line1,line2,fig]

ani = animation.FuncAnimation(fig, update, len(t0s), fargs=[x_axis, y, red_line, white_line, ax, fig, t0s],interval=10, blit=True, repeat_delay=2000,)

ax.legend()
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)

ax.set_xlim(KD_beginning,KD_end)
ax.set_ylim(0,120)

ax.grid()
ax.set_ylabel(r"[Compound] (µM)", fontsize=12)
ax.set_xlabel(r"Compound K$_\mathrm{D}$ (µM)", fontsize=12)

plt.tight_layout(rect=(0,0,1,1))
ani.save("anim.mp4", writer=writer,dpi=150,)
plt.show()
