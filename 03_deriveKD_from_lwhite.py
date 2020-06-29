"""
Derive KD from lred

Derive KD for target-ligand complex using microdialysis derived lwhite values-
compound concentration in the white chamber (without protein),
"""

from microdialysis_equations import *
# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.  If uM is used for concentrations, then ul should be
# used as the standard unit for volume etc.


#  - l0 is total compound, or ligand concentration over the entire volume
#     if no protein was present (25 equal to starting in the red chamber
#     with a conc of 100 µM if chambers are 100 and 300 ul)
#  - t0 is target, or protein concentration in the red chamber, kdpc is the Kd
#     of pc complex, redvol and whitevol are volumes of respective chambers.
#  - redvol is the volume of the red chamber (the protein-containing chamber)
#  - whitevol is the volume of the white chamber (the no protein chamber)
#  - pc is the partition coefficient of the red and white chambers in a control
#     experiment in the absence of protein.  A perfectly equlibrating compound
#     will have a pc value of 1.0
#  - lwhite is the experimentally determined concentration of compound in the
#     white chamber

t0=80
l0=50
redvol=100
whitevol=300
pc=1.0
lwhite=39.2

system_parameters={v:eval(v) for v in ["t0","l0","redvol","whitevol","pc","lwhite"]}

print("Finding KD for the system in",system_parameters)
print(f"Kd={qud_Kd_from_lwhite(**system_parameters):.4f}")
