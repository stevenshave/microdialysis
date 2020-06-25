"""
Simulate a microdialysis experiment

Simple example simulafting a microdialysis experiment with single point
values, once with a perfectly equilibrating compound, and then again with
one which does not.
"""


from microdialysis_equations import *
# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.  If uM is used for concentrations, then ul should be
# used as the standard unit for volume etc.


# Conc over the entire volume if no protein present (25 equal to starting in
# the red chamber with a conc of 100 µM if chambers are 100 and 300 ul)
experimental_parameters = {
    # l0 is total compound, or ligand concentration over the entire volume
    #  if no protein was present (25 equal to starting in the red chamber
    #  with a conc of 100 µM if chambers are 100 and 300 ul)
    # t0 is target, or protein concentration in the red chamber, kdpc is the Kd
    #  of pc complex, redvol and whitevol are volumes of respective chambers.
    'l0': 50,
    't0': 80,  # protein concentration
    'kdtl': 500,
    'redvol': 100,
    'whitevol': 300,
    'pcvalue':1.0
}
print("Perfectly equilibrating compounds")
print("*********************************")
print(
    f"Red volume concentration of {qud_lred (**experimental_parameters):.4f}")
print(
    f"White volume concentration of {qud_lwhite(**experimental_parameters):.4f}")
print(
    f"Pt-value = {qud_ptvalue(**experimental_parameters):.4f}")
    
print()
print("Poorly equilibrating compounds")
print("******************************")
experimental_parameters['pcvalue']=1.2
print(f"Red volume concentration of {qud_lred(**experimental_parameters):.4f}")
print(f"White volume concentration of {qud_lwhite(**experimental_parameters):.4f}")
print(f"Pt-value = {qud_ptvalue(**experimental_parameters):.4f}")

print()
print("Pt-value to Kd example")
print("*******************")
experimental_observations={
    'l0':50,
    't0':80,
    'redvol': 100,
    'whitevol': 300,
    'pcvalue':1.08,
    'ptvalue':1.43,
}
print(f"{experimental_observations=}")
print(f"Kd={qud_Kd_from_ptvalue(**experimental_observations):.4f}")
