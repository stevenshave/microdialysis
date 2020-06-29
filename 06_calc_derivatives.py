"""
Differentiate redconc, whiteconc, pt

Simulate the rate of change of redconc, whiteconc and pts, also KD with
respect to these values.  Requires the autograd package. We redefine all the
equations present in microdialysis_equations.py so that sqrt is suppled by the
autograd package.
"""
import sys
from autograd import grad
from autograd.numpy import sqrt
t0 = 80.0
l0 = 50.0
redvol = 100.0
whitevol = 300.0
pc = 1.0
KDS_to_simulate = [1., 100., 200., 300., 400., 500.]

def qud_lred(t0: float, l0: float, kdtl: float, redvol: float, whitevol: float, pc: float):
    """Calculate the compound concentration in the red chamber in a partially equlibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        kdtl (float): Kd of target-ligand interaction
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pc (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Ligand concentration in the red chamber
    """
    return (2*l0*pc**2*redvol**2 + l0*pc*redvol*whitevol +
            kdtl*pc*redvol*whitevol + 2*l0*pc**2*redvol*whitevol + pc*t0*redvol*whitevol +
            kdtl*whitevol**2 + l0*pc*whitevol**2 -
            sqrt((-2*l0*pc**2*redvol**2 - l0*pc*redvol*whitevol - kdtl*pc*redvol*whitevol -
                  2*l0*pc**2*redvol*whitevol - pc*t0*redvol*whitevol - kdtl*whitevol**2 -
                  l0*pc*whitevol**2)**2 -
                 4*(pc**2*redvol**2 + pc*redvol*whitevol) *
                 (l0**2*pc**2*redvol**2 + l0*kdtl*pc*redvol*whitevol +
                  2*l0**2*pc**2*redvol*whitevol + l0*pc*t0*redvol*whitevol +
                  l0*kdtl*pc*whitevol**2 + l0**2*pc**2*whitevol**2 + l0*pc*t0*whitevol**2)))/(2.*(pc**2*redvol**2 + pc*redvol*whitevol))


def qud_lwhite(t0: float, l0: float, kdtl: float, redvol: float, whitevol: float, pc: float):
    """Calculate the compound concentration in the white chamber in a partially equlibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        kdtl (float): Kd of target-ligand interaction
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pc (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Ligand concentration in the white chamber
    """
    return (l0*pc*redvol - kdtl*pc*redvol - pc*t0*redvol -
            kdtl*whitevol + l0*pc*whitevol +
            sqrt(-4*(-(l0*kdtl*redvol) - l0*kdtl*whitevol)*(pc**2*redvol + pc*whitevol) +
                 (-(l0*pc*redvol) + kdtl*pc*redvol + pc*t0*redvol + kdtl*whitevol -
                  l0*pc*whitevol)**2))/(2.*(pc**2*redvol + pc*whitevol))


def qud_pt(t0: float, l0: float, kdtl: float, redvol: float, whitevol: float, pc: float):
    """Calculate the pt value in a partially equlibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        kdtl (float): Kd of target-ligand interaction
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pc (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: pt value
    """    
    return (kdtl*pc*redvol - l0*pc*redvol + pc*redvol*t0 - kdtl*whitevol - l0*pc*whitevol + sqrt((-(kdtl*pc*redvol) + l0*pc*redvol - pc*redvol*t0 + kdtl*whitevol + l0*pc*whitevol)**2 - 4*kdtl*redvol*(-(l0*pc**2*redvol) - kdtl*pc*whitevol - l0*pc**2*whitevol - pc*t0*whitevol)))/(2.*kdtl*redvol)


def qud_Kd_from_pt(pt: float, t0: float, l0: float, redvol: float, whitevol: float, pc: float):
    """Calculate the protein-ligand interaction Kd from Pt in a partially equilibrated system

    Args:
        pt (float): Pt value (lred/lwhite)
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pt (float): Pt - Ligand partition coefficient in the presence of protein
        pc (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Kd of the target-ligand interaction
    """
    return (-(l0*pc**2*redvol) + l0*pc*pt*redvol - pc*t0*pt*redvol -
            l0*pc**2*whitevol - pc*t0*whitevol + l0*pc*pt*whitevol)/((pc - pt)*(pt*redvol + whitevol))


def qud_Kd_from_lred(lred: float, t0: float, l0: float, redvol: float, whitevol: float, pc: float):
    """Calculate the protein-ligand interaction Kd from ligand in red chamber in a partially equilibrated system

    Args:
        lred (float): Ligand concentration in the red chamber
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pc (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Kd of the target-ligand interaction
    """
    return ((-l0 ** 2)*pc ** 2*redvol ** 2 + 2*l0*lred*pc ** 2*redvol ** 2 - lred ** 2*pc ** 2*redvol ** 2 +
            l0*lred*pc*redvol*whitevol - lred ** 2*pc*redvol*whitevol - 2*l0 ** 2*pc ** 2*redvol*whitevol +
            2*l0*lred*pc ** 2*redvol*whitevol - l0*pc*redvol*t0*whitevol + lred*pc*redvol*t0*whitevol + l0*lred*pc*whitevol ** 2 -
            l0 ** 2*pc ** 2*whitevol ** 2 - l0*pc*t0*whitevol ** 2)/(whitevol*(l0*pc*redvol - lred*pc*redvol -
                                                                                     lred*whitevol + l0*pc*whitevol))

def qud_Kd_from_lwhite(lwhite: float, t0: float, l0: float, redvol: float, whitevol: float, pc: float):
    """Calculate the protein-ligand interaction Kd from ligand in white chamber in a partially equilibrated system

    Args:
        lwhite (float): Ligand concentration in the white chamber
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pc (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Kd of the target-ligand interaction
    """
    return -((lwhite*(l0*pc*redvol - lwhite*pc**2*redvol - pc*redvol*t0 + l0*pc*whitevol - lwhite*pc*whitevol))/(l0*redvol - lwhite*pc*redvol + l0*whitevol - lwhite*whitevol))













g_lred_wrt_kd = grad(lambda kd: qud_lred(t0,l0,kd,redvol,whitevol,1.0))
g_lwhite_wrt_kd = grad(lambda kd: qud_lwhite(t0,l0,kd,redvol,whitevol,1.0))
g_pt_wrt_kd = grad(lambda kd: qud_pt(t0,l0,kd,redvol,whitevol,1.0))

print(f"{'KD':>10},{'lred':>10},{'dlreddKD':>10},{'lwhite':>10},{'dlwhitedKD':>10},{'pt':>10},{'dlptdKD':>10}")
for kd in KDS_to_simulate:
    lwhite = qud_lwhite(t0, l0, kd, redvol, whitevol, 1.0)
    lred = qud_lred(t0, l0, kd, redvol, whitevol, 1.0)
    pt = qud_pt(t0, l0, kd, redvol, whitevol, 1.0)
    print(f"{kd:>10.0f},{lred:>10.4f},{g_lred_wrt_kd(kd):>10.4f},{lwhite:>10.4f},{g_lwhite_wrt_kd(kd):>10.4f},{pt:>10.4f},{g_pt_wrt_kd(kd):>10.4f}")
print()


g_kd_wrt_lred = grad(lambda lred: qud_Kd_from_lred(lred,t0,l0,redvol,whitevol,1.0))
g_kd_wrt_lwhite = grad(lambda lwhite: qud_Kd_from_lwhite(lwhite,t0,l0,redvol,whitevol,1.0))
g_kd_wrt_ptval = grad(lambda pt: qud_Kd_from_pt(pt,t0,l0,redvol,whitevol,1.0))
print(f"{'KD':>10},{'lred':>10},{'dKDdlred':>10},{'lwhite':>10},{'dlKDdlwhite':>10},{'pt':>10},{'dKDdlpt':>10}")
for kd in KDS_to_simulate:
    lred = qud_lred(t0, l0, kd, redvol, whitevol, 1.0)
    lwhite = qud_lwhite(t0, l0, kd, redvol, whitevol, 1.0)
    pt = qud_pt(t0, l0, kd, redvol, whitevol, 1.0)
    print(f"{kd:>10.0f},{lred:>10.4f},{g_kd_wrt_lred(lred):>10.4f},{lwhite:>10.4f},{g_kd_wrt_lwhite(lwhite):>10.4f},{pt:>10.4f},{g_kd_wrt_ptval(pt):>10.4f}")
