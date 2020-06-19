#from numpy import sqrt
from autograd.numpy import sqrt

def qud_lred(t0:float, l0:float, kdtl:float, redvol:float, whitevol:float, pcvalue:float):
    """Calculate the compound concentration in the red chamber in a partially equlibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        kdtl (float): Kd of target-ligand interaction
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pcvalue (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Ligand concentration in the red chamber
    """
    return (2*l0*pcvalue**2*redvol**2 + l0*pcvalue*redvol*whitevol +
            kdtl*pcvalue*redvol*whitevol + 2*l0*pcvalue**2*redvol*whitevol + pcvalue*t0*redvol*whitevol +
            kdtl*whitevol**2 + l0*pcvalue*whitevol**2 -
            sqrt((-2*l0*pcvalue**2*redvol**2 - l0*pcvalue*redvol*whitevol - kdtl*pcvalue*redvol*whitevol -
                  2*l0*pcvalue**2*redvol*whitevol - pcvalue*t0*redvol*whitevol - kdtl*whitevol**2 -
                  l0*pcvalue*whitevol**2)**2 -
                 4*(pcvalue**2*redvol**2 + pcvalue*redvol*whitevol) *
                 (l0**2*pcvalue**2*redvol**2 + l0*kdtl*pcvalue*redvol*whitevol +
                  2*l0**2*pcvalue**2*redvol*whitevol + l0*pcvalue*t0*redvol*whitevol +
                  l0*kdtl*pcvalue*whitevol**2 + l0**2*pcvalue**2*whitevol**2 + l0*pcvalue*t0*whitevol**2)))/(2.*(pcvalue**2*redvol**2 + pcvalue*redvol*whitevol))


def qud_lwhite(t0:float, l0:float, kdtl:float, redvol:float, whitevol:float, pcvalue:float):
    """Calculate the compound concentration in the white chamber in a partially equlibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        kdtl (float): Kd of target-ligand interaction
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pcvalue (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Ligand concentration in the white chamber
    """
    return (l0*pcvalue*redvol - kdtl*pcvalue*redvol - pcvalue*t0*redvol -
            kdtl*whitevol + l0*pcvalue*whitevol +
            sqrt(-4*(-(l0*kdtl*redvol) - l0*kdtl*whitevol)*(pcvalue**2*redvol + pcvalue*whitevol) +
                 (-(l0*pcvalue*redvol) + kdtl*pcvalue*redvol + pcvalue*t0*redvol + kdtl*whitevol -
                  l0*pcvalue*whitevol)**2))/(2.*(pcvalue**2*redvol + pcvalue*whitevol))


def qud_ptvalue(t0:float, l0:float, kdtl:float, redvol:float, whitevol:float, pcvalue:float):
    return (kdtl*pcvalue*redvol - l0*pcvalue*redvol + pcvalue*redvol*t0 - kdtl*whitevol - l0*pcvalue*whitevol + sqrt((-(kdtl*pcvalue*redvol) + l0*pcvalue*redvol - pcvalue*redvol*t0 + kdtl*whitevol + l0*pcvalue*whitevol)**2 - 4*kdtl*redvol*(-(l0*pcvalue**2*redvol) - kdtl*pcvalue*whitevol - l0*pcvalue**2*whitevol - pcvalue*t0*whitevol)))/(2.*kdtl*redvol)

def qud_Kd_from_ptvalue(ptvalue:float, t0:float, l0:float, redvol:float, whitevol:float, pcvalue:float):
    """Calculate the protein-ligand interaction Kd from Pt in a partially equilibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        ptvalue (float): Pt - Ligand partition coefficient in the presence of protein
        pcvalue (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Kd of the target-ligand interaction
    """
    return (-(l0*pcvalue**2*redvol) + l0*pcvalue*ptvalue*redvol - pcvalue*t0*ptvalue*redvol -
            l0*pcvalue**2*whitevol - pcvalue*t0*whitevol + l0*pcvalue*ptvalue*whitevol)/((pcvalue - ptvalue)*(ptvalue*redvol + whitevol))


def qud_Kd_from_lred(lred:float, t0:float, l0:float, redvol:float, whitevol:float, pcvalue:float):
    """Calculate the protein-ligand interaction Kd from ligand in red chamber in a partially equilibrated system

    Args:
        lred (float): Ligand concentration in the red chamber
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pcvalue (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Kd of the target-ligand interaction
    """
    return (-(l0**2*pcvalue**2*redvol**2) + 2*l0*lred*pcvalue**2*redvol**2 - lred**2*pcvalue**2*redvol**2 + l0*lred*pcvalue*redvol*whitevol - lred**2*pcvalue*redvol*whitevol -
            -       2*l0**2*pcvalue**2*redvol*whitevol + 2*l0*lred*pcvalue**2*redvol*whitevol - l0*pcvalue*redvol*t0*whitevol + lred*pcvalue*redvol*t0*whitevol + l0*lred*pcvalue*whitevol**2 - l0**2*pcvalue**2*whitevol**2 -
            -       l0*pcvalue*t0*whitevol**2)/(whitevol*(l0*pcvalue*redvol - lred*pcvalue*redvol - lred*whitevol + l0*pcvalue*whitevol))


def qud_Kd_from_lwhite(lwhite:float, t0:float, l0:float, redvol:float, whitevol:float, pcvalue:float):
    """Calculate the protein-ligand interaction Kd from ligand in white chamber in a partially equilibrated system

    Args:
        lwhite (float): Ligand concentration in the white chamber
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        pcvalue (float): Pc - Ligand partition coefficient in the absence of protein (control)

    Returns:
        float: Kd of the target-ligand interaction
    """
    return -((lwhite*(l0*pcvalue*redvol - lwhite*pcvalue**2*redvol - pcvalue*redvol*t0 + l0*pcvalue*whitevol - lwhite*pcvalue*whitevol))/(l0*redvol - lwhite*pcvalue*redvol + l0*whitevol - lwhite*whitevol))
