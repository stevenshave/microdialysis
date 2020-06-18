from numpy import sqrt


def ud_red_l_conc(t0, l0, kdtl, redvol, whitevol):
    """Calculate the compound concentration in the red chamber in a fully
    equilibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        kdtl (float): Kd of target-ligand interaction
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber

    Returns:
        float: Ligand concentration in the red chamber
    """    
    return (2*l0*redvol**2 + 3*l0*redvol*whitevol + kdtl*redvol*whitevol + t0*redvol*whitevol + l0*whitevol**2 + kdtl*whitevol**2 -
            sqrt((-2*l0*redvol**2 - 3*l0*redvol*whitevol - kdtl*redvol*whitevol - t0*redvol*whitevol - l0*whitevol**2 - kdtl*whitevol**2)**2 -
                 4*(redvol**2 + redvol*whitevol) * (l0**2*redvol**2 + 2*l0**2*redvol*whitevol + l0*kdtl*redvol*whitevol + l0*t0*redvol*whitevol +
                                                    l0**2*whitevol**2 + l0*kdtl*whitevol**2 + l0*t0*whitevol**2))) / (2.*(redvol**2 + redvol*whitevol))


def ud_red_l_conc_nonEq(t0, l0, kdtl, redvol, whitevol, pcvalue):
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


def ud_white_l_conc(t0, l0, kdtl, redvol, whitevol):
    """Calculate the compound concentration in the white chamber in a fully
    equilibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        kdtl (float): Kd of target-ligand interaction
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber

    Returns:
        float: Ligand concentration in the white chamber
    """
    return (l0*redvol - kdtl*redvol - t0*redvol + l0*whitevol - kdtl*whitevol +
            sqrt((-(l0*redvol) + kdtl*redvol + t0*redvol - l0*whitevol + kdtl*whitevol)**2 -
                 4*(redvol + whitevol)*(-(l0*kdtl*redvol) - l0*kdtl*whitevol)))/(2.*(redvol + whitevol))


def ud_white_l_conc_nonEq(t0, l0, kdtl, redvol, whitevol, pcvalue):
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


def ud_Kd_from_ptvalue(t0, l0, redvol, whitevol, ptvalue):
    """Calculate the protein-ligand interaction Kd from Pt in a fully equilibrated system

    Args:
        t0 (float): Target concentration (in the red chamber)
        l0 (float): Ligand concentration, over the entire volume of red and white chambers when fully equilibrated.
        redvol (float): Volume of the red chamber
        whitevol (float): Volume of the white chamber
        ptvalue (float): Pt - Ligand partition coefficient in the presence of protein

    Returns:
        float: Kd of the target-ligand interaction
    """ 
    return (l0*redvol - l0*ptvalue*redvol + t0*ptvalue*redvol + l0*whitevol + t0*whitevol -
            l0*ptvalue*whitevol)/((-1 + ptvalue)*(ptvalue*redvol + whitevol))


def ud_Kd_from_ptvalue_nonEq(t0, l0, redvol, whitevol, ptvalue, pcvalue):
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
