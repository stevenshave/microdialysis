from numpy import sqrt


def ud_red_c_conc(p, c, kdpc, redvol, whitevol):
    return (2*c*redvol**2 + 3*c*redvol*whitevol + kdpc*redvol*whitevol + p*redvol*whitevol + c*whitevol**2 + kdpc*whitevol**2 -
            sqrt((-2*c*redvol**2 - 3*c*redvol*whitevol - kdpc*redvol*whitevol - p*redvol*whitevol - c*whitevol**2 - kdpc*whitevol**2)**2 -
                 4*(redvol**2 + redvol*whitevol) * (c**2*redvol**2 + 2*c**2*redvol*whitevol + c*kdpc*redvol*whitevol + c*p*redvol*whitevol +
                                                    c**2*whitevol**2 + c*kdpc*whitevol**2 + c*p*whitevol**2))) / (2.*(redvol**2 + redvol*whitevol))


def ud_red_c_conc_nonEq(p, c, kdpc, redvol, whitevol, no_protein_pval):
    return (2*c*no_protein_pval**2*redvol**2 + c*no_protein_pval*redvol*whitevol +
            kdpc*no_protein_pval*redvol*whitevol + 2*c*no_protein_pval**2*redvol*whitevol + no_protein_pval*p*redvol*whitevol +
            kdpc*whitevol**2 + c*no_protein_pval*whitevol**2 -
            sqrt((-2*c*no_protein_pval**2*redvol**2 - c*no_protein_pval*redvol*whitevol - kdpc*no_protein_pval*redvol*whitevol -
                  2*c*no_protein_pval**2*redvol*whitevol - no_protein_pval*p*redvol*whitevol - kdpc*whitevol**2 -
                  c*no_protein_pval*whitevol**2)**2 -
                 4*(no_protein_pval**2*redvol**2 + no_protein_pval*redvol*whitevol) *
                 (c**2*no_protein_pval**2*redvol**2 + c*kdpc*no_protein_pval*redvol*whitevol +
                  2*c**2*no_protein_pval**2*redvol*whitevol + c*no_protein_pval*p*redvol*whitevol +
                  c*kdpc*no_protein_pval*whitevol**2 + c**2*no_protein_pval**2*whitevol**2 + c*no_protein_pval*p*whitevol**2)))/(2.*(no_protein_pval**2*redvol**2 + no_protein_pval*redvol*whitevol))


def ud_white_c_conc(p, c, kdpc, redvol, whitevol):
    return (c*redvol - kdpc*redvol - p*redvol + c*whitevol - kdpc*whitevol +
            sqrt((-(c*redvol) + kdpc*redvol + p*redvol - c*whitevol + kdpc*whitevol)**2 -
                 4*(redvol + whitevol)*(-(c*kdpc*redvol) - c*kdpc*whitevol)))/(2.*(redvol + whitevol))


def ud_white_c_conc_nonEq(p, c, kdpc, redvol, whitevol, no_protein_pval):
    return (c*no_protein_pval*redvol - kdpc*no_protein_pval*redvol - no_protein_pval*p*redvol -
            kdpc*whitevol + c*no_protein_pval*whitevol +
            sqrt(-4*(-(c*kdpc*redvol) - c*kdpc*whitevol)*(no_protein_pval**2*redvol + no_protein_pval*whitevol) +
                 (-(c*no_protein_pval*redvol) + kdpc*no_protein_pval*redvol + no_protein_pval*p*redvol + kdpc*whitevol -
                  c*no_protein_pval*whitevol)**2))/(2.*(no_protein_pval**2*redvol + no_protein_pval*whitevol))


def ud_Kd_from_pval(p, c, redvol, whitevol, pval):
    return (c*redvol - c*pval*redvol + p*pval*redvol + c*whitevol + p*whitevol -
            c*pval*whitevol)/((-1 + pval)*(pval*redvol + whitevol))


def ud_Kd_from_pval_nonEq(p, c, redvol, whitevol, pval, no_protein_pval):
    return (-(c*no_protein_pval**2*redvol) + c*no_protein_pval*pval*redvol - no_protein_pval*p*pval*redvol -
            c*no_protein_pval**2*whitevol - no_protein_pval*p*whitevol + c*no_protein_pval*pval*whitevol)/((no_protein_pval - pval)*(pval*redvol + whitevol))
