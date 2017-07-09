"""General utilities."""

import math

import numpy as np
from astropy.table import Table


def pivot_table(t, value_col_name, colfmts, values, colfmts_replace=None,
                values_replace=None):
    """Pivot the table. Similar to gather in tidyr but allows groups of
    columns.

    Parameters
    ----------
    t : astropy.table.Table
    value_col_name : str
    colfmts : iterable of str
        e.g., ['{}mag', 'e_{}mag']
    values : iterable of str
        e.g., ['B', 'V', 'R', 'I']
    colfmts_replace : iterable of str, optional
        If given, new columns will have these names instead of `colfmts`
        with an empty string.
    values_replace : iterable of str, optional
        If given, values will be replaced with these 
    """

    # all column names we will remove.
    remove = [colfmt.format(v) for v in values for colfmt in colfmts]

    if values_replace is not None:
        assert len(values_replace) == len(values)
    else:
        values_replace = values

    # replacement column names.
    if colfmts_replace is not None:
        add = colfmts_replace
    else:
        add = [colfmt.format('') for colfmt in colfmts]

    add_dtypes = [t[colfmt.format(values[0])].dtype for colfmt in colfmts]

    # existing columns to keep.
    keep = list(filter(lambda name: name not in remove, t.colnames))

    # create new table
    colnames = keep + [value_col_name] + add
    nrows = len(values) * len(t)
    data = [np.repeat(t[name], len(values)) for name in keep]
    data.append(np.tile(values_replace, len(t)))
    data.extend([np.empty(nrows, dtype=dt) for dt in add_dtypes])
    reshaped = Table(data, names=colnames, copy=False, masked=t.masked)

    # fill new empty columns
    for i, v in enumerate(values):
        for colfmt, addname in zip(colfmts, add):
            removecol = colfmt.format(v)
            reshaped[addname][i::len(values)] = t[removecol]

    return reshaped


def hms_to_deg(h, m, s):
    return 15. * (h + m / 60. + s / 3600.)


def sxhr_to_deg(s):
    """sexagesimal hours to degrees"""
    h, m, s = s.split(':')
    return hms_to_deg(int(h), int(m), float(s))


def sdms_to_deg(sign, d, m, s):
    sign = 1. - 2. * (sign == '-')
    return sign * (d + m / 60. + s / 3600.)


def sx_to_deg(s):
    """sexagesimal to degrees. Sign must be the first character"""
    sign = s[0]
    d, m, s = s.split(':')
    return sdms_to_deg(sign, abs(int(d)), int(m), float(s))


def jd_to_mjd(t):
    return t - 2400000.5


def mag_to_flux(m, me, zp):
    """Convert magnitude and magnitude error to flux, given a zeropoint."""
    
    f = 10.**(0.4 * (zp - m))
    fe = math.log(10.) * 0.4 * me * f

    return f, fe


def radec_to_xyz(ra, dec):
    x = math.cos(np.deg2rad(dec)) * math.cos(np.deg2rad(ra))
    y = math.cos(np.deg2rad(dec)) * math.sin(np.deg2rad(ra))
    z = math.sin(np.deg2rad(dec))

    return np.array([x, y, z], dtype=np.float64)


def cmb_dz(ra, dec):
    """See http://arxiv.org/pdf/astro-ph/9609034
     CMBcoordsRA = 167.98750000 # J2000 Lineweaver
     CMBcoordsDEC = -7.22000000
    """

    # J2000 coords from NED
    CMB_DZ = 371000. / 299792458.
    CMB_RA = 168.01190437
    CMB_DEC = -6.98296811
    CMB_XYZ = radec_to_xyz(CMB_RA, CMB_DEC)

    coords_xyz = radec_to_xyz(ra, dec)
    
    dz = CMB_DZ * np.dot(CMB_XYZ, coords_xyz)

    return dz


def helio_to_cmb(z, ra, dec):
    """Convert from heliocentric redshift to CMB-frame redshift.
    
    Parameters
    ----------
    z : float
        Heliocentric redshift.
    ra, dec: float
        RA and Declination in degrees (J2000).
    """

    dz = -cmb_dz(ra, dec)
    one_plus_z_pec = math.sqrt((1. + dz) / (1. - dz))
    one_plus_z_CMB = (1. + z) / one_plus_z_pec

    return one_plus_z_CMB - 1.


def cmb_to_helio(z, ra, dec):
    """Convert from CMB-frame redshift to heliocentric redshift.
    
    Parameters
    ----------
    z : float
        CMB-frame redshift.
    ra, dec: float
        RA and Declination in degrees (J2000).
    """

    dz = -cmb_dz(ra, dec)
    one_plus_z_pec = math.sqrt((1. + dz) / (1. - dz))
    one_plus_z_helio = (1. + z) * one_plus_z_pec

    return one_plus_z_helio - 1.
