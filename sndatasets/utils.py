from __future__ import print_function

import math
import os
from os.path import join
import sys
from six.moves.urllib.request import urlopen
from collections import OrderedDict

import numpy as np
from astropy.io import ascii
from astropy.table import Table


# CDS_PREFIX = "ftp://cdsarc.u-strasbg.fr/pub/cats/J/"
# example postfix: ApJ/686/749/table10.[dat,fit]
# but FITS download through FTP seems broken, so we use http here.
CDS_PREFIX = "http://cdsarc.u-strasbg.fr/vizier/ftp/cats/"
#CDS_PREFIX = "http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/"
CACHE_DIR = "cache"

def info(s):
    """Print a message with some formatting."""

    print("\033[1m\033[34mINFO: {}\033[0m".format(s))
    sys.stdout.flush()


def download_file(url, subdir):
    """Download a file from url, save to CACHE_DIR/subdir"""

    fname = url.split('/')[-1]
    destdir = join(CACHE_DIR, subdir)
    dest = join(destdir, fname)

    if os.path.exists(dest):
        return

    os.makedirs(destdir, exist_ok=True)

    info("get " + url)
    r = urlopen(url)
    with open(dest, "wb") as f:
        f.write(r.read())


def query_ned_position(name):
    """Return RA, Dec (J2000) for named objects, queried from NED.

    Parameters
    ----------
    names : str
        list of SN names, e.g., "SN 1999aa" or "SN1999aa".
    
    Returns
    -------
    ra, dec : float
        RA, Dec position in degrees.
    """

    
    url = ("http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname={}"
           "&out_csys=Equatorial&out_equinox=J2000.0&of=ascii_bar"
           "&list_limit=5&img_stamp=NO").format(name)
    info("fetching {} position from ned.ipac.caltech.edu".format(name))
    r = urlopen(url)
    lastline = r.read().decode('utf-8').split('\n')[-2]
    items = lastline.split('|')
    if items[0] != '1':
        raise RuntimeError("more than one position returned for {}"
                           .format(name))

    return float(items[2]), float(items[3])


def fetch_sn_positions(names, dest):
    """Fetch SN positions from NED and save to a target csv file.

    Parameters
    ----------
    names : list
        Names such as '1999aa' (no "SN" prefix). "SN" will be prepended in
        the query.
    """

    rows = ["name,ra,dec"]
    for name in names:
        ra, dec = query_ned_position('SN' + name)
        rows.append("{},{},{}".format(name, ra, dec))

    # Postpone writing file until after all the queries have succeeded.
    with open(dest, 'w') as f:
        for row in rows:
            f.write(row)
            f.write('\n')


def pivot_table(t, valuecolname, colpatterns, values):
    """Pivot the table. Similar to gather in tidyr but allows groups of
    columns.

    Parameters
    ----------
    t : astropy.table.Table
    valuecolname : str
    colpatterns : iterable of str
        e.g., ['{}mag', 'e_{}mag']
    values : iterable of str
        e.g., ['B', 'V', 'R', 'I']
    """

    # all column names we will remove.
    remove = [p.format(v) for v in values for p in colpatterns]

    # replacement column names.
    add = [p.format('') for p in colpatterns]
    add_dtypes = [t[p.format(values[0])].dtype for p in colpatterns]

    # existing columns to keep.
    keep = list(filter(lambda name: name not in remove, t.colnames))

    # create new table
    colnames = keep + [valuecolname] + add
    nrows = len(values) * len(t)
    data = [np.repeat(t[name], len(values)) for name in keep]
    data.append(np.tile(values, len(t)))
    data.extend([np.empty(nrows, dtype=dt) for dt in add_dtypes])
    reshaped = Table(data, names=colnames, copy=False, masked=t.masked)

    # fill new empty columns
    for i, v in enumerate(values):
        for p in colpatterns:
            addcol = p.format('')
            removecol = p.format(v)
            reshaped[addcol][i::len(values)] = t[removecol]

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
    x = math.cos(np.deg2rad(ra)) * math.cos(np.deg2rad(ra))
    y = math.cos(np.deg2rad(dec)) * math.sin(np.deg2rad(ra))
    z = math.sin(np.deg2rad(dec))

    return np.array([x, y, z], dtype=float64)


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
