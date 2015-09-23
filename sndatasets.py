"""
sndatasets : Download and normalize published SN photometric data.
"""

from __future__ import print_function

import math
import os
from os.path import join
from six.moves.urllib.request import urlopen
from collections import OrderedDict

import numpy as np
from astropy.io import ascii
from astropy.table import Table

# CDS FTP:
# CDS_PREFIX = "ftp://cdsarc.u-strasbg.fr/pub/cats/J/"
# example postfix: ApJ/686/749/table10.[dat,fit]
# but FITS download through FTP seems broken, so we use http here.
CDS_PREFIX = "http://cdsarc.u-strasbg.fr/vizier/ftp/cats/"
#CDS_PREFIX = "http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/"
CACHE_DIR = "cache"

def _download_file(url, dst):
    """Download a file from url, save to CACHE_DIR/dst"""

    dst = join(CACHE_DIR, dst)
    if os.path.exists(dst):
        return

    # Ensure target directory exists.
    os.makedirs(os.path.dirname(dst), exist_ok=True)

    r = urlopen(url)
    with open(dst, "wb") as f:
        f.write(r.read())

def hms_to_deg(h, m, s):
    return 15. * (h + m / 60. + s / 3600.)

def sdms_to_deg(sign, d, m, s):
    sign = 1. - 2. * (sign == '-')
    return sign * (d + m / 60. + s / 3600.)

def jd_to_mjd(t):
    return t - 2400000.5

def mag_to_flux(m, me, zp):
    """Convert magnitude and magnitude error to flux, given a zeropoint."""
    
    f = 10.**(0.4 * (zp - m))
    fe = math.log(10) * 0.4 * me * f

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
    # dv = CMB_DZ * np.dot(CMB_XYZ, coords_xyz) * 299792.458
    # TODO: is dv needed?

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


def fetch_nearby99():
    """
    Nearby 99 set from Kowalski et al 2008
    http://adsabs.harvard.edu/abs/2008ApJ...686..749K
    """
    _download_file(CDS_PREFIX + "J/ApJ/686/749/ReadMe",
                   join("k08", "ReadMe"))
    _download_file(CDS_PREFIX + "J/ApJ/686/749/table1.dat",
                   join("k08", "table1.dat"))
    _download_file(CDS_PREFIX + "J/ApJ/686/749/table10.dat",
                   join("k08","table10.dat"))

    # Parse SN coordinates and redshifts
    meta = ascii.read("cache/k08/table1.dat", format='cds',
                      readme="cache/k08/ReadMe")
    ra = hms_to_deg(meta['RAh'], meta['RAm'], meta['RAs'])
    dec = sdms_to_deg(meta['DE-'], meta['DEd'], meta['DEm'], meta['DEs'])

    data = ascii.read("cache/k08/table10.dat", format='cds',
                      readme="cache/k08/ReadMe")

    # pivot photometry table.
    data['Tel'] = np.char.replace(data['Tel'], ' ', '_')
    names = ('SN', 'JD', 'band', 'mag', 'e_mag')
    dtypes = ('U9', 'f8', 'U16', 'f8', 'f8')
    reshaped = Table([np.empty(4 * len(data), dtype=dt) for dt in dtypes],
                     names=names, copy=False)
    for i, band in enumerate(('B', 'V', 'R', 'I')):
        bandsuffix = '_' + band
        reshaped['SN'][i::4] = data['SN']
        reshaped['JD'][i::4] = data['JD']
        reshaped['band'][i::4] = np.char.add(data['Tel'], '_' + band)
        reshaped['mag'][i::4] = data[band + 'mag']
        reshaped['e_mag'][i::4] = data['e_' + band + 'mag']

    # Eliminate missing value rows
    reshaped = reshaped[reshaped['mag'] != 0.]

    # Split up table into one table per SN and add metadata.
    sne = OrderedDict()
    for i in range(len(meta)):
        name = meta['SN'][i]
        sndata = reshaped[reshaped['SN'] == name]
        snmeta = OrderedDict([('name', name),
                              ('dataset', 'nearby99'),
                              ('z_helio', meta['z'][i]),
                              ('ra', ra[i]),
                              ('dec', dec[i])])
        zp = 29. * np.ones_like(sndata['mag'])
        zpsys = len(sndata) * ['vega']
        flux, fluxerr = mag_to_flux(sndata['mag'], sndata['e_mag'], zp)
        sne[name] = Table([jd_to_mjd(sndata['JD']), sndata['band'],
                           flux, fluxerr, zp, zpsys],
                          names=('time', 'band', 'flux', 'fluxerr', 'zp',
                                 'zpsys'),
                          meta=snmeta)
        # TODO: correct descriptions on table columns.

    return sne


# testing
if __name__ == "__main__":
    sne = fetch_nearby99()

