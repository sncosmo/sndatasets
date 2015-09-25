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

def _download_file(url, subdir):
    """Download a file from url, save to CACHE_DIR/subdir"""

    fname = url.split('/')[-1]
    destdir = join(CACHE_DIR, subdir)
    dest = join(destdir, fname)

    if os.path.exists(dest):
        return

    os.makedirs(destdir, exist_ok=True)

    r = urlopen(url)
    with open(dest, "wb") as f:
        f.write(r.read())


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


def fetch_kowalski08():
    """
    Nearby 99 set from Kowalski et al 2008
    http://adsabs.harvard.edu/abs/2008ApJ...686..749K
    """
    _download_file(CDS_PREFIX + "J/ApJ/686/749/ReadMe", "k08")
    _download_file(CDS_PREFIX + "J/ApJ/686/749/table1.dat", "k08")
    _download_file(CDS_PREFIX + "J/ApJ/686/749/table10.dat", "k08")

    # Parse SN coordinates and redshifts
    meta = ascii.read("cache/k08/table1.dat", format='cds',
                      readme="cache/k08/ReadMe")
    ra = hms_to_deg(meta['RAh'], meta['RAm'], meta['RAs'])
    dec = sdms_to_deg(meta['DE-'], meta['DEd'], meta['DEm'], meta['DEs'])

    data = ascii.read("cache/k08/table10.dat", format='cds',
                      readme="cache/k08/ReadMe")
    data = data.filled(0.)  # convert from masked table

    data = pivot_table(data, 'band', ['{}mag', 'e_{}mag'],
                       ['B', 'V', 'R', 'I'])
    data = data[data['mag'] != 0.]  # eliminate missing values

    # Join telescope and band into one column
    data['band'] = np.char.add(np.char.replace(data['Tel'], ' ', '_'),
                               np.char.add('_', data['band']))
    del data['Tel']

    # Split up table into one table per SN and add metadata.
    sne = OrderedDict()
    for i in range(len(meta)):
        name = meta['SN'][i]
        sndata = data[data['SN'] == name]
        snmeta = OrderedDict([('name', name),
                              ('dataset', 'kowalski08'),
                              ('z_helio', meta['z'][i]),
                              ('ra', ra[i]),
                              ('dec', dec[i])])
        zp = 29. * np.ones(len(sndata), dtype=np.float64)
        zpsys = len(sndata) * ['vega']
        flux, fluxerr = mag_to_flux(sndata['mag'], sndata['e_mag'], zp)
        sne[name] = Table([jd_to_mjd(sndata['JD']), sndata['band'],
                           flux, fluxerr, zp, zpsys],
                          names=('time', 'band', 'flux', 'fluxerr', 'zp',
                                 'zpsys'),
                          meta=snmeta)
        # TODO: correct descriptions on table columns.

    return sne

def fetch_hamuy96():
    """Hamuy et al. 1996 AJ 112 2408 "Calan Tololo" sample
    http://adsabs.harvard.edu/abs/1996AJ....112.2408H

    Photometry has been corrected to Bessell filters.
    """

    _download_file(CDS_PREFIX + "J/AJ/112/2408/ReadMe", "h96")
    _download_file(CDS_PREFIX + "J/AJ/112/2408/table4.dat", "h96")

    # TODO authoritative source for this metadata?
    # NOTE: commented-out lines are SNe not in the phtometric data table.
    #                  ra             dec            z_helio
    meta = {'1990O':  ('17:15:35.92', '+16:19:25.8', 0.0303),
            '1990T':  ('19:59:02.28', '-56:15:30.0', 0.0404),
            '1990Y':  ('03:37:22.64', '-33:02:40.1', 0.0391),
            '1990af': ('21:34:58.12', '-62:44:07.4', 0.0506),
            '1991S':  ('10:29:27.79', '+22:00:46.4', 0.0546),
            '1991U':  ('13:23:22.20', '-26:06:28.7', 0.0317),
            '1991ag': ('20:00:08.65', '-55:22:03.4', 0.0141),
            '1992J':  ('10:09:00.30', '-26:38:24.4', 0.0446),
            '1992K':  ('13:10:04.20', '-46:26:30.3', 0.0103),
            #'1992O':  ('19:23:42.29', '-62:49:30.1', 0.037),
            '1992P':  ('12:42:48.95', '+10:21:37.5', 0.0252),
            '1992ae': ('21:28:17.66', '-61:33:00.0', 0.0752),
            '1992ag': ('13:24:10.12', '-23:52:39.3', 0.0249),
            #'1992ai': ('01:29:08.04', '-32:16:30.0', -1.0),
            '1992al': ('20:45:56.45', '-51:23:40.0', 0.0146),
            '1992aq': ('23:04:34.76', '-37:20:42.1', 0.1018),
            '1992au': ('00:10:40.27', '-49:56:43.3', 0.0614),
            '1992bc': ('03:05:17.28', '-39:33:39.7', 0.0202),
            '1992bg': ('07:41:56.53', '-62:31:08.8', 0.0352),
            '1992bh': ('04:59:27.55', '-58:49:44.2', 0.0450),
            '1992bk': ('03:43:01.90', '-53:37:56.8', 0.0581),
            '1992bl': ('23:15:13.25', '-44:44:34.5', 0.0437),
            '1992bo': ('01:21:58.44', '-34:12:43.5', 0.0189),
            '1992bp': ('03:36:37.95', '-18:21:13.7', 0.0793),
            '1992br': ('01:45:44.83', '-56:05:57.9', 0.0882),
            '1992bs': ('03:29:27.20', '-37:16:18.9', 0.0637),
            '1993B':  ('10:34:51.38', '-34:26:30.0', 0.0696),
            '1993H':  ('13:52:50.34', '-30:42:23.3', 0.0239),
            '1993M':  ('19:13:01.53', '-64:17:28.3', 0.090),
            '1993O':  ('13:31:07.87', '-33:12:50.5', 0.0510),
            #'1993T':  ('23:10:54.09', '-44:58:48.6', 0.088),
            #'1993af': ('05:08:00.71', '-37:29:18.0', 0.0034),
            '1993ag': ('10:03:35.00', '-35:27:47.6', 0.0490),
            '1993ah': ('23:51:50.27', '-27:57:47.0', 0.0297)}

    data = ascii.read("cache/h96/table4.dat", format='cds',
                      readme="cache/h96/ReadMe")
    data = data.filled(0.)

    data = pivot_table(data, 'band', ['{}mag', 'e_{}mag'],
                       ['B', 'V', 'R', 'I'])
    data = data[data['mag'] != 0.]  # eliminate missing values

    # Split up table into one table per SN and add metadata.
    sne = OrderedDict()
    for name in meta:
        snmeta = OrderedDict([('name', name),
                              ('dataset', 'hamuy96'),
                              ('z_helio', meta[name][2]),
                              ('ra', sxhr_to_deg(meta[name][0])),
                              ('dec', sx_to_deg(meta[name][1]))])

        sndata = data[data['SN'] == name]

        time = jd_to_mjd(sndata['HJD'])
        band = np.char.add('bessell', sndata['band'])
        zp = 29. * np.ones(len(sndata), dtype=np.float64)
        zpsys = len(sndata) * ['vega']
        flux, fluxerr = mag_to_flux(sndata['mag'], sndata['e_mag'], zp)

        sne[name] = Table([time, band, flux, fluxerr, zp, zpsys],
                          names=('time', 'band', 'flux', 'fluxerr', 'zp',
                                 'zpsys'),
                          meta=snmeta)

    return sne

    

# testing
if __name__ == "__main__":
    sne = fetch_kowalski08()

    sne = fetch_hamuy96()
