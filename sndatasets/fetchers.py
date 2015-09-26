"""
sndatasets : Download and normalize published SN photometric data.
"""

from __future__ import print_function

from collections import OrderedDict
import os

import numpy as np
from astropy.io import ascii
from astropy.table import Table

from .utils import (download_file, hms_to_deg, sdms_to_deg, pivot_table,
                    mag_to_flux, jd_to_mjd, sxhr_to_deg, sx_to_deg,
                    fetch_sn_positions, CACHE_DIR, CDS_PREFIX)

__all__ = ["fetch_kowalski08", "fetch_hamuy96", "fetch_krisciunas"]


def fetch_kowalski08():
    """
    Nearby 99 set from Kowalski et al 2008
    http://adsabs.harvard.edu/abs/2008ApJ...686..749K
    """
    download_file(CDS_PREFIX + "J/ApJ/686/749/ReadMe", "kowalski08")
    download_file(CDS_PREFIX + "J/ApJ/686/749/table1.dat", "kowalski08")
    download_file(CDS_PREFIX + "J/ApJ/686/749/table10.dat", "kowalski08")

    # Parse SN coordinates and redshifts
    meta = ascii.read("cache/kowalski08/table1.dat", format='cds',
                      readme="cache/kowalski08/ReadMe")
    ra = hms_to_deg(meta['RAh'], meta['RAm'], meta['RAs'])
    dec = sdms_to_deg(meta['DE-'], meta['DEd'], meta['DEm'], meta['DEs'])

    data = ascii.read("cache/kowalski08/table10.dat", format='cds',
                      readme="cache/kowalski08/ReadMe")
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

    Position and heliocentric redshift metadata is hard-coded.
    """

    download_file(CDS_PREFIX + "J/AJ/112/2408/ReadMe", "hamuy96")
    download_file(CDS_PREFIX + "J/AJ/112/2408/table4.dat", "hamuy96")

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

    data = ascii.read("cache/hamuy96/table4.dat", format='cds',
                      readme="cache/hamuy96/ReadMe")
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
        band = np.char.add('bessell', np.char.lower(sndata['band']))
        zp = 29. * np.ones(len(sndata), dtype=np.float64)
        zpsys = len(sndata) * ['vega']
        flux, fluxerr = mag_to_flux(sndata['mag'], sndata['e_mag'], zp)

        sne[name] = Table([time, band, flux, fluxerr, zp, zpsys],
                          names=('time', 'band', 'flux', 'fluxerr', 'zp',
                                 'zpsys'),
                          meta=snmeta)

    return sne


def fetch_krisciunas():
    """Load the following SNe:

    1999aa 2000ApJ...539..658K Table 2
    1999cl 2000ApJ...539..658K Table 4
    1999cp 2000ApJ...539..658K Table 4

    Photometry has been corrected to Bessell filters.
    """

    from pkg_resources import resource_stream

    # Metadata
    z_helio = OrderedDict([("1999aa", 0.014443),
                           ("1999cc", 0.031328),
                           ("1999cl", 0.007609),
                           ("1999cp", 0.009480),
                           ("1999da", 0.012695),
                           ("1999dk", 0.014960),
                           ("1999ek", 0.017522),
                           ("1999gp", 0.026745),
                           ("2000bh", 0.022809),
                           ("2000bk", 0.025444),
                           ("2000ca", 0.023616),
                           ("2000ce", 0.016305),
                           ("2000cf", 0.036425),
                           ("2001ba", 0.029557),
                           ("2001bt", 0.014637),
                           ("2001cn", 0.015154),
                           ("2001cz", 0.015489),
                           ("2001el", 0.003896),
                           ("2002bo", 0.004240)])

    # Fetch positions to local file in cache.
    posfile = os.path.join(CACHE_DIR, "krisciunas/positions.csv")
    if not os.path.exists(posfile):
        os.makedirs(os.path.dirname(posfile), exist_ok=True)
    fetch_sn_positions(list(z_helio.keys()), posfile)

    f = resource_stream(__name__, 'data/krisciunas00/table2.txt')
