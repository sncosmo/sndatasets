"""Utilities for downloading and cacheing files."""

import os
from os.path import join
import sys
from six.moves.urllib.request import urlopen
from collections import OrderedDict

CACHE_DIR = "cache"


def info(s):
    """Print a message with some formatting."""

    print("\033[1m\033[34mINFO: {}\033[0m".format(s))
    sys.stdout.flush()


def download_file(url, cachedir)
    """Download a file from url, cache it, return path to cached file."""

    destdir = join(CACHE_DIR, subdir)
    destpath = join(destdir, url.split('/')[-1])

    if not os.path.exists(destpath):
        if not os.path.exists(destdir):
            os.makedirs(destdir)

        info("get " + url)
        r = urlopen(url)
        with open(destpath, "wb") as f:
            f.write(r.read())

    return destpath


def query_ned_position(name):
    """Return RA, Dec (J2000) for named objects, queried from NED.

    Parameters
    ----------
    names : str
        list of SN names, e.g., "1999aa" or "1999aa".
        "SN" is prefixed to query.
    
    Returns
    -------
    ra, dec : float
        RA, Dec position in degrees.
    """

    name = "SN" + name
    
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


def download_sn_positions(names, cachename):
    """Fetch SN positions from NED and save to a target csv file.

    Parameters
    ----------
    names : list
        Names such as '1999aa' (no "SN" prefix). "SN" will be prepended in
        the query.

    Returns
    -------
    destpath : str
        Path to downloaded file.
    """

    destpath = join(CACHE_DIR, cachename)

    if not os.path.exists(destpath):
        rows = ["name,ra,dec"]
        for name in names:
            ra, dec = query_ned_position(name)
            rows.append("{},{},{}".format(name, ra, dec))

        # Postpone writing file until after all the queries have succeeded.
        destdir = os.path.dirname(destpath)
        if not os.path.exists(destdir):
            os.makedirs(destdir)

        with open(destpath, 'w') as f:
            for row in rows:
                f.write(row)
                f.write('\n')

    return destpath
