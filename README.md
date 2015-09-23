# sndatasets

Download and normalize published supernova photometric data

*Work in progress* based on code by David Rubin.

## Current API

`sndatasets.fetch_nearby99()`: download and parse Nearby set from Kowalski et al 2008. Returns dictionary of astropy Tables, keyed by SN name.

Uses `sndatasets.CACHE_DIR` to cache downloaded files. Default is `'cache'`.