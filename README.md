# sndatasets

Download and normalize published supernova photometric data

*Work in progress* based on code by David Rubin.

## Current API

The following functions return a dictionary of `astropy.table.Table`s,
with each table holding the photometric data for a supernova. The
dictionary is keyed by the supernova name.

- `sndatasets.fetch_kowalski08()`: Nearby 99 sample from Kowalski et al 2008.
- `sndatasets.fetch_hamuy96()`: Calan-Tololo sample from Hamuy et al1996.

Uses `sndatasets.CACHE_DIR` to cache downloaded files. Default is `'cache'`.
