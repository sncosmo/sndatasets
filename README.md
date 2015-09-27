# sndatasets

Download and normalize published supernova photometric data

*Work in progress* based on code by David Rubin.

Other than making it easy to load SN datasets in a consistent fashion,
this package aims to be maximally reproducible. Therefore, data are
fetched from online sources whenever possible. In cases where data are
not available online (e.g., it had to be copy-pasted from the text of
a published paper), the data are included alongside the source code,
sometimes with hand-editing to make the data machine-readable. In addition,
some small metadata are currently hard-coded in the source code.

## Current API

The following functions return a dictionary of `astropy.table.Table`s,
with each table holding the photometric data for a supernova. The
dictionary is keyed by the supernova name.

- `sndatasets.fetch_kowalski08()`: Nearby 99 sample from Kowalski et al 2008.
- `sndatasets.fetch_hamuy96()`: Calan-Tololo sample from Hamuy et al1996.

Uses `sndatasets.CACHE_DIR` to cache downloaded files. Default is `'cache'`.
