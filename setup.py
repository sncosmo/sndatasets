#!/usr/bin/env python
import os
import re
from distutils.core import setup


# Need a recursive glob to find all package data files if there are
# subdirectories
import fnmatch
def recursive_glob(basedir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(basedir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches


# Synchronize version from code.
version = re.findall(r"__version__ = \"(.*?)\"",
                     open(os.path.join("sndatasets", "__init__.py")).read())[0]

setup(name="sndatasets", 
      version=version,
      description="Download and normalize published supernova photometric data",
      long_description="",
      classifiers = ["Development Status :: 3 - Alpha",
                     "Programming Language :: Python :: 2",
                     "Programming Language :: Python :: 3",
                     "License :: OSI Approved :: MIT License",
                     "Topic :: Scientific/Engineering",
                     "Topic :: Scientific/Engineering :: Astronomy",
                     "Intended Audience :: Science/Research"],
      packages=["sndatasets"],
      package_data={"sndatasets": ["data/*/*"]},
      url="http://github.com/kbarbary/sndatasets",
      author="Kyle Barbary, David Rubin",
      author_email="kylebarbary@gmail.com")
