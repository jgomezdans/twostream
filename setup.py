#!/usr/bin/env python
"""Setup script for building 2stream's python bindings"""

if __name__ == "__main__":
    from numpy.distutils.core import setup, Extension
    # Global variables for this extension:
    name         = "twostream"  # name of the generated python extension (.so)
    description  = "Python wrappers for the 2stream radiative transfer model"
    long_description = "Python wrappers for the 2stream radiative transfer model from B Pinty"
    author       = "J Gomez-Dans/NCEO & University College London"
    author_email = "j.gomez-dans@ucl.ac.uk"
    url = "http://github.com/jgomezdans/twostream"
    long_description = open ("README.rst", 'r' ).read()
    setup(name=name,
        description=description, \
        long_description=long_description, \
        author=author, \
        author_email = author_email, \
        url=url, version="1.0.0",
        ext_modules=[Extension(name='twostream', \
          sources=['twostream/twostream.f90'],
          f2py_options=['only:']+['twostream_solver']+[':'])])    
    