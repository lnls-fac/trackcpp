#!/usr/bin/env python-sirius

from setuptools import setup

with open('trackcpp\\VERSION','r') as _f:
    __version__ = _f.read().strip()

setup(
    name='trackcpp',
    version=__version__,
    author='lnls-fac',
    description='trackcpp python package',
    url='https://github.com/lnls-fac/trackcpp',
    download_url='https://github.com/lnls-fac/trackcpp',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['trackcpp'],
    package_data={'trackcpp': ['_trackcpp.pyd', 'VERSION']},
    zip_safe=False
)
