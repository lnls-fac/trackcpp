#!/usr/bin/env python-sirius

from setuptools import setup
import pathlib


def get_abs_path(relative):
    """."""
    return str(pathlib.Path(__file__).parent / relative)


with open(get_abs_path("../python_package/README.md"), "r") as _f:
    _long_description = _f.read().strip()

# with open(get_abs_path("VERSION"), "r") as _f:
#     __version__ = _f.read().strip()
__version__ = '0.0.1'

with open(get_abs_path("../python_package/requirements.txt"), "r") as _f:
    _requirements = _f.read().strip().split("\n")

setup(
    name='trackpybind',
    version=__version__,
    author='lnls-fac',
    description='trackcpp python package',
    long_description=_long_description,
    url='https://github.com/lnls-fac/trackcpp',
    download_url='https://github.com/lnls-fac/trackcpp',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    install_requires=_requirements,
    packages=['trackpybind'],
    package_data={'trackpybind': ['trackcpp_pybind.so', 'VERSION']},
    zip_safe=False
)
