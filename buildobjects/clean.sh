#!bin/bash

# Cleaning first
cd objects
rm -rf *.obj *.lib *.exe
cd ../pythonbuild
rm -rf *.obj *.cxx *.lib *.exp __pycache__
cd ../trackcpp
rm "trackcpp.py" "_trackcpp.pyd"
cd ..
