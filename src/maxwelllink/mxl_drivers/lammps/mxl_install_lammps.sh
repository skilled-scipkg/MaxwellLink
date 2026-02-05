#!/bin/bash

#--------------------------------------------------------------------------------------#
# Copyright (c) 2026 MaxwellLink                                                       #
# This file is part of MaxwellLink. Repository: https://github.com/TaoELi/MaxwellLink  #
# If you use this code, always credit and cite arXiv:2512.06173.                       #
# See AGENTS.md and README.md for details.                                             #
#--------------------------------------------------------------------------------------#

if [ ! -d build ]; then
    mkdir build
fi

cd build/

echo "Installing modified LAMMPS with MaxwellLink..."
wget -c https://github.com/lammps/lammps/releases/download/stable_29Aug2024_update1/lammps-src-29Aug2024_update1.tar.gz
tar -xvf lammps-src-29Aug2024_update1.tar.gz
cd lammps-29Aug2024
mkdir build 
cd build

# copy modified files
cp ../../../src/maxwelllink/mxl_drivers/lammps/fix_maxwelllink.* ../src/MISC

# build LAMMPS with no GPU and minimal packages and install
cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake -D PKG_GPU=off ../cmake
make -j4

# copy the lmp executable to a known location
location=$(which mxl_driver.py)
# at the end of location is /mxl_driver.py, remove it to get the directory
dir=${location%/*}
echo "Copying lmp executable to ${dir}/lmp_mxl"
cp lmp ${dir}/lmp_mxl
