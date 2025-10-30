#!/bin/bash

cd lmp_input/

# This script prepares LAMMPS input files for MaxwellLink tutorials.
PORT=$1
echo "Preparing LAMMPS input files with port $PORT..."

# modify in_mxl.lmp to set the correct port number
cp in_mxl.lmp in.lmp
sed -i -e "s/PORT_NUMBER/$PORT/g" in.lmp

# check whether lmp_mxl command is available
if ! command -v lmp_mxl &> /dev/null
then
    echo "lmp_mxl could not be found. Please ensure LAMMPS with MaxwellLink is installed and lmp_mxl is in your PATH."
    exit 1
fi

# launch lammps from bash script
lmp_mxl < in.lmp

cd ..