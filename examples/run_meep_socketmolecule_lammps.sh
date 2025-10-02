#!/bin/bash

python meep_socketmolecule_lammps.py &

echo "After launching MEEP, sleeping 10 seconds to wait for the socket server to start..."
sleep 10

lmp_mxl < lammps_input/in.lmp