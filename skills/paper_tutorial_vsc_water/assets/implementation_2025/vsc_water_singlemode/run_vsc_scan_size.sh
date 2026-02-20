#!/bin/bash

for size in 1 2 4 8
do

    DIR="size_${size}_cube"
    mkdir -p ${DIR}
    cd ${DIR}

    # copy input files from template
    cp ../template/* .

    # modify input parameters to accomodate different system sizes
    # 1. modify lammps input file to enlarge the water molecular size
    cp ../lmp_input/in_mxl.lmp in_mxl.lmp
    cp ../lmp_input/data.lmp data.lmp

    sed -i -e "s/replicate 1 1 1/replicate ${size} ${size} ${size}/g" in_mxl.lmp
    # 2. modify maxwellmd python file to adjust the coupling strength accordingly
    sed -i -e "s/enlarge_factor = 1.0/enlarge_factor = ${size}/g" vsc_water_singlemode.py
    # 3. modify the driver submission script to change the number of MPI tasks accordingly
    ntasks=$(( size * size * size ))
    # 3.1 estimate the required nodes (each node has 128 cores)
    if (( ntasks % 128 == 0 )); then
        nodes=$(( ntasks / 128 ))
    else
        nodes=$(( ntasks / 128 + 1 ))
    fi
    # set the upper limit of nodes so that we do not abuse the HPC resource
    if (( nodes > 16 )); then
        nodes=16
        ntasks=$(( nodes * 128 ))
    fi
    echo "ntasks: ${ntasks}, nodes: ${nodes}"
    sed -i -e "s/--ntasks=1/--ntasks=${ntasks}/g" submit_driver.sh
    sed -i -e "s/--nodes=1/--nodes=${nodes}/g" submit_driver.sh
    # also switch from the shared queue to the wholenode queue for larger jobs
    if (( nodes > 1 )); then
        sed -i -e "s/#SBATCH -p shared/#SBATCH -p wholenode/g" submit_driver.sh
    fi

    sed -i -e "s/-J vsc_water/-J vscd_${size}/g" submit_driver.sh
    sed -i -e "s/-J vsc_water/-J vscm_${size}/g" submit_main.sh

    # 4. submit jobs
    ./submit_all.sh

    cd ..

done
