#!/bin/bash

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=128       # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J plasmon_nomol      # Job name
#SBATCH -o main.o%j      # Name of stdout output file
#SBATCH -e main.e%j      # Name of stderr error file
#SBATCH -p shared  # Queue (partition) name

np=128;

for a in 2.79; do
    mpirun -np ${np} python -u emitter.py -empty -aa ${a} > flux0_a${a}.out;
    grep flux1: flux0_a${a}.out |cut -d , -f2- > flux0_a${a}.dat;
done;



