#!/bin/bash

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=128       # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J plasmon_nomol      # Job name
#SBATCH -o main.o%j      # Name of stdout output file
#SBATCH -e main.e%j      # Name of stderr error file
#SBATCH -p shared  # Queue (partition) name

np=128;
nmol=1;

for a in 2.79; do
    for r_frac in 0.4; do
        r=$(printf "%0.2f" `echo "${a}*${r_frac}" |bc`);           
        mpirun -np ${np} python -u emitter.py -dielectric -aa ${a} -rr ${r} -nmol ${nmol} > flux_a${a}_r${r}.out;
        grep flux1: flux_a${a}_r${r}.out |cut -d , -f2- > flux_a${a}_r${r}.dat;
    done;
done;



