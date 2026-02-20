#!/bin/bash
# this job should be run after submit_main.sh and after the tcp host and port info file is created in file

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=128        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J driver_ase      # Job name
#SBATCH -o driver_ase.o%j      # Name of stdout output file
#SBATCH -e driver_ase.e%j      # Name of stderr error file
#SBATCH -p shared # Queue (partition) name

# 0. wait for a few seconds to ensure the main job is running and the tcp host and port info file is created
sleep 10s

# 1. read tcp host and port info from file
HOST_PORT_FILE="tcp_host_port_info.txt"
if [[ ! -f "$HOST_PORT_FILE" ]]; then
    echo "Error: Host and port info file '$HOST_PORT_FILE' not found!"
    exit 1
fi
HOST=$(sed -n '1p' "$HOST_PORT_FILE")
PORT=$(sed -n '2p' "$HOST_PORT_FILE")

# 2. launch drivers
# nmol should equal to slurm ntasks
nmol=128
if [ "$nmol" -gt 1 ]; then
    echo "Launching mxl_driver for multiple molecules (nmol=$nmol)..."
    # run nmol-1 drivers in the background
    for ((i=1; i<nmol; i++)); do
        echo "Starting background mxl_driver instance $i..."
        mxl_driver --model ase --address $HOST --port $PORT --param  "atoms=../HCN_benchmark/hcn.xyz, calculator=psi4, calc_kwargs=method=b3lyp, basis=cc-pvdz, memory=2GB, num_threads=1, charges=[-0.2467948  -0.00296554  0.24976034]" &
    done
    # run the last driver in the foreground
    echo "Starting foreground mxl_driver instance $nmol..."
    mxl_driver --model ase --address $HOST --port $PORT --param  "atoms=../HCN_benchmark/hcn.xyz, calculator=psi4, calc_kwargs=method=b3lyp, basis=cc-pvdz, memory=2GB, num_threads=1, charges=[-0.2467948  -0.00296554  0.24976034]" 
else
    echo "Launching mxl_driver for single molecule (nmol=$nmol)..."
    mxl_driver --model ase --address $HOST --port $PORT --param  "atoms=../HCN_benchmark/hcn.xyz, calculator=psi4, calc_kwargs=method=b3lyp, basis=cc-pvdz, memory=2GB, num_threads=1, charges=[-0.2467948  -0.00296554  0.24976034]" --verbose
fi
