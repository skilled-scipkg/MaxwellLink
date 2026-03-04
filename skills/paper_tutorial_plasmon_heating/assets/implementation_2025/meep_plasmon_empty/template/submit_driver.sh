#!/bin/bash
# this job should be run after submit_main.sh and after the tcp host and port info file is created in file

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=16        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J driver_tls      # Job name
#SBATCH -o driver_tls.o%j      # Name of stdout output file
#SBATCH -e driver_tls.e%j      # Name of stderr error file
#SBATCH -p shared  # Queue (partition) name

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
        mxl_driver --model tls --address $HOST --port $PORT --param "omega=0.0157904, mu12=0.05, orientation=1, pe_initial=0.0" &
    done
    # run the last driver in the foreground
    echo "Starting foreground mxl_driver instance $nmol..."
    mxl_driver --model tls --address $HOST --port $PORT --param "omega=0.0157904, mu12=0.05, orientation=1, pe_initial=0.0" 
else
    echo "Launching mxl_driver for single molecule (nmol=$nmol)..."
    mxl_driver --model tls --address $HOST --port $PORT --param "omega=0.0157904, mu12=0.05, orientation=1, pe_initial=0.0" --verbose
fi
