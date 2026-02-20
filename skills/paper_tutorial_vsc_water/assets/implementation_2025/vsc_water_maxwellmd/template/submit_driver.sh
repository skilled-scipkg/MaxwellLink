#!/bin/bash
# this job should be run after submit_main.sh and after the tcp host and port info file is created in file

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=1        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J vsc_water      # Job name
#SBATCH -o driver.o%j      # Name of stdout output file
#SBATCH -e driver.e%j      # Name of stderr error file
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

# 2. modify lammps template input file
cp in_mxl.lmp in.lmp
sed -i -e "s/HOST/$HOST/g" in.lmp
sed -i -e "s/PORT/$PORT/g" in.lmp

# check whether lmp_mxl command is available
if ! command -v lmp_mxl &> /dev/null
then
    echo "lmp_mxl could not be found. Please ensure LAMMPS with MaxwellLink is installed and lmp_mxl is in your PATH."
    exit 1
fi

# launch lammps from bash script
srun -n $SLURM_NTASKS  lmp_mxl -in in.lmp
