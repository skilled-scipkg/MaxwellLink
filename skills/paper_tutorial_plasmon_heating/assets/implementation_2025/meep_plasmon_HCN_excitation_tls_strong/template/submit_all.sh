#!/bin/bash

# simulation with no maxwelllink to generate spectruw
# submit_nomol.sh

# the driver code should be submitted after the main code is submitted and running

job_main_id=$(sbatch submit_main.sh | awk '{print $4}')
echo "Submitted main job with Job ID: ${job_main_id}"

nmol=256;
if [ ${nmol} -le 128 ]; then
    sbatch --dependency=after:${job_main_id} submit_driver.sh
else
    # submit multiple drivers jobs if nmol > 128
    ndrivers=$(( (nmol + 127) / 128 ))
    for ((j=0; j<ndrivers; j++)); do
        sbatch --dependency=after:${job_main_id} submit_driver.sh
    done
fi
