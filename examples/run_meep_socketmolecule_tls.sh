#!/bin/bash

python meep_socketmolecule_tlsrelaxation.py &

echo "After launching MEEP, sleeping 20 seconds to wait for the socket server to start..."
sleep 20

c