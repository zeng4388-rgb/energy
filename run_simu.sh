#!/bin/bash
# run_simu.sh - Run simulation and recovery of E_iso function
# Usage: bash run_simu.sh

export OMP_NUM_THREADS=4

echo "=== FRB E_iso Function Simulation & Recovery ==="

# Create simulation directory
mkdir -p simu

# Generate mock data for two CHIME-like surveys
echo "Generating mock FRB data..."
for phis in 1e3 1e4
do
    echo "Simulating phis=${phis}, survey 1 (wide FoV)..."
    python3 simufrb.py -ns 100 -phis ${phis} -alpha -1.5 -logeisomax 41.0 -logeiso0 38.0 \
        -ga 1.4 -bw 400 -ts 25 -fov 200 -sn0 10 -out simu/simdat_${phis}_200.txt

    echo "Simulating phis=${phis}, survey 2 (narrow FoV)..."
    python3 simufrb.py -ns 100 -phis ${phis} -alpha -1.5 -logeisomax 41.0 -logeiso0 38.0 \
        -ga 0.7 -bw 400 -ts 50 -fov 20 -sn0 10 -out simu/simdat_${phis}_20.txt
done

# Run recovery
echo "Running E_iso function recovery..."
galaxy_type=ALG_YMW16
for phis in 1e3 1e4
do
    fout1=simdat_${phis}
    fout2=simdat_${phis}_upper
    fin1=./simu/simdat_${phis}_200.txt
    fin2=./simu/simdat_${phis}_20.txt

    echo "mpiexec -n 4 python3 nest_simu.py -f1 $fin1 -f2 $fin2 -o $fout1 -g ${galaxy_type}"
    mpiexec -n 4 python3 nest_simu.py -f1 $fin1 -f2 $fin2 -o $fout1 -g ${galaxy_type}

    echo "mpiexec -n 4 python3 nest_simu.py -upper 1 -f1 $fin1 -f2 $fin2 -o $fout2 -g ${galaxy_type}"
    mpiexec -n 4 python3 nest_simu.py -upper 1 -f1 $fin1 -f2 $fin2 -o $fout2 -g ${galaxy_type}
done

echo "=== Done ==="
