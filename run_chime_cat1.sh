#!/bin/bash
# run_chime_cat1.sh - Run E_iso function analysis on CHIME Catalog 1 data
# Usage: bash run_chime_cat1.sh

export OMP_NUM_THREADS=4

echo "=== CHIME/FRB Catalog 1 E_iso Function Analysis ==="
echo "Using CHIME Catalog 1 (600 FRBs) with baseband-informed parameters"

# All bursts with different galaxy DM models
for fgtype in ALG_NE2001 ALG_YMW16 ETG LTG_NE2001 LTG_YMW16
do
    echo "Running with galaxy type: ${fgtype}"
    echo "mpiexec -n 4 python3 nest_samp.py -fc chime_cat1.json -g ${fgtype} -upper 1 -o cat1_${fgtype}_upper"
    mpiexec -n 4 python3 nest_samp.py -fc chime_cat1.json -g ${fgtype} -upper 1 -o cat1_${fgtype}_upper
done

echo "=== Done ==="
