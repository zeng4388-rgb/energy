#!/bin/bash
# run_chime_cat2.sh - Run E_iso function analysis on CHIME Catalog 2 data
# Usage: bash run_chime_cat2.sh

export OMP_NUM_THREADS=4

echo "=== CHIME/FRB Catalog 2 E_iso Function Analysis ==="
echo "Using CHIME Catalog 2 (5045 FRBs)"

for fgtype in ALG_NE2001 ALG_YMW16
do
    echo "Running with galaxy type: ${fgtype}"
    # All bursts
    echo "mpiexec -n 4 python3 nest_samp.py -fc chime_cat2.json -g ${fgtype} -upper 1 -o cat2_${fgtype}_upper"
    mpiexec -n 4 python3 nest_samp.py -fc chime_cat2.json -g ${fgtype} -upper 1 -o cat2_${fgtype}_upper

    # Only one-off bursts
    echo "mpiexec -n 4 python3 nest_samp.py -fc chime_cat2.json -g ${fgtype} -upper 1 -oneoff -o cat2_oneoff_${fgtype}_upper"
    mpiexec -n 4 python3 nest_samp.py -fc chime_cat2.json -g ${fgtype} -upper 1 -oneoff -o cat2_oneoff_${fgtype}_upper
done

echo "=== Done ==="
