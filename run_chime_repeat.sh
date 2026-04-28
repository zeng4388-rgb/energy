#!/bin/bash
# run_chime_repeat.sh - Run E_iso function analysis on CHIME repeating FRBs
# Usage: bash run_chime_repeat.sh

export OMP_NUM_THREADS=4

echo "=== CHIME/FRB Repeating FRB E_iso Function Analysis ==="
echo "Using only repeating bursts from CHIME Catalog 1 & 2"

for fgtype in ALG_NE2001 ALG_YMW16
do
    echo "Running repeaters with galaxy type: ${fgtype}"
    # Catalog 1 repeaters
    echo "mpiexec -n 4 python3 nest_samp.py -fc chime_cat1.json -g ${fgtype} -upper 1 -repeaters -o cat1_repeaters_${fgtype}_upper"
    mpiexec -n 4 python3 nest_samp.py -fc chime_cat1.json -g ${fgtype} -upper 1 -repeaters -o cat1_repeaters_${fgtype}_upper

    # Catalog 2 repeaters
    echo "mpiexec -n 4 python3 nest_samp.py -fc chime_cat2.json -g ${fgtype} -upper 1 -repeaters -o cat2_repeaters_${fgtype}_upper"
    mpiexec -n 4 python3 nest_samp.py -fc chime_cat2.json -g ${fgtype} -upper 1 -repeaters -o cat2_repeaters_${fgtype}_upper
done

echo "=== Done ==="
