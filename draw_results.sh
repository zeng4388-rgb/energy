#!/bin/bash
# draw_results.sh - Plot posterior distributions from MultiNest output
# Usage: bash draw_results.sh

mkdir -p plots/samp plots/simu

echo "=== Plotting Posterior Distributions ==="

# CHIME Catalog 1 results
for fgtype in ALG_NE2001_upper ALG_YMW16_upper ETG_upper LTG_NE2001_upper LTG_YMW16_upper
do
    if [ -f "./nest_out/samp/cat1_${fgtype}" ]; then
        echo "Plotting cat1_${fgtype}..."
        python3 pltpost.py -f ./nest_out/samp/cat1_${fgtype} \
            -o ./plots/samp/cat1_${fgtype}.eps \
            -title "CHIME Cat1 E_iso (${fgtype})" -up 1 -bo 1
    fi
done

# CHIME Catalog 2 results
for fgtype in ALG_NE2001_upper ALG_YMW16_upper
do
    if [ -f "./nest_out/samp/cat2_${fgtype}" ]; then
        echo "Plotting cat2_${fgtype}..."
        python3 pltpost.py -f ./nest_out/samp/cat2_${fgtype} \
            -o ./plots/samp/cat2_${fgtype}.eps \
            -title "CHIME Cat2 E_iso (${fgtype})" -up 1 -bo 1
    fi

    if [ -f "./nest_out/samp/cat2_oneoff_${fgtype}" ]; then
        echo "Plotting cat2_oneoff_${fgtype}..."
        python3 pltpost.py -f ./nest_out/samp/cat2_oneoff_${fgtype} \
            -o ./plots/samp/cat2_oneoff_${fgtype}.eps \
            -title "CHIME Cat2 One-off E_iso (${fgtype})" -up 1 -bo 1
    fi

    if [ -f "./nest_out/samp/cat2_repeaters_${fgtype}" ]; then
        echo "Plotting cat2_repeaters_${fgtype}..."
        python3 pltpost.py -f ./nest_out/samp/cat2_repeaters_${fgtype} \
            -o ./plots/samp/cat2_repeaters_${fgtype}.eps \
            -title "CHIME Cat2 Repeaters E_iso (${fgtype})" -up 1 -bo 1
    fi
done

# Simulation recovery results
for phis in 1e3 1e4
do
    for suffix in "" "_upper"
    do
        if [ -f "./nest_out/simu/simdat_${phis}${suffix}" ]; then
            echo "Plotting simdat_${phis}${suffix}..."
            python3 pltpost.py -f ./nest_out/simu/simdat_${phis}${suffix} \
                -o ./plots/simu/simdat_${phis}${suffix}.eps \
                -title "Mock data (phis=${phis})" -up ${suffix:+1} -bo 1
        fi
    done
done

echo "=== Done ==="
