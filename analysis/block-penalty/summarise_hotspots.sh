#!/bin/bash

# summarise LDhelmet data with different combinations
# windowsizes and flank sizes for downstream hotspot detection power analysis

w=$1
flank=$2

echo "running with windowsize ${w} and flank size ${flank}"

sleep 3

for block in 5 10 50 100; do
    mkdir -p data/macs-runs/ldhelmet_${w}_${flank}/block_${block};
    for rho in 0.0001 0.001 0.01 0.1 1.0 2.5; do
        echo "currently on rho ${rho} for block ${block}";
        for run in {0..9}; do
            python3.5 analysis/block-penalty/find_hotspots.py \
            --input data/macs-runs/ldhelmet/block_${block}/finals/haplo_rho${rho}_10${run}.txt \
            --out data/macs-runs/ldhelmet_${w}_${flank}/block_${block}/haplo_rho${rho}_10${run}.txt \
            --chr sim \
            --block ${w} \
            --flank ${flank} ;
        done
    done
done
