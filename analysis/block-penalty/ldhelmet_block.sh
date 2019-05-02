count=1
filecount=$(ls data/macs-runs/haplotypes | wc -l)
block=$1
mkdir -p data/macs-runs/ldhelmet/block_${block}
mkdir -p data/macs-runs/ldhelmet/block_${block}/finals
outdir=block_$1

for fname in data/macs-runs/haplotypes/*fa; do
    base=$(basename $fname .fa)
    echo "Currently on ${base}"
    echo "File ${count} of ${filecount}"

    time ./bin/ldhelmet find_confs \
    --num_threads 10 \
    --window_size 50 \
    --output_file data/macs-runs/ldhelmet/${outdir}/${base}.conf ${fname}

    sleep 1

    echo "table_gen for ${base}"

    time ./bin/ldhelmet table_gen \
    --num_threads 10 \
    --conf_file data/macs-runs/ldhelmet/${outdir}/${base}.conf \
    --theta 0.03 \
    --rhos 0.0 0.1 10.0 1.0 100.0 \
    --output_file data/macs-runs/ldhelmet/${outdir}/${base}.lk > table_gen 2> table_gen2

    rm -v table_gen*

    sleep 1

    time ./bin/ldhelmet pade \
    --num_threads 10 \
    --conf_file data/macs-runs/ldhelmet/${outdir}/${base}.conf \
    --theta 0.03 \
    --output_file data/macs-runs/ldhelmet/${outdir}/${base}.pade

    sleep 1

    time ./bin/ldhelmet rjmcmc \
    --num_threads 30 \
    --window_size 50 \
    --seq_file ${fname} \
    --lk_file data/macs-runs/ldhelmet/${outdir}/${base}.lk \
    --pade_file data/macs-runs/ldhelmet/${outdir}/${base}.pade \
    --num_iter 1000000 \
    --burn_in 100000 \
    --block_penalty ${block} \
    --mut_mat_file data/mut_mat \
    --output_file data/macs-runs/ldhelmet/${outdir}/${base}.post

    sleep 1

    time ./bin/ldhelmet post_to_text \
    --mean \
    --perc 0.025 \
    --perc 0.50 \
    --perc 0.975 \
    --output_file data/macs-runs/ldhelmet/${outdir}/finals/${base}.txt \
    data/macs-runs/ldhelmet/${outdir}/${base}.post

    sleep 3

    echo "Removing temp files..."
    rm -v data/macs-runs/ldhelmet/${outdir}/${base}.* ;

    (( count ++ ))

done 




