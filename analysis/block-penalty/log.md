
simulations to determine optimal block penalty

1. block penalty for background rates (100 for Singhal)
2. block penalty for hotspot detection (5 for Singhal)

## 25/3/2019

- generate 1 Mb for each of 19 individuals
- repeat simulation for different background recombination rates
    - Singhal - 0.0001, 0.001, 0.01, 0.1, 1, 2.5
- values based off preliminary recombination maps - we should use our own prelim LDhelmet data for this

what is the range of prelim values we have?

checking range at block = 50 (`data/ldhelmet_files/genome50`)

```R
> fnames <- dir_ls('.', regexp = '^chromosome_[0-9]{1,2}_50.txt$')
> f_cols <- c('left_snp', 'right_snp', 'mean', 'p025', 'p975')
> d <- map_dfr(fnames, read_delim, delim = ' ', skip = 3,
+ col_types = cols(), col_names = f_cols, .id = 'name')
> dim(d)
[1] 7333919       6
> d
# A tibble: 7,333,919 x 6
                   name left_snp right_snp       mean       p025      p975
                  <chr>    <int>     <int>      <dbl>      <dbl>     <dbl>
 1 chromosome_10_50.txt      796       801 0.00120500 0.00101350 0.0012328
 2 chromosome_10_50.txt      801       814 0.00100800 0.00100000 0.0010491
 3 chromosome_10_50.txt      814       823 0.00090795 0.00059584 0.0010000
 4 chromosome_10_50.txt      823       824 0.00100000 0.00100000 0.0010000
 5 chromosome_10_50.txt      824       839 0.00109140 0.00100000 0.0012339
 6 chromosome_10_50.txt      839       840 0.00100000 0.00100000 0.0010000
 7 chromosome_10_50.txt      840       841 0.00107200 0.00100000 0.0011541
 8 chromosome_10_50.txt      841       877 0.00130660 0.00100000 0.0016031
 9 chromosome_10_50.txt      877       880 0.00104530 0.00100000 0.0010948
10 chromosome_10_50.txt      880       887 0.00118030 0.00100000 0.0019334
# ... with 7,333,909 more rows
> summary(d$mean)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
0.0000696 0.0018844 0.0029856 0.0034082 0.0044192 2.1555000
```

these don't look too different from the Singhal values at all -
we'll stick to the same ones. 

## 9/4/2019

theta and mutation matrix

theta = 0.03 (Rob 2015)

mutation matrix = ?

looking through the files in `mutability/data/mutation_prediction` - could
get transition probs from this? check with Rob

other things to do:
- get on rewriting that intro!
- change CO density to crossovers

## 11/4/2019

today: mutation matrix estimation

from LDhelmet manual:

The mutation matrix specifies the transition probabilities from an allele to another allele. 
The order of alleles is (A, C, G, T) for both the rows and the columns. 
The file format is a line for every row in the mutation matrix, 
where each line contains 4 transition probabilities (floating point numbers).

For example,

0.1 0.2 0.3 0.4
0.4 0.3 0.2 0.1
0.2 0.2 0.2 0.4
0.1 0.1 0.1 0.7

quick R script to tally:

```R
library(readr)
library(dplyr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

fname <- args[1]
outfile <- args[2]

d <- read_tsv(fname, 
    col_types = cols_only(
        anc_base = 'c',
        mut_base = 'c',
        mut_type = 'c')) %>%
  filter(mut_type == 'SNP') %>%
  select(contains('base')) %>%
  group_by(anc_base, mut_base) %>%
  tally()

write_tsv(d, outfile)
```

in bash:

```bash
for i in {1..17}; do
    Rscript analysis/block-penalty/mutation_counts.R \
    data/mutation-matrix/chromosome_${i}.mutation_prediction.txt \
    data/mutation-matrix/transition-counts/chromosome_${i}.counts;
done
```

## 18/4/2019

things we need for this estimation:
- f_i for all i in [ACGT] (# of bases)
- f_ij for all bases where i != j (in count files above)
- M - frequency of mutation at 'most mutable base'
    - needs to be estimated from f_ij values

combining the files in `mutation-matrix/transition-counts/` above:

```R
library(dplyr)
library(magrittr)
library(purrr)
library(fs)

fnames <- dir_ls('.') # while in dir
d <- map_dfr(fnames, read_tsv, col_types = cols())
d %<>% group_by(anc_base, mut_base) %>% summarise(n = sum(n))
write_tsv(d, 'all_counts.tsv')
```

what is the 'most mutable base'?

```R
> d %>% group_by(anc_base) %>% summarise(m = sum(n)) %>% arrange(desc(m))
# A tibble: 4 x 2
  anc_base     m
     <chr> <int>
1        C  2356
2        G  2281
3        T   561
4        A   512
```

looks to be C. 

what about the # of each base? ie f_i for all bases?

quick python script:

```python
>>> from Bio import SeqIO
>>> chrom_names = ['chromosome_' + str(i) for i in range(1,18)]
>>> chroms = dict.fromkeys(chrom_names, dict.fromkeys(['A', 'C', 'G', 'T']))
>>> fname = 'chlamy.5.3.w_organelles_mtMinus.fasta'
>>> from tqdm import tqdm
>>> chroms = dict.fromkeys(chrom_names)
>>> for record in tqdm(SeqIO.parse(fname, 'fasta')):
  2     if 'chromosome' in record.id:
  3         chroms[record.id] = {}
  4         for base in ['A', 'T', 'C', 'G']:
  5             chroms[record.id][base] = record.seq.count(base)
>>> for chrom in chroms:
  2     print(chroms[chrom])
{'A': 587277, 'T': 584704, 'G': 1029065, 'C': 1025169}
{'A': 1376787, 'T': 1373629, 'G': 2434482, 'C': 2457889}
{'A': 720370, 'T': 714444, 'G': 1298441, 'C': 1289548}
{'A': 1241126, 'T': 1250242, 'G': 2217532, 'C': 2228488}
{'A': 1595704, 'T': 1589010, 'G': 2885377, 'C': 2867247}
{'A': 1115178, 'T': 1121290, 'G': 2027597, 'C': 2021374}
{'A': 1375486, 'T': 1383567, 'G': 2452827, 'C': 2453964}
{'A': 652437, 'T': 656057, 'G': 1150839, 'C': 1144622}
{'A': 716502, 'T': 713839, 'G': 1245855, 'C': 1246372}
{'A': 1354048, 'T': 1358236, 'G': 2459037, 'C': 2455849}
{'A': 1566973, 'T': 1564447, 'G': 2814011, 'C': 2816441}
{'A': 310969, 'T': 313726, 'G': 532922, 'C': 533403}
{'A': 1611074, 'T': 1607956, 'G': 2843698, 'C': 2849477}
{'A': 1697111, 'T': 1689021, 'G': 3034091, 'C': 3047984}
{'A': 901319, 'T': 909454, 'G': 1612824, 'C': 1621231}
{'A': 858184, 'T': 860355, 'G': 1533324, 'C': 1527248}
{'A': 1166521, 'T': 1170171, 'G': 2092101, 'C': 2092599}
>>> total_counts = dict.fromkeys(['A', 'T', 'C', 'G'], 0)
>>> for key in tqdm(total_counts):
  2     for chrom in chroms:
  3         total_counts[key] += chroms[chrom][key]
100%|███████████████████████████████████████| 4/4 [00:00<00:00, 57653.66it/s]
>>> total_counts
{'A': 18847066, 'T': 18860148, 'G': 33664023, 'C': 33678905}
```

taking stock - so far we have f_i for all i (A, C, G, T)
as well as all f_ij for all j (A, C, G, T where j != i)

next up we need M - the highest frequency of mutation away
from a site. 

```python
>>> mut_counts = {'C': 2356, 'G': 2281, 'T': 561, 'A': 512}
>>> freqs = {}
>>> for base in mut_counts:
  2     freqs[base] = mut_counts[base] / total_counts[base]
>>> freqs
{'T': 2.974525968725166e-05, 'C': 6.995476842254818e-05, 'G': 6.775779591167699e-05, 'A': 2.7166032102821734e-05}
>>> max(freqs.values())
6.995476842254818e-05 # C is the most mutable
```

so M = 6.995 * 10^-5

## 19/4/2019

the mutation matrix should be in the order A, C, G, T

```python
mut_mat = {}
for base in ['A', 'T', 'C', 'G']:
    mut_mat[base] = {}
total_counts = {'A': 18847066, 'T': 18860148, 
                'G': 33664023, 'C': 33678905}
M = 6.995e-05


fname = 'data/mutation-matrix/all_counts.tsv'

# i != j
with open(fname, 'r') as f:
    for line in f:
        if not line.startswith('anc_base'):
            anc, mut, n = line.rstrip('\n').split('\t')
            n = int(n)
            mut_mat[anc][mut] = n / (total_counts[anc] * M)

# i == j
for base in ['A', 'T', 'C', 'G']:
    mut_mat[base][base] = 1 - sum(mut_mat[base].values())

```

this works out nicely:

```python
>>> for i in ['A', 'C', 'G', 'T']:
  2     for j in ['A', 'C', 'G', 'T']:
  3         print(i, j, round(mut_mat[i][j], 2))
A A 0.61
A C 0.11
A G 0.19
A T 0.09
C A 0.26
C C -0.0
C G 0.22
C T 0.52
G A 0.5
G C 0.22
G G 0.03
G T 0.25
T A 0.11
T C 0.2
T G 0.11
T T 0.57
```

and so the file will be (in the order ACGT)

```
0.61 0.11 0.19 0.09
0.26 0 0.22 0.52
0.5 0.22 0.03 0.25
0.11 0.2 0.11 0.57
```

alright - now for the sims:

- generate 1 Mb for each of 19 individuals
- repeat simulation for different background recombination rates
    - Singhal - 0.0001, 0.001, 0.01, 0.1, 1, 2.5
- theta = 0.03
- mutation matrix = data/mut_mat

sim structure: pseudocode

```
for rho in [0.0001, 0.001, 0.01, 0.1, 1, 2.5]:
    generate 19 indivs, 1 Mb of sequence each, rcmb = $rho
    for hotspot_heat in [10x, 20x, 40x, 60x]:
        add 2 hotspots w/ relative heat $hotspot heat

repeat each parameter set 12 times (we'll just do 10)

remove singletons from simulated haps

for block in 5, 10, 50, 100, 500:
    run LDhelmet on haps w/ block = $block
```

from Singhal's supplementary:

For these simulation results, we asked two questions. First, which block penalty provides the
most accurate general picture of recombination rates? To address this question, we calculated how
much inferred rates from LDhelmet deviated from known recombination rates across the entire
length of the 1 Mb of simulated sequence. These results suggested that block penalty 100 provided the most accurate estimation of background recombination rate across a range of recombination rates (Fig. S30). Second, we asked which block penalty provides the best power to identify
hotspots. To answer this question, we identified putative hotspots in the inferred recombination
maps, identifying regions 2 kb or greater that had 5× or greater recombination rate than their 80
kb of surrounding sequence, as was done with the empirical data. From this, we calculated the
power to identify hotspots at each given parameter set, finding that block penalty 5 provided the
most power (Fig. S31).

sample macs run - from the docs:

./macs 100 1e6 -T -t .001 -r .001 -h 1e2 -R example_input/hotspot.txt -F example_input/ascertainment.txt 0 2>trees.txt | ./msformatter > haplotypes.txt

The command line arguments above says:
Simulate 100 sequences on a region 1e6 basepairs long.  The per base pair mutation and recombination rate scaled at 4N is .001.  The h parameter approximates the Markovian order by instructing the program to include all geneologies from the current point to 1e2 base pairs to the left if one considers simulation proceeding from the left end to the right end of the entire sequence.  Any branches unique to the geneology that is beyond 1e2 base pairs is pruned from the ARG. -T tells MACS to output the local trees in Newick format similar to MS output.

in our command:

- 1 Mb
- theta (-t) = 0.03
- rho (-r) = 0.0001
- hotspots -> -R file
    - hotspots are 2 kb
    - evenly spread - `starts = np.linspace(50e3, 1e6 - 50e3, 10)`
        - ends are each value + 2000
        - need to divide start and end by 1e6
        - Singhal script seems to lock diffs to 5

## 23/4/2019

let's try a 'naive' macs run:

first, making the tab separated hotspot file (`test_hotspot`)

```
0   0.3 0.57
```

and then

```bash
./bin/macs 19 1e6 -t 0.03 -r 0.0001 -R data/macs-runs/test_hotspot | bin/msformatter
```

pulling in the Singhal script - need to figure out the equilibrium base freqs and
the modified mutation matrix (ie for all i != j, Pij := Pij / (1 - Pii))


total freqs:

```python
>>> x = {'A': 18847066, 'T': 18860148, 'G': 33664023, 'C': 33678905}
>>> total = sum(x.values())
>>> total
105050142
>>> for k, v in x.items():
  2     print(k, round(v / total, 3))
A 0.179
C 0.321
T 0.18
G 0.32
```

remainder of the matrix - getting relative, non-diag frequencies: 

```
0.61 0.11 0.19 0.09
0.26 0 0.22 0.52
0.5 0.22 0.03 0.25
0.11 0.2 0.11 0.57
```

or alternatively:

```
d = {
'AA': 0.61,
'AC': 0.11,
'AG': 0.19,
'AT': 0.09,
'CA': 0.26,
'CC': -0.0,
'CG': 0.22,
'CT': 0.52,
'GA': 0.5,
'GC': 0.22,
'GG': 0.03,
'GT': 0.25,
'TA': 0.11,
'TC': 0.2,
'TG': 0.11,
'TT': 0.57
}
```

and so:

```python
>>> d = {
  2 'AA': 0.61,
  3 'AC': 0.11,
  4 'AG': 0.19,
  5 'AT': 0.09,
  6 'CA': 0.26,
  7 'CC': -0.0,
  8 'CG': 0.22,
  9 'CT': 0.52,
 10 'GA': 0.5,
 11 'GC': 0.22,
 12 'GG': 0.03,
 13 'GT': 0.25,
 14 'TA': 0.11,
 15 'TC': 0.2,
 16 'TG': 0.11,
 17 'TT': 0.57
 18 }
>>> d
{'TA': 0.11, 'TG': 0.11, 'AT': 0.09, 'GG': 0.03, 'GC': 0.22, 'GT': 0.25, 'CA': 0.26, 'CC': -0.0, 'CT': 0.52, 'CG': 0.22, 'TT': 0.57, 'AG': 0.19, 'TC': 0.2, 'AC': 0.11, 'AA': 0.61, 'GA': 0.5}
>>> non_diag_sums = dict.fromkeys(['A', 'C', 'G', 'T'], 0.0)
>>> for k, v in d.items():
  2     anc, mut = k[0], k[1]
  3     if anc != mut:
  4         non_diag_sums[anc] += v
>>> non_diag_sums
{'A': 0.39, 'C': 1.0, 'G': 0.97, 'T': 0.42000000000000004}
>>> for k, v in d.items():
  2     anc, mut = k[0], k[1]
  3     if anc != mut:
  4         d[k] = v / non_diag_sums[anc]
>>> d
{'TA': 0.26190476190476186, 'TG': 0.26190476190476186, 'AT': 0.23076923076923075, 'GG': 0.03, 'GC': 0.2268041237113402, 'GT': 0.2577319587628866, 'CA': 0.26, 'CC': -0.0, 'CT': 0.52, 'CG': 0.22, 'TT': 0.57, 'AG': 0.48717948717948717, 'TC': 0.47619047619047616, 'AC': 0.28205128205128205, 'AA': 0.61, 'GA': 0.5154639175257733}
>>> # checking
>>> d['TA'] + d['TC'] + d['TG']
0.9999999999999999
>>> for k, v in d.items():
  2     anc, mut = k[0], k[1]
  3     if anc != mut:
  4         print(anc, mut, round(v, 3))
T A 0.262
T G 0.262
A T 0.231
G C 0.227
G T 0.258
C A 0.26
C T 0.52
C G 0.22
A G 0.487
T C 0.476
A C 0.282
G A 0.515
```

these can now be manually added in to the script - and the script is now ready to run in full!

```bash
time python3.5 analysis/block-penalty/simulations_hotspot_power.py \
--outdir data/macs-runs
```

## 29/4/2019

so that took a while:

```bash
real    2157m56.549s
user    2031m27.540s
sys     128m59.895s
```

but now we're good to start some LDhelmet runs - there are 100 haps to run LDhelmet on

this is probably best done with a shell script that iterates through the haps,
runs the indiv LDhelmet steps on them, and deletes temp files as soon as they're no longer needed

```bash
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

    time ./bin/ldhelmet table_gen \    
    --num_threads 10 \  
    --conf_file data/macs-runs/ldhelmet/${outdir}/${base}.conf \   
    --theta 0.03 \  
    --rhos 0.0 0.1 10.0 1.0 100.0 \ 
    --output_file data/macs-runs/ldhelmet/${outdir}/${base}.lk > table_gen 2> table_gen2   

    rm -v table_gen* 

    time ./bin/ldhelmet pade \ 
    --num_threads 10 \  
    --conf_file data/macs-runs/ldhelmet/${outdir}/${base}.conf \   
    --theta 0.03 \  
    --output_file data/macs-runs/ldhelmet/${outdir}/${base}.pade   

     time ./bin/ldhelmet rjmcmc \   
    --num_threads 10 \  
    --window_size 50 \  
    --seq_file ${filename} \    
    --lk_file data/macs-runs/ldhelmet/${outdir}/${base}.lk \   
    --pade_file data/macs-runs/ldhelmet/${outdir}/${base}.pade \   
    --num_iter 1000000 \    
    --burn_in 100000 \  
    --block_penalty ${block} \    
    --mut_mat_file data/mut_mat
    --output_file data/recombination-ldhelmet/intermediate-files/${base}.post   

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
```

followed by

```bash
time bash analysis/block-penalty/ldhelmet_block.sh 50
```

to start at block = 50. 


## 5/5/2019

so rho = 0.1 is agonizingly slow - while this is going
on, let's queue a temp script up for block = 5 at rho = 0.0001

```bash
# same as above file, but fixed parameters
# block = 5, rho = 0.0001
time bash ldhelmet_block_temp_2.sh
```

## 6/5/2019

so this went pretty fast by comparison (done in 7 hours)
which tells me LDhelmet likes lower rho values

now getting rho = 0.001 going for block = 5

update: even faster - done in 144 min

now on rho = 0.01 for block = 5

update: 149 min for that - running rho = 0.1 and 1.0 
simultaneously - it seems the MCMC step is quite fast at block = 5

## 7/5/2019

so both of these were done in the wee hours of the morning

running the final one (rho = 2.5) for block = 5
as well as the first one (rho = 0.0001) for block = 10

## 8/5/2019

the rho = 1.0 block = 50 script is finally done
after 6400 m! (that's 4 days!) - wonder what slowed
that down so much

meanwhile - block 5 is completely done for all rhos

continuing block = 10 (0.0001 is done) over 0.001 and 0.01 and
also queuing up the final block = 50 (rho = 2.5)

update - the block 10s are done? queuing up rho = 0.1 and 1.0

## 9/5/2019

queuing up rho = 2.5 for block = 10 - and with that, block 10 is done
(as will be block 50 when the rho = 2.5 run is done)

last but not least, block = 100. queuing up rho = 0.0001

also reminder to self to run rho = 2.5 replicate 100 for block 50 separately! 
that was not included in the temp script somehow

## 10/5/2019

running rho = 2.5 replicate 100 on its own

rho = 0.0001 for block = 100 is done (took 14 hours) - will
queue up the others after the current three jobs are done

## 11/5/2019

queuing up rho = 0.001 for block = 100

block = 10 is done - queuing up rho = 0.01 for block = 100

## 12/5/2019

0.001 and 0.01 for block = 100 are done

queuing up 0.1 and 1.0 next - almost done!

interesting though that block = 50 takes the longest...

## 15/5/2019

block 50 is done (with rho = 2.5 taking 6 straight days)

queuing up 100-104 and 105-109 for rho = 2.5 and block = 100 separately - and
then we'll be done!

## 25/5/2019

today: adapt `find_hotspots.py` for use with simulation outputs

this should be a more generalized script - will also be using
it to summarise the actual recombination estimates when LDhelmet
is run on chlamy itself

need to compare:
- accuracy of mean background rho estimate across different block penalties
    - try this including + excluding hotspots
- accuracy of hotspot estimates

grabbing a quick 'test dataset':

```
cp -v data/macs-runs/ldhelmet/block_50/finals/haplo_rho0.01_100.txt .
mv -v haplo_rho0.01_100.txt rho_test.txt
```

update: let's just use `find_hotspots.py` - after adding in chlamy chromosome
lengths and removing some of the hardcoded paths it should be good
for use as is

## 26/5/2019

script to run `find_hotspots.py` over all sims

added a chr called 'sim' that's 1e6 bp in length

initial run over 0.0001, block = 5:

```bash
mkdir -p data/macs-runs/ldhelmet_2k
mkdir -p data/macs-runs/ldhelmet_2k/block_5

for run in {0..9}; do
    time python3.5 analysis/block-penalty/find_hotspots.py \
    --input data/macs-runs/ldhelmet/block_5/finals/haplo_rho0.0001_10${run}.txt \
    --out data/macs-runs/ldhelmet_2k/block_5/haplo_rho0.0001_10${run}.txt \
    --chr sim \
    --block 2000 \
    --flank 40000 ;
done

```

## 27/5/2019

looks good - doing the same for the other blocks:

```bash
for block in 5 10 50 100; do
    mkdir -p data/macs-runs/ldhelmet_2k/block_${block};
    for rho in 0.0001 0.001 0.01 0.1 1.0 2.5; do
        echo "currently on rho ${rho} for block ${block}";
        for run in {0..9}; do
            python3.5 analysis/block-penalty/find_hotspots.py \
            --input data/macs-runs/ldhelmet/block_${block}/finals/haplo_rho${rho}_10${run}.txt \
            --out data/macs-runs/ldhelmet_2k/block_${block}/haplo_rho${rho}_10${run}.txt \
            --chr sim \
            --block 2000 \
            --flank 40000 ;
        done
    done
done

```

## 29/5/2019

looking at these files in `analysis/block-penalty/comparisons.Rmd`

block penalty = 100 seems to be the best bet for background rho!

looks like power is relatively poor with 2 kb windows and 40 kb flanks though - below 0.5 for all rates/blocks

## 30/5/2019

creating a shell script (`summarise_hotspots.sh`) to run different combinations of windows/flanks:

```bash
w=$1
flank=$2

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
```

also renaming `ldhelmet_2k` to `ldhelmet_2k_40k` to include flank size as well

run with 4 kb windows and 40 kb flank sizes:

```bash
time bash analysis/block_penalty/summarise_hotspots.sh 4000 40000
```


















