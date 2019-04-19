
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

