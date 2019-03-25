
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
