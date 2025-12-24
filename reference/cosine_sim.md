# Efficient calculation of cosine similarity for sparse or dense matrices

Efficient calculation of cosine similarity for sparse or dense matrices

## Usage

``` r
cosine_sim(x, sparse = FALSE, transpose = FALSE)
```

## Arguments

- x:

  A data matrix suitable for clustering

- sparse:

  Whether to return a sparse or dense matrix (default dense)

- transpose:

  When `F` (the default) calculates the cosine distance between columns.
  When `T` calculates the distance between rows.

## Value

A square matrix, dense unless `sparse=TRUE` and `x` is sparse.

## Details

This is much faster than e.g. `lsa::cosine` with the added benefit that
it works for sparse input matrices

## Examples

``` r
da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
cosine_sim(am)
#>            1796817841 1797505019 1796818119 1827516355  818983130
#> 1796817841 1.00000000 0.09631851 0.00000000  0.0000000 0.03509041
#> 1797505019 0.09631851 1.00000000 0.00000000  0.0000000 0.09507157
#> 1796818119 0.00000000 0.00000000 1.00000000  0.0000000 0.07046832
#> 1827516355 0.00000000 0.00000000 0.00000000  1.0000000 0.20027384
#> 818983130  0.03509041 0.09507157 0.07046832  0.2002738 1.00000000

if (FALSE) { # \dontrun{
fam_pnkc2=flywire_adjacency_matrix2(
  inputids = 'class:ALPN_R', outputids = 'class:Kenyon_Cell_R',
  sparse=T , threshold = 2)
kckc.cos=cosine_sim(fam_pnkc2)
pnpn.cos=cosine_sim(fam_pnkc2, transpose=T)
} # }
```
