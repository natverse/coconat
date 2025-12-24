# Add dataframe columns detailing ordering and groups from a dendrogram

Add dataframe columns detailing ordering and groups from a dendrogram

## Usage

``` r
add_cluster_info(
  df,
  dend,
  h = NULL,
  k = NULL,
  colnames = NULL,
  idcol = c("key", "id", "root_id", "bodyid")
)
```

## Arguments

- df:

  A dataframe e.g. as returned by `coconatfly::cf_meta`

- dend:

  A [`dendrogram`](https://rdrr.io/r/stats/dendrogram.html) or
  [`hclust`](https://rdrr.io/r/stats/hclust.html) object

- h:

  A height to cut the dendrogram

- k:

  A number of clusters to cut the dendrogram

- colnames:

  The names of the two new columns

- idcol:

  The name of the column containing id information. If you do not
  provide this argument then the function will choose the first of three
  default columns present in `df` (with a warning when \>1 column
  present).

## Value

A copy of `df` with one or two extra columns called `dendid` and e.g.
`group_h2`.
