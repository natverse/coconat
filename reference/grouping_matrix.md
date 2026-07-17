# Sparse 0/1 aggregation matrix mapping ids to groups

Construct a sparse matrix with one row per id and one column per group,
containing a 1 where an id belongs to a group. Right-multiplying a
connectivity matrix by `grouping_matrix(colnames(M), group)` sums the
columns of `M` within each group (e.g. collapsing partner neurons to
cell types).

## Usage

``` r
grouping_matrix(ids, group)
```

## Arguments

- ids:

  Character vector (or coercible) of identifiers, typically the row or
  column names of a matrix to be aggregated.

- group:

  Grouping labels. Either a vector parallel to `ids` or a *named* vector
  that will be looked up by `ids`. Ids whose group is `NA` are dropped
  (contribute no 1s), so they vanish from any aggregation.

## Value

A sparse [`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html) of
dimension `length(ids)` x `number of groups` with ids as row names and
group labels as column names.

## See also

[`effective_connectivity`](https://natverse.org/coconat/reference/effective_connectivity.md)

## Examples

``` r
library(Matrix)
M <- Matrix(c(1,0,2, 0,3,0), nrow=2, byrow=TRUE,
  dimnames=list(c('a','b'), c('x','y','z')))
# collapse columns x,y -> g1 and z -> g2
G <- grouping_matrix(colnames(M), c(x='g1', y='g1', z='g2'))
M %*% G
#> 2 x 2 Matrix of class "dgeMatrix"
#>   g1 g2
#> a  1  2
#> b  3  0
```
