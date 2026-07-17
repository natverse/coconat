# Effective connectivity through multi-step pathways

Compute an estimate of the effective connectivity between the inputs of
the first matrix and the outputs of the last matrix in a chain, passing
through one or more intermediate neuron layers. This follows the
standard approach of input-normalising each step (so the inputs to every
postsynaptic cell sum to 1) and then multiplying the matrices together.

## Usage

``` r
effective_connectivity(matrices, group = NULL, normalise = TRUE)
```

## Arguments

- matrices:

  A list of sparse (or dense) adjacency matrices to chain, each with
  presynaptic neurons as rows and postsynaptic neurons as columns and
  with row/column names to align successive matrices.

- group:

  Optional grouping for the columns (outputs) of the final matrix,
  passed to
  [`grouping_matrix`](https://natverse.org/coconat/reference/grouping_matrix.md) -
  either parallel to the final columns or a named vector indexed by
  them.

- normalise:

  Whether to column-normalise (input-normalise) each matrix with
  [`colScaleM`](https://natverse.org/coconat/reference/colScaleM.md)
  before multiplication (default `TRUE`).

## Value

A sparse matrix of effective connectivity with the inputs of the first
matrix as rows and the outputs (or output groups) of the last matrix as
columns.

## Details

Each supplied matrix must be oriented with presynaptic (input) neurons
as rows and postsynaptic (output) neurons as columns (the convention of
[`partner_summary2adjacency_matrix`](https://natverse.org/coconat/reference/partner_summary2adjacency_matrix.md)).
Consecutive matrices are chained by name: the columns of matrix `k` are
matched to the rows of matrix `k+1`, and any output of matrix `k` that
is absent from the rows of matrix `k+1` contributes a zero (dead-end)
path.

When `normalise=TRUE` (the default) every matrix is column-normalised
with [`colScaleM`](https://natverse.org/coconat/reference/colScaleM.md)
*before* multiplication, implementing the assumption that an interneuron
conveys information about its inputs in proportion to their synaptic
weight. This is the procedure described by Schlegel et al. (2021)
[doi:10.7554/eLife.62576](https://doi.org/10.7554/eLife.62576) .

`group` collapses the final (output) dimension to groups such as cell
types. Because grouping is applied *after* normalisation and
multiplication, and matrix multiplication is associative, grouping the
final product is identical to grouping the last matrix's columns before
the earlier multiplications - the per-neuron normalisation is always
preserved.

**Iterative use**. Passing every hop as one list is convenient but
materialises all the adjacency matrices at once, which can be a memory
hog for deep or highly connected walks. The more scalable (and often
more common) pattern is to advance one hop at a time, keeping only the
running product and the current hop in memory: normalise each new hop
yourself with
[`colScaleM`](https://natverse.org/coconat/reference/colScaleM.md) and
call `effective_connectivity(list(running, step), normalise = FALSE)`.
Here `normalise = FALSE` is essential - `running` is already a product
of normalised matrices and must not be rescaled, while `step` you have
normalised in advance. This also lets you inspect, threshold or prune
the intermediate frontier between hops (e.g. dropping weakly connected
partners before fetching the next layer, as `coconatfly`'s multihop
clustering does) rather than committing to the whole chain up front.

## References

Schlegel et al. (2021) *eLife*
[doi:10.7554/eLife.62576](https://doi.org/10.7554/eLife.62576)

## See also

[`colScaleM`](https://natverse.org/coconat/reference/colScaleM.md),
[`grouping_matrix`](https://natverse.org/coconat/reference/grouping_matrix.md),
[`partner_summary2adjacency_matrix`](https://natverse.org/coconat/reference/partner_summary2adjacency_matrix.md)

## Examples

``` r
library(Matrix)
set.seed(42)
# query (Q) -> interneurons (I) -> targets (T)
A1 <- Matrix(rbinom(6,5,0.4), nrow=2,
  dimnames=list(c('q1','q2'), c('i1','i2','i3')))
A2 <- Matrix(rbinom(9,5,0.4), nrow=3,
  dimnames=list(c('i1','i2','i3'), c('t1','t2','t3')))
effective_connectivity(list(A1, A2))
#> 2 x 3 Matrix of class "dgeMatrix"
#>           t1     t2        t3
#> q1 0.4583333 0.4375 0.4642857
#> q2 0.5416667 0.5625 0.5357143

# equivalent iterative form: advance one hop at a time, normalising each new
# hop yourself and keeping only the running product (memory friendly, and you
# can prune the frontier between hops)
running <- colScaleM(A1)
running <- effective_connectivity(list(running, colScaleM(A2)),
                                  normalise = FALSE)
running
#> 2 x 3 Matrix of class "dgeMatrix"
#>           t1     t2        t3
#> q1 0.4583333 0.4375 0.4642857
#> q2 0.5416667 0.5625 0.5357143
```
