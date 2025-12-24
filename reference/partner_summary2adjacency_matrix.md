# Convert a partner summary table into an adjacency matrix

Convert a partner summary table into an adjacency matrix

## Usage

``` r
partner_summary2adjacency_matrix(
  x,
  sparse = TRUE,
  inputcol = "pre_id",
  outputcol = "post_id",
  inputids = NULL,
  outputids = NULL,
  standardise_input = T
)
```

## Arguments

- x:

  dataframe as produced by `flywire_partner_summary`,
  `neuprint_connection_table` or equivalent.

- sparse:

  Whether to return a sparse matrix (default `TRUE` in case you are
  making a big one)

- inputcol, outputcol:

  Character vector specifying the columns containing input and output
  identifiers. See **details** section for more information.

- inputids, outputids:

  Optional vectors of input/output ids to ensure that these are present
  in the output matrix. Alternatively these may contain a function that
  takes the dataframe `x` as input and returns a grouping vector. See
  **details** section for more information.

- standardise_input:

  whether to standardise the column names/types in the input dataframe.
  The default should work for flywire `fafbseg` and `neuprintr` input
  and ensure that we identify appropriate `pre_id`/`post_id` columns.

## Value

A matrix with inputs as rows and outputs (downstream neurons) as columns

## Details

The `inputcol` and `outputcol` arguments can name columns containing
other values besides the unique numerical identifiers for neurons. For
example you can refer to a cell type column, thereby generating a
*grouped* connectivity matrix. This is very useful for bringing together
neurons with similar connectivity patterns across brain hemispheres and
individuals.

Passing a function to `inputids` and/or `outputids` allows partner
neurons to be grouped with maximum flexibility. The input to the
function will be the dataframe `x` (after standardisation if this has
been requested). The output must be a single vector which can be
interpreted as a factor to group partner neurons.

## Examples

``` r
if (FALSE) { # \dontrun{
da2ds=neuprintr::neuprint_connection_table('DA2_lPN', details=TRUE, partners='out', threshold=5)
am=partner_summary2adjacency_matrix(da2ds, inputcol = 'bodyid', outputcol = 'partner')
library(Matrix)
image(am, ylab='DA2 PNs', xlab='outputs')

amg=partner_summary2adjacency_matrix(da2ds, inputcol = 'bodyid', outputcol = 'type')
image(amg, ylab='DA2 PNs', xlab='output cell types')
} # }
```
