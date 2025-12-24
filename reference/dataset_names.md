# Return dataset names either all or those matching a query

Return dataset names either all or those matching a query

## Usage

``` r
dataset_names(
  query = NULL,
  return.short = FALSE,
  match = c("both", "long", "short"),
  namespace = "default"
)
```

## Arguments

- query:

  A character vector partially matched against dataset names

- return.short:

  Whether to return the long or short name

- match:

  Whether the query should match against long or short forms of the
  dataset name.

- namespace:

  Optional character vector specifying a namespace used to organise
  datasets (advanced use only).

## Value

A character vector of names

## Examples

``` r
dataset_names()
#> character(0)
```
