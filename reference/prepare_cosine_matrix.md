# Cosine matrix utility functions

These functions are intended for use by package authors rather than end
users.

## Usage

``` r
prepare_cosine_matrix(
  x,
  partners = c("inputs", "outputs"),
  action = c("zero", "drop")
)
```

## Arguments

- x:

  A matrix or a named list of input/output matrices

- partners:

  Whether to select input or output matrices when both are available

- action:

  Whether to zero out or drop any NA values in the cosine matrix (these
  may be present when some columns have no entries)

## Value

A matrix. When both inputs and outputs are used these will be weighted
by the total number of input and output synapses.
