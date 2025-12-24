# Convert any id into a character vector

Convert any id into a character vector

## Usage

``` r
id2char(x)
```

## Arguments

- x:

  A character or numeric vector of ids

## Value

A character vector

## Examples

``` r
id2char(1000)
#> [1] "1000"
id2char("1000")
#> [1] "1000"
id2char(1e5)
#> [1] "100000"
```
