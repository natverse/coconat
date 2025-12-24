# identifiers-and-datasets

``` r
library(coconat)
```

## datasets

coconat needs to be able to handle a range of datasets. For each dataset
we record both a long name and a short name (2 characters \[a-z0-9\])
that we use to refer to them.

You can tell coconat about a new dataset by doing:

``` r
register_dataset('flywire', shortname = 'fw', species = 'Drosophila melanogaster', sex = 'F')
```

## identifiers and keys

Each dataset will have unique identifiers for the neurons that it
contains. All connectomics data sources that we currently work with have
identifiers within the 64 bit integer range, but there is no reason that
they cannot contain additional characters. To ensure interoperability,
we will

1.  default to handling all identifiers as character vectors (strings)
    even if the native form in a package providing access to a dataset
    is some king of numeric.
2.  introduce the notion of a *key* in the form `"<dataset>:<id>"`

The second point is vital since the same numeric identifier can be used
in different datasets to refer to different neurons.

We provide a number of functions to assist with this. At the most basic
level
[`id2char()`](https://natverse.org/coconat/reference/id2char.md)will
convert numeric ids to a character representation taking care of issues
like the fact that `as.character(1e5)`=1e+05 which can cause trouble.

``` r
id2char(1000)
#> [1] "1000"
id2char(1e5)
#> [1] "100000"
id2char(NA)
#> [1] NA
id2char(bit64::as.integer64('12345678901234'))
#> [1] "12345678901234"
```
