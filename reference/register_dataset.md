# Register a dataset for use with coconat and related packages

Register a dataset for use with coconat and related packages

## Usage

``` r
register_dataset(
  name,
  shortname = NULL,
  species = NULL,
  sex = c("F", "M", "H", "U"),
  age = NULL,
  idfun = NULL,
  metafun = NULL,
  partnerfun = NULL,
  namespace = "default",
  inherits = NULL,
  ...
)
```

## Arguments

- name:

  The long name of the dataset; this must be unique

- shortname:

  An abbreviation for the dataset - will be used to construct keys and
  plot labels etc

- species:

  The binomial name for the species

- sex:

  Female, Male, Hermaphrodite or Uncertain

- age:

  A description of the stage (e.g. adult or L1 or P19)

- idfun:

  A function or function name to fetch ids based on a query
  specification.

- metafun:

  A function or function name to fetch metadata about a set of ids for
  this dataset.

- partnerfun:

  A function or function name to fetch connections for ids in this
  dataset.

- namespace:

  Expert use only. Can be used to define separate namespaces across
  which dataset names and keys do not have to be unique. Currently only
  used by the coconatfly package.

- inherits:

  The name of a previously registered dataset. See details.

- ...:

  Additional named arguments specifying properties of the dataset

## Value

No return value. Called for its side effect.

## Details

You can use the `inherits` argument to simplify adding an additional
handler for an existing dataset. For example imagine you have some of
your own annotations that you would like to supplement publicly released
ones for the banc dataset. You can register a new dataset `bancx` and
inherit from the `banc` definition and only replace the `metafun`
argument keeping everything else the same.

## Examples

``` r
if (FALSE) { # \dontrun{
# partial example. metafun and partnerfun are pretty important to specify
register_dataset("flywire", shortname='fw', species='Drosophila melanogaster', sex='F', age="adult")
} # }
```
