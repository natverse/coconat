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
  used for testing purposes.

- ...:

  Additional named arguments specifying properties of the dataset

## Value

No return value. Called for its side effect.

## Examples

``` r
if (FALSE) { # \dontrun{
# partial example. metafun and partnerfun are pretty important to specify
register_dataset("flywire", shortname='fw', species='Drosophila melanogaster', sex='F', age="adult")
} # }
```
