.datasets <- new.env()

#' Register a dataset for use with coconat and related packages
#'
#' @param name The long name of the dataset; this must be unique
#' @param shortname An abbreviation for the dataset - will be used to construct
#'   keys and plot labels etc
#' @param sex Female, Male, Hermaphrodite or Uncertain
#' @param species The binomial name for the species
#' @param age A description of the stage (e.g. adult or L1 or P19)
#' @param idfun A function or function name to fetch ids based on a query
#'   specification.
#' @param metafun A function or function name to fetch metadata about a set of
#'   ids for this dataset.
#' @param partnerfun A function or function name to fetch connections for ids in
#'   this dataset.
#' @param namespace Expert use only. Can be used to define separate namespaces
#'   across which dataset names and keys do not have to be unique. Currently
#'   only used for testing purposes.
#' @param ... Additional named arguments specifying properties of the dataset
#'
#' @return
#' @export
#'
#' @examples
register_dataset <- function(name, shortname=NULL, species=NULL,
                             sex=c("F", "M", "H", "U"), age=NULL,
                             idfun=NULL, metafun=NULL, partnerfun=NULL,
                             namespace='default', ...) {

  ns=dataset_namespace(namespace)
  nn=ls(ns)
  if(name %in% nn) {
    warning("Dataset with name: ", name, " already registered. Overwriting!")
    ns[[name]]=NULL
  }

  if(is.null(shortname))
    shortname=unname(abbreviate(name, minlength = 2))

  sex=match.arg(sex)

  if(is.null(idfun))
    idfun <- function(ids, integer64=FALSE) {default_id_fun(ids, metafun = metafun, integer64 = integer64)}

  ns[[name]]=list(
    name=name,
    shortname=shortname,
    species=species,
    sex=sex,
    age=age,
    idfun=idfun,
    metafun=metafun,
    partnerfun=partnerfun,
    ...
  )
  invisible()
}

remove_namespace <- function(namespace) {
  .datasets[[namespace]] <- NULL
  invisible()
}

dataset_namespace <- function(namespace='default') {
  ns=.datasets[[namespace]]
  if(is.null(ns)) {
    ns=new.env()
    .datasets[[namespace]]=ns
  }
  ns
}

#' Return dataset names either all or those matching a query
#'
#' @param query A character vector partially matched against dataset names
#' @param short Whether to match (and return) short names
#' @param namespace
#'
#' @return
#' @export
#'
#' @examples
dataset_names <- function(query=NULL, return.short=FALSE, match=c("both", "long", "short"), namespace='default') {

  ns=dataset_namespace(namespace)
  ln=ls(ns)
  sn=dataset_shortnames(namespace)

  if(is.null(query)) {
    res=if(return.short) unname(sn) else ln
    return(res)
  }
  match=match.arg(match)
  if(match=='both') match=c("long", "short")
  for (m in match) {
    table=if(m=='short') sn else ln
    idx=pmatch(query, table, duplicates.ok = T)
    if(!all(is.na(idx))){
      res <- if(return.short) sn[idx] else ln[idx]
      return(unname(res))
    }
  }
  stop("Unable to match query: ", query, " to known datasets!")
}

dataset_details <- function(query, namespace='default') {
  n=dataset_names(query, namespace = namespace)
  stopifnot(length(n)==1)
  ns=dataset_namespace(namespace = namespace)
  ns[[n]]
}

dataset_shortnames <- function(namespace='default') {
  ns=dataset_namespace(namespace)
  l=sapply(ls(ns), function(n) ns[[n]][['shortname']])
  if(length(l)==0) character() else unlist(l)
}

