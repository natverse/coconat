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
#'   only used by the coconatfly package.
#' @param inherits The name of a previously registered dataset. See details.
#' @param ... Additional named arguments specifying properties of the dataset
#'
#' @return No return value. Called for its side effect.
#' @export
#' @details You can use the \code{inherits} argument to simplify adding an
#' additional handler for an existing dataset. For example imagine you have some
#' of your own annotations that you would like to supplement publicly released
#' ones for the banc dataset. You can register a new dataset \code{bancx} and
#' inherit from the \code{banc} definition and only replace the \code{metafun}
#' argument keeping everything else the same.
#'
#' @examples
#' \dontrun{
#' # partial example. metafun and partnerfun are pretty important to specify
#' register_dataset("flywire", shortname='fw', species='Drosophila melanogaster', sex='F', age="adult")
#' }
register_dataset <- function(name, shortname=NULL, species=NULL,
                             sex=c("F", "M", "H", "U"), age=NULL,
                             idfun=NULL, metafun=NULL, partnerfun=NULL,
                             namespace='default', inherits=NULL, ...) {

  ns=dataset_namespace(namespace)
  nn=ls(ns)
  if(name %in% nn) {
    warning("Dataset with name: ", name, " already registered. Overwriting!")
    ns[[name]]=NULL
  }

  baselist <- if(!is.null(inherits)) {
    if(!inherits %in% nn)
      stop("You have asked to inherit from a non-existent dataset:",
           inherits,
           "\nNB hard-coded coconatfly datasets do not yet support this mechanism.")
    dataset_details(inherits, namespace = namespace)
  } else NULL

  if(is.null(shortname))
    shortname=unname(abbreviate(name, minlength = 2))

  sex=match.arg(sex)

  if(is.null(idfun))
    idfun <- function(ids, integer64=FALSE) {default_id_fun(ids, metafun = metafun, integer64 = integer64)}

  mf <- match.call(expand.dots = FALSE)
  newlist <- list(
    name=name,
    shortname=shortname,
    species=species,
    sex=sex,
    age=age,
    idfun=idfun,
    metafun=metafun,
    partnerfun=partnerfun,
    call=mf,
    ...
  )
  if(!is.null(baselist)) {
    newlist=newlist[!sapply(newlist, is.null)]
    baselist[names(newlist)]=newlist
    newlist <- baselist

  }
  ns[[name]]=newlist

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
#' @param namespace Optional character vector specifying a namespace used to
#'   organise datasets (advanced use only).
#' @param return.short Whether to return the long or short name
#' @param match Whether the query should match against long or short forms of
#'   the dataset name.
#'
#' @return A character vector of names
#' @export
#'
#' @examples
#' dataset_names()
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

dataset_summaries <- function(namespace='default') {
  check_package_available("dplyr")
  dns=dataset_names(namespace = namespace)
  dd=sapply(dns, simplify = F, function(dn) {
    dd=dataset_details(dn, namespace = namespace)
    ndd=names(dd)[sapply(dd, is.atomic)]
    cns=names(dd$call)[-1]
    cvs=as.character(dd$call)[-1]
    l=as.list(cvs)
    names(l)=cns
    ddsel=c(dd[ndd], l[setdiff(names(l), ndd)])
    as.data.frame(ddsel)
  })
  dplyr::bind_rows(dd)
}
