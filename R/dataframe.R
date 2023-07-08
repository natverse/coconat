#' Add dataframe columns detailing ordering and groups from a dendrogram
#'
#' @param df A dataframe e.g. as returned by \code{coconatfly::cf_meta}
#' @param dend A \code{\link{dendrogram}} or \code{\link{hclust}} object
#' @param h A height to cut the dendrogram
#' @param k A number of clusters to cut the dendrogram
#' @param colnames The names of the two new columns
#' @param idcol The name of the column containing id information. If you do not
#'   provide this argument then the function will choose the first of three
#'   default columns present in \code{df} (with a warning when >1 column
#'   present).
#'
#' @return A copy of \code{df} with one or two extra columns called
#'   \code{dendid} and e.g. \code{group_h2}.
#' @export
add_cluster_info <- function(df, dend, h=NULL, k=NULL, colnames=NULL,
                             idcol=c('key','id', 'root_id', 'bodyid')) {
  stopifnot(inherits(dend, 'hclust') || inherits(dend, 'dendrogram'))
  stopifnot(inherits(df, 'data.frame'))
  if(missing(idcol)) {
    idcol=intersect(idcol, colnames(df))
    if(length(idcol)==0)
      stop("None of the standard id columns are present in ",
           deparse(substitute(df)))
    if(length(idcol)>1)
      warning("Multiple standard id columns are present in ",
           deparse(substitute(df)),
           "\nChoosing ", idcol <- idcol[1])
  }
  stopifnot(idcol %in% colnames(df))
  drop_group_col <- is.null(k) && is.null(h)
  if(drop_group_col) k=1

  check_package_available('dendroextras')
  check_package_available('dplyr')
  gg=dendroextras::slice(dend, h=h, k=k)

  # some juggling as we need ids to be character for matching purposes but we
  # can keep them e.g. as int64 if that's how they arrive
  saved_ids=df[[idcol]]
  df[[idcol]]=as.character(df[[idcol]])
  ggdf=data.frame(id=names(gg), dendid=seq_len(length(gg)), group=unname(gg))

  if(is.null(colnames))
    colnames="dendid"
  if(length(colnames)==1) {
    gname=if(is.null(k)) glue::glue('group_h{h}') else glue::glue('group_k{k}')
    colnames=c(colnames, gname)
  } else checkmate::check_character(colnames, len = 2, unique=TRUE)
  colnames(ggdf)=c(idcol, colnames)
  if(drop_group_col)
    ggdf=ggdf[-3L]

  res=dplyr::left_join(df, ggdf, by=idcol)
  res[[idcol]]=saved_ids
  res
}
