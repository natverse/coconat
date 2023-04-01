custom_interactive_heatmap <- function(hm) {
  check_package_available('InteractiveComplexHeatmap', repo = 'Bio')
  check_package_available('kableExtra')
  check_package_available('shiny')

  ui = shiny::fluidPage(
    InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(output_ui = shiny::htmlOutput("info")),
  )

  click_action = function(df, output) {
    output[["info"]] = shiny::renderUI({
      if(!is.null(df)) {
        shiny::HTML(glue::glue("<p style='background-color:#FF8080;color:white;padding:5px;'>You have clicked on heatmap {df$heatmap}, row {df$row_index}, column {df$column_index}</p>"))
      }
    })
  }
  brush_action = function(df, output) {
    row_index = unique(unlist(df$row_index))
    column_index = unique(unlist(df$column_index))
    output[["info"]] = shiny::renderUI({
      if(!is.null(df)) {
        shiny::HTML(kableExtra::kable_styling(kableExtra::kbl(hm@matrix[row_index, column_index, drop = FALSE], digits = 2, format = "html"), full_width = FALSE, position = "left"))
      }
    })
  }

  server = function(input, output, session) {
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, hm,
                                                             click_action = click_action, brush_action = brush_action)
  }

  shiny::shinyApp(ui, server)
}

# private function to draw a cosine heatmap using either the basic stats::heatmap
# or InteractiveComplexHeatmap
#
#' @importFrom stats heatmap as.dist hclust
#' @importFrom grDevices hcl.colors
cosine_heatmap <- function(x, labRow=rownames(x), interactive=FALSE,
                           heatmap=TRUE, col=hcl.colors(12, "YlOrRd", rev = TRUE),
                           method=c("ward.D", "single", "complete", "average",
                                    "mcquitty", "median", "centroid", "ward.D2"),
                                   ...) {
  method=match.arg(method)
  FUN=stats::heatmap
  if(!is.logical(heatmap)) {
    FUN=try(match.fun(heatmap))
    if(inherits(FUN, 'try-error'))
      stop("The heatmap argument should either be T/F or define a function in a way that match.fun() understands")
    heatmap=TRUE
  }

  if(interactive) {
    check_package_available('ComplexHeatmap', repo = 'Bio')
    hm <- ComplexHeatmap::Heatmap(
      x,
      row_labels=labRow,
      col=col,
      cluster_rows=function(x,...) hclust(as.dist(1-x), method=method,...),
      cluster_columns=function(x,...) hclust(as.dist(1-x), method=method,...),
      ...
    )
    custom_interactive_heatmap(hm)
  } else if(isTRUE(heatmap)) {
    FUN(x,
        distfun = function(x) as.dist(1-x),
        hclustfun = function(...) hclust(..., method=method),
        symm = T, keep.dendro = T, labRow=labRow, col=col, ...)
  } else {
    hclust(as.dist(1-x), method = method, ...)
  }
}
