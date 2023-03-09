custom_interactive_heatmap <- function(hm) {
  check_package_available('InteractiveComplexHeatmap', repo = 'Bio')
  check_package_available('kableExtra')
  check_package_available('shiny')
  ht = ComplexHeatmap::draw(hm)

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
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, ht,
                                                             click_action = click_action, brush_action = brush_action)
  }

  shiny::shinyApp(ui, server)
}

# private function to prepare a cosine heatmap using either
#' @importFrom stats heatmap as.dist hclust
#' @importFrom grDevices hcl.colors
prepare_cosine_heatmap <- function(x, labRow, interactive=FALSE,
                                   col=hcl.colors(12, "YlOrRd", rev = TRUE),
                                   method=c("ward.D", "single", "complete", "average",
                                            "mcquitty", "median", "centroid", "ward.D2"),
                                   ...) {
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
  } else heatmap(x,
               distfun = function(x) as.dist(1-x),
               hclustfun = function(...) hclust(..., method=method),
               symm = T, keep.dendro = T, labRow=labRow, col=col, ...)
}
