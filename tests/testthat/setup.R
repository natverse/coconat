# Route implicit graphics output to an in-memory/null PDF device during tests.
# Without this, heatmap tests can create tests/testthat/Rplots.pdf.
options(device = function(...) grDevices::pdf(file = NULL))
