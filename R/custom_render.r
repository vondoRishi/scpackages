report_QC <- function(obj_scRNA,title) {
  library(rmarkdown)
  render("QC_report.Rmd",output_file = paste(title,format(Sys.time(), '%d_%B_%Y'),sep = "_"), params = list(set_title=title))
}

