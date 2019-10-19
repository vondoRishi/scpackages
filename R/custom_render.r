#' Title
#'
#' @param obj_scRNA
#' @param title
#'
#' @return
#' @export
#'
#' @import rmarkdown
#' @examples
report_QC <- function(obj_scRNA,title) {
      render(
        system.file("rmd_template", "QC_report.Rmd", package = "scpackages"),
        output_file = paste(title, format(Sys.time(), '%d_%B_%Y'), sep = "_"),
        params = list(set_title = title), output_dir = getwd()
    )
    saveRDS(obj_scRNA,paste(title, format(Sys.time(), '%d_%B_%Y'),"rds", sep = "."))
}

