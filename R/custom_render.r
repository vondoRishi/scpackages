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
        params = list(set_title = title), output_dir = getwd(), quiet = TRUE
    )
    print(paste("Saving data", paste(title, format(Sys.time(), '%d_%B_%Y'),"rds", sep = ".")))
    saveRDS(obj_scRNA,paste(title, format(Sys.time(), '%d_%B_%Y'),"rds", sep = "."))
}

#' Title
#'
#' @param obj_scRNA
#' @param title
#'
#' @return
#' @export
#'
#' @import DT
#'
#' @examples
report_Cluster <- function(obj_scRNA,title) {
    render(
        system.file("rmd_template", "CC_Cluster_report.Rmd", package = "scpackages"),
        output_file = paste(title, format(Sys.time(), '%d_%B_%Y'), sep = "_"),
        params = list(set_title = title), output_dir = getwd() , quiet = TRUE
    )
    print(paste("Saving data", paste(title, format(Sys.time(), '%d_%B_%Y'),"rds", sep = ".")))
    saveRDS(obj_scRNA,paste(title, format(Sys.time(), '%d_%B_%Y'),"rds", sep = "."))
}

