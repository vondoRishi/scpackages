#' Title
#'
#' @param df
#' @param digits
#'
#' @return
#' @export
#'
#' @examples
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}

tmpfun <- function(a,b,...) {
  print((match.call()))
  print("2nd")
  print(as.list(match.call()))
  print(as.list(match.call(expand.dots=FALSE)))
}
