myAncestors <- function (id) 
{
  pagesize <- 20
  stopifnot(inherits(id, "Term"))
  url0 <- id@links$hierarchicalAncestors[[1]]
  if (is.null(url0)) {
    message("No ancestor terms.")
    return(NULL)
  }
  url <- paste0(url0, "?page=0&size=", pagesize)
  x <- httr::GET(url)
  httr::stop_for_status(x)
  cx <- httr::content(x)
  if (cx$page$totalElements > pagesize) {
    pagesize <- cx$page$totalElements
    url <- paste0(url0, "?page=0&size=", pagesize)
    x <- httr::GET(url)
    httr::warn_for_status(x)
    cx <- httr::content(x)
  }
  ans <- lapply(cx[["_embedded"]][[1]], rols:::makeTerm)
  names(ans) <- sapply(ans, rols::termId)
  rols:::Terms(x = ans)
}