#' Add counts to the levels of a factor
#'
#' @param f a factor variable
#'
#' @export
add_counts <- function(f) {
    levels(f) <- paste0(levels(f), " (", table(f), ")")
    f
}


#' perform an operation on each line of a rowwise data frame
#'
#' A \code{\link{do}} operation given multiple arguments creates a rowwise
#' data frame. Often one wants to follow this with an operation on each line
#' of the data frame, but this requires regrouping (which leads to an error),
#' then performing an awkward expression with lots of \code{.$column[[1]]}
#' in it. This is a simple wrapper of \code{do} that makes that easier.
#'
#' @param x a tbl
#' @param ... Additional arguments passed on to \code{do}
#' @param groups Vector of columns to group by. If missing,
#' will attempt to guess based on all columns that are not lists.
#'
#' @name do_row
#'
#' @examples
#'
#' library(dplyr)
#' library(broom)
#'
#' g <- mtcars %>% group_by(cyl) %>% do(m = lm(mpg ~ wt, .))
#' g %>% do_row(tidy(m))
#' g %>% do_row(augment(m))
#' g %>% do_row(tidied = tidy(m), augmented = augment(m))
#'
#' @import dplyr
#'
#' @export
do_row <- function(x, ..., groups = NULL) {
    if (is.null(groups)) {
        # guess how to group: group on all items that aren't lists
        classes <- sapply(x, class)
        groups <- names(classes)[classes != "list"]
    }
    # group x, suppressing warning about "strips rowwise nature"
    g <- suppressWarnings(group_by_(x, .dots = groups))
    args <- dplyr:::dots(...)
    named <- dplyr:::named_args(args)
    subargs <- lapply(args, function(a) substitute((with(lapply(., function(z) z[[1]]), e)), list(e = a)))

    #print(g)
    ret <- do.call(do, c(list(g), subargs))
    #print(ret)

    if (named) {
        # add columns to original data
        x <- x %>% inner_join(ret, by = groups)
        return(x)
    }

    ret
}



#' Expand a dataset to include all factorial combinations of one or more
#' variables
#'
#' @param .data a tbl
#' @param ... arguments
#' @param stringsAsFactors logical specifying if character vectors are
#' converted to factors.
#'
#' @return A data.frame, grouped by the arguments in ...
#'
#' @import dplyr
#'
#' @export
inflate <- function(.data, ..., stringsAsFactors = FALSE) {
    ret <- expand.grid(..., stringsAsFactors = stringsAsFactors)
    ret %>% group_by_(colnames(ret)) %>% do(.data)
}
