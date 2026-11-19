#' @title Throw a retrogroup error
#' @noRd
#' @return Nothing. Throws an error of class `"retrogroup_error"`.
#' @param message Character vector of the error message, optionally with
#'   rlang bullet names such as `"x"` and `"i"`.
#' @param ... Named fields passed to [rlang::abort()].
#' @examples
#'   try(throw_error(c("Something went wrong.", i = "Try again.")))
throw_error <- function(message, ...) {
  rlang::abort(message, ..., class = "retrogroup_error", call = NULL)
}

#' @title Throw a retrogroup warning
#' @noRd
#' @return `NULL` (invisibly). Throws a warning of class
#'   `"retrogroup_warning"`.
#' @param message Character vector of the warning message, optionally with
#'   rlang bullet names such as `"x"` and `"i"`.
#' @param ... Named fields passed to [rlang::warn()].
#' @examples
#'   throw_warning("Something looks odd.")
throw_warning <- function(message, ...) {
  rlang::warn(message, ..., class = "retrogroup_warning")
}

#' @title Validate a character string argument
#' @noRd
#' @return `NULL` (invisibly). Throws an error if `x` is not a single
#'   non-missing, non-empty character string.
#' @param x Object to check.
#' @param name Character string, name of the argument for error messages.
#' @examples
#'   assert_string("time", "time")
#'   try(assert_string(1, "time"))
assert_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    throw_error(
      paste(name, "must be a single non-missing non-empty character string.")
    )
  }
  invisible()
}

#' @title Validate column name arguments
#' @noRd
#' @return `NULL` (invisibly). Throws an error if any column name is not a
#'   single string, two arguments name the same column, or a column is not
#'   in `data`.
#' @param data Data frame.
#' @param columns Named list of column name arguments. Names are the
#'   argument names, values are the column names. `NULL` values are skipped.
#' @examples
#'   data <- group_example_patient()
#'   assert_columns(data, list(arm = "arm", time = "time"))
#'   try(assert_columns(data, list(arm = "arm", time = "arm")))
assert_columns <- function(data, columns) {
  columns <- Filter(Negate(is.null), columns)
  for (name in names(columns)) {
    assert_string(columns[[name]], name)
  }
  values <- unlist(columns, use.names = FALSE)
  if (anyDuplicated(values)) {
    throw_error(
      paste(
        paste(names(columns), collapse = ", "),
        "must name different columns."
      )
    )
  }
  missing <- setdiff(values, colnames(data))
  if (length(missing)) {
    throw_error(
      paste0(
        "column(s) not in data: ",
        paste0("\"", missing, "\"", collapse = ", "),
        "."
      )
    )
  }
  invisible()
}

#' @title Disclosure note about adjusted patient counts
#' @noRd
#' @return Character string, a note reminding the user that the patient
#'   counts define the reconstructed subgroups, and that choosing counts
#'   different from the publication is their decision and must be
#'   disclosed.
#' @examples
#'   message_count_adjustment()
message_count_adjustment <- function() {
  paste(
    "The patient counts define the reconstructed subgroups",
    "and need not match the publication verbatim.",
    "Choosing counts that add up to the total number of patients",
    "is a decision about the reconstruction,",
    "and it is your responsibility.",
    "Please disclose it transparently wherever the reconstructed",
    "data or downstream results are reported."
  )
}
