# Render the computationally intensive vignettes.
#
# Each vignettes/*.Rmd.source file is knitted into a static vignettes/*.Rmd
# file with all results baked in. R CMD check and pkgdown only render the
# static *.Rmd files, so they never repeat the expensive computation.
#
# Workflow:
#   1. Commit and push the code changes to GitHub.
#   2. Install retrogroup from GitHub so packageDescription() records the
#      commit hash in the GithubSHA1 field:
#        pak::pak("wlandau/retrogroup")
#   3. Run render.R from any working directory, e.g.:
#        Rscript vignettes/render.R
#   4. Commit the generated *.Rmd files and *_figures/ directories.
#
# The vignettes render in parallel, one mirai daemon per vignette. Each
# daemon is a separate R process, so the vignettes do not share global
# variables, loaded packages, or random number streams. Each vignette
# writes to its own *.Rmd file and *_figures/ directory, so the processes
# never write to the same file.
render <- function(source, directory) {
  # knit() resolves fig.path and other relative paths against the working
  # directory, so knit from vignettes/ and restore the directory on exit.
  withr::local_dir(directory)
  output <- sub("\\.source$", "", source)
  # Stop on chunk errors instead of baking the error message into the output.
  knitr::opts_chunk$set(error = FALSE)
  knitr::knit(input = source, output = output, quiet = TRUE)
  # knit() evaluates inline R in the YAML header, so restore fields that
  # must stay relative to the vignettes/ directory.
  source_lines <- readLines(source)
  output_lines <- readLines(output)
  for (field in c("bibliography", "csl")) {
    pattern <- sprintf("^%s:", field)
    source_line <- grep(pattern, source_lines)[1L]
    output_line <- grep(pattern, output_lines)[1L]
    if (!is.na(source_line) && !is.na(output_line)) {
      output_lines[output_line] <- source_lines[source_line]
    }
  }
  writeLines(output_lines, output)
  output
}

# Find the package root by searching upward from this script, so the
# working directory does not matter. The script path comes from Rscript
# (--file=) or source() (ofile). If neither is available, search upward
# from the working directory instead.
script_directory <- function() {
  file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) {
      file <- frame$ofile
      break
    }
  }
  if (length(file) == 1L && nzchar(file)) dirname(normalizePath(file)) else "."
}
directory <- rprojroot::find_package_root_file(
  "vignettes",
  path = script_directory()
)
sources <- list.files(directory, pattern = "\\.Rmd\\.source$")
if (!length(sources)) {
  stop("No *.Rmd.source files found in ", directory, call. = FALSE)
}
mirai::daemons(length(sources))
message("Rendering in parallel: ", paste(sources, collapse = ", "))
tasks <- mirai::mirai_map(
  sources,
  render,
  .args = list(directory = directory)
)
# [.stop] waits for every task and stops with an error if any task failed.
outputs <- tasks[.stop]
mirai::daemons(0L)
message("Rendered: ", paste(unlist(outputs), collapse = ", "))
