## Shared helpers for the per-panel figure scripts.
##
## Each Fig*.R script in this directory is a thin layer on top of the step1--step4
## pipeline: the pipeline does the analysis and saves intermediates under
## pseudobulk_clean_output/, and each panel script reads one intermediate and
## emits (a) the panel and (b) its JCI Insight source-data table.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(rlang)
})

## Repo root, whether the script is run from the root or from figure_panels/.
findRepoRoot <- function() {
  for (p in c(".", "..")) {
    if (dir.exists(file.path(p, "pseudobulk_clean_output"))) return(normalizePath(p))
  }
  stop("Run this from the repo root (or figure_panels/); pseudobulk_clean_output/ not found.")
}

REPO_ROOT <- findRepoRoot()
OUT_DIR   <- file.path(REPO_ROOT, "pseudobulk_clean_output")
FIG_DIR   <- file.path(OUT_DIR, "figures")
SRC_DIR   <- file.path(OUT_DIR, "figure_source_data")
for (d in c(FIG_DIR, SRC_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## Published colours and the Ctrl/DS labelling used in the figures. "TC" is the
## internal study-group code; the manuscript calls that group "Ctrl".
GROUP_COLORS <- c(TC = "blue", DS = "red")
GROUP_LABELS <- c(TC = "Ctrl", DS = "DS")

savePanel <- function(plot, panel, width = 6.5, height = 4) {
  ggsave(file.path(FIG_DIR, sprintf("%s.pdf", panel)), plot, width = width, height = height)
  ggsave(file.path(FIG_DIR, sprintf("%s.png", panel)), plot, width = width, height = height,
         dpi = 300, bg = "white")
  cat(sprintf("panel:       %-42s %s\n", sprintf("%s.pdf/.png", panel), "figures/"))
  invisible(plot)
}

## Per-panel source tables, plus a generated manifest.
##
## Each table is written as a CSV (one source of truth per table, diffable in
## git) and registered in `_manifest.csv`. build_supporting_data_values.R reads
## that manifest to assemble the single Excel workbook JCI Insight requires --
## one tab per figure panel -- so the manifest is never hand-maintained.
##
##   panel  the figure panel this table belongs to, e.g. "Fig 5A". Becomes the
##          worksheet name, so keep it short (Excel allows 31 characters).
##   role   "plotted"    the values behind the graph itself (required by JCI)
##          "supporting" intermediate values a test is computed from
##          "statistics" the test results reported on the panel
MANIFEST   <- file.path(SRC_DIR, "_manifest.csv")
ROLE_ORDER <- c("plotted", "supporting", "statistics")

saveSource <- function(dat, name, description, panel = NA_character_,
                       role = c("plotted", "supporting", "statistics")) {
  role <- match.arg(role)
  path <- file.path(SRC_DIR, sprintf("%s.csv", name))
  write_csv(dat, path, na = "")

  row <- tibble::tibble(panel = as.character(panel), role = role,
                        file = basename(path), description = description,
                        nRows = nrow(dat), nCols = ncol(dat))
  man <- if (file.exists(MANIFEST)) {
    readr::read_csv(MANIFEST, show_col_types = FALSE, col_types = readr::cols(.default = "c")) %>%
      dplyr::mutate(nRows = as.integer(nRows), nCols = as.integer(nCols))
  } else row[0, ]
  man <- man %>%
    dplyr::filter(.data$file != row$file) %>%
    dplyr::bind_rows(row) %>%
    dplyr::arrange(.data$panel, match(.data$role, ROLE_ORDER), .data$file)
  write_csv(man, MANIFEST, na = "")

  cat(sprintf("source data: %-42s %4d rows  [%s/%s]  (%s)\n",
              basename(path), nrow(dat), panel, role, description))
  invisible(path)
}
