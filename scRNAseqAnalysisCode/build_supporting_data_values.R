## Build the JCI Insight Supporting Data Values sheets.
##
## JCI Insight requires the values behind every graph and every reported mean,
## with a separate tab per figure panel, in one Excel file separate from other
## supplemental material (https://www.jci.org/kiosks/publish/supporting-data).
##
## On this manuscript the corresponding author keeps that master workbook, so
## the default here is ONE FILE PER PANEL, each containing a single worksheet
## already named as the tab should be named in the master ("Fig 5A"). Copying it
## across -- right-click the sheet tab > Move or Copy > choose the master -- then
## carries the tab name with it, and nothing has to be renamed by hand.
##
## Input:  pseudobulk_clean_output/figure_source_data/_manifest.csv  (+ the CSVs)
## Output: pseudobulk_clean_output/supporting_data_values/Fig5A.xlsx, Fig5B.xlsx, ...
##         and, if OUTPUT_MODE includes it, a single combined workbook.
##
## Run:  Rscript figure_panels/build_supporting_data_values.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
suppressPackageStartupMessages(library(openxlsx))

## "per_panel" (default), "combined", or "both".
OUTPUT_MODE <- "per_panel"

SDV_DIR      <- file.path(OUT_DIR, "supporting_data_values")
COMBINED_XLSX <- file.path(OUT_DIR, "Supporting_Data_Values.xlsx")

if (!file.exists(MANIFEST))
  stop("no manifest at ", MANIFEST, " -- run the figure_panels/Fig*.R scripts first")

man <- read_csv(MANIFEST, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(panel), nzchar(panel))
if (!nrow(man)) stop("manifest has no panel-assigned tables")

## Panels in figure order: figure number, then panel letter.
panelKey <- function(p) {
  fig <- suppressWarnings(as.numeric(sub("^\\D*(\\d+).*$", "\\1", p)))
  ltr <- sub("^.*?(\\d+)\\s*", "", p)
  sprintf("%03d-%s", ifelse(is.na(fig), 999, fig), ltr)
}
panels <- man %>% dplyr::distinct(panel) %>%
  dplyr::arrange(panelKey(panel)) %>% dplyr::pull(panel)

ROLE_HEADINGS <- c(plotted    = "Values plotted in the panel",
                   supporting = "Donor-level values behind the reported test",
                   statistics = "Statistics reported on the panel")

## Excel sheet names: <= 31 characters and none of []:*?/\
sheetName <- function(p) {
  ## A regex class starting "[:" is read as a POSIX class name, so replace the
  ## forbidden characters one at a time, literally.
  for (ch in c("[", "]", ":", "*", "?", "/", "\\")) p <- gsub(ch, "-", p, fixed = TRUE)
  substr(p, 1, 31)
}

## Write one panel's blocks onto a worksheet that already exists in `wb`.
writePanel <- function(wb, sheet, pn) {
  tabs <- man %>% dplyr::filter(panel == pn) %>%
    dplyr::arrange(match(role, ROLE_ORDER), file)

  writeData(wb, sheet, pn, startRow = 1, startCol = 1)
  addStyle(wb, sheet, createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
  row <- 3L

  for (i in seq_len(nrow(tabs))) {
    tb  <- tabs[i, ]
    dat <- read_csv(file.path(SRC_DIR, tb$file), show_col_types = FALSE, na = "")

    writeData(wb, sheet, sprintf("%s -- %s", ROLE_HEADINGS[[tb$role]], tb$description),
              startRow = row, startCol = 1)
    addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = row, cols = 1)
    row <- row + 1L

    writeData(wb, sheet, as.data.frame(dat), startRow = row, startCol = 1,
              headerStyle = createStyle(textDecoration = "bold", border = "bottom"),
              keepNA = FALSE)
    row <- row + nrow(dat) + 3L
  }

  setColWidths(wb, sheet, cols = 1:12, widths = "auto")
  freezePane(wb, sheet, firstActiveRow = 4)
  invisible(sum(tabs$nRows))
}

## --- one file per panel ------------------------------------------------------
if (OUTPUT_MODE %in% c("per_panel", "both")) {
  dir.create(SDV_DIR, showWarnings = FALSE, recursive = TRUE)
  for (pn in panels) {
    sheet <- sheetName(pn)
    wb <- createWorkbook()
    addWorksheet(wb, sheet)
    n <- writePanel(wb, sheet, pn)
    path <- file.path(SDV_DIR, sprintf("%s.xlsx", gsub("[^A-Za-z0-9]", "", pn)))
    saveWorkbook(wb, path, overwrite = TRUE)
    cat(sprintf("%-14s sheet \"%s\"  %4d data rows\n", basename(path), sheet, n))
  }
  cat(sprintf("\n%d file(s) in %s\n", length(panels),
              sub(paste0("^", REPO_ROOT, "/"), "", SDV_DIR)))
  cat("Copy into the master workbook: right-click the sheet tab > Move or Copy >\n",
      "select the master > tick \"Create a copy\". The tab keeps its name.\n", sep = "")
}

## --- single combined workbook ------------------------------------------------
if (OUTPUT_MODE %in% c("combined", "both")) {
  wb <- createWorkbook()
  for (pn in panels) { addWorksheet(wb, sheetName(pn)); writePanel(wb, sheetName(pn), pn) }
  saveWorkbook(wb, COMBINED_XLSX, overwrite = TRUE)
  cat(sprintf("\ncombined: %s (%d panels)\n", basename(COMBINED_XLSX), length(panels)))
} else if (file.exists(COMBINED_XLSX)) {
  ## Don't leave a stale combined workbook lying around next to the per-panel files.
  unlink(COMBINED_XLSX)
  cat("\nremoved stale ", basename(COMBINED_XLSX), "\n", sep = "")
}

cat("\nMethods callout required by JCI Insight:\n",
    "  \"Values for all data points in graphs are reported in the\n",
    "   Supporting Data Values file.\"\n", sep = "")
