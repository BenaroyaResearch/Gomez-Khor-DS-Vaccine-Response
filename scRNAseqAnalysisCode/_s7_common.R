## Shared loader for the Supplemental Figure 7 panels that come from step 2.
##
## Supplemental figure numbering changed from S5 to S7 between manuscript
## versions; these panels were called SupplFig5B/C/D in step 2's own output.
step2Object <- function() {
  p <- file.path(OUT_DIR, "step2_cells_chains.rds")
  if (!file.exists(p)) stop("missing ", p, " -- knit step 1 then step 2 first")
  s2 <- readRDS(p)
  if (is.null(s2$libAnno))
    stop("this step2_cells_chains.rds predates the libAnno addition; re-knit step 2")
  s2
}

## Donors ordered by group, then numerically within group (P635-2 before P635-10).
donorOrder <- function(project) order(as.numeric(sub("^P635-", "", project)))
donorLevelsFrom <- function(df) {
  d <- df %>% dplyr::distinct(project, studyGroup) %>%
    dplyr::mutate(num = as.numeric(sub("^P635-", "", project))) %>%
    dplyr::arrange(studyGroup, num)
  d$project
}
