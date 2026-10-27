## Main text Figure 5A --- IgH isotype frequencies per donor
##
## "Frequency of each IgH isotype in BCRseq analyses of spike-specific B cells
## from participants with and without DS" (n = 8 per cohort).
##
## Input:  pseudobulk_clean_output/step2_cell_metadata.csv
##         (one row per analysis cell, written by P635_pseudobulk_step2_BCR_repertoire.Rmd;
##          run step 1 then step 2 first)
## Output: pseudobulk_clean_output/figures/Fig5A_isotype_frequencies.pdf/.png
##         pseudobulk_clean_output/figure_source_data/Fig5A_isotype_frequencies.csv
##         pseudobulk_clean_output/figure_source_data/Fig5A_statistics.csv
##
## Run:  Rscript figure_panels/Fig5A_isotype_frequencies.R

## Locate _shared.R whether this is run from the repo root or from figure_panels/.
local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})

## ---------------------------------------------------------------------------
## Parameters
## ---------------------------------------------------------------------------

## Isotypes plotted, in the published left-to-right order.
FIG5A_ISOTYPES <- c("IGHD", "IGHM", "IGHG", "IGHA")
ISOTYPE_LABELS <- c(IGHD = "IgD", IGHM = "IgM", IGHG = "IgG", IGHA = "IgA")

## Which summary line is drawn, and which test is annotated on the panel.
## Both tests and both summaries are always written to the statistics table, so
## switching these changes only what is drawn, never what is reported.
##
##   CENTER_LINE  "mean" reproduces the submitted panel; "median" matches the
##                house style used by the other figure legends in this paper.
##   TEST         "welch_t" reproduces the submitted panel's annotations;
##                "wilcoxon" matches the Figure 5 legend, which says "(A) Wilcoxon".
##   P_ADJUST     "none" reproduces the submitted panel. "BH" corrects across the
##                four isotypes of this panel; they are compositional (they sum to
##                100% within a donor), so the four tests are positively dependent
##                and BH remains valid.
CENTER_LINE <- "mean"
TEST        <- "wilcoxon"
P_ADJUST    <- "none"

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------

cells <- read_csv(file.path(OUT_DIR, "step2_cell_metadata.csv"), show_col_types = FALSE)

stopifnot(!any(cells$is_doublet))

## Denominator: the cells in that donor carrying an assignable IgH isotype. Cells
## with no heavy chain, or a heavy chain with no parseable C gene, are excluded
## from numerator and denominator alike. The denominator is kept in the source
## table so every percentage can be recomputed from counts.
fig5aDenominator <- cells %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::count(studyGroup, project, donorId, name = "nCellsWithIsotype")

fig5a <- cells %>%
  dplyr::filter(!is.na(isotype), isotype %in% FIG5A_ISOTYPES) %>%
  dplyr::count(studyGroup, project, donorId, isotype, name = "nCells") %>%
  ## every donor x isotype combination, so a donor with none of an isotype is 0
  tidyr::complete(tidyr::nesting(studyGroup, project, donorId),
                  isotype = FIG5A_ISOTYPES, fill = list(nCells = 0)) %>%
  dplyr::left_join(fig5aDenominator, by = c("studyGroup", "project", "donorId")) %>%
  dplyr::mutate(
    pctIsotype   = 100 * nCells / nCellsWithIsotype,
    studyGroup   = factor(studyGroup, levels = c("TC", "DS")),
    isotype      = factor(isotype, levels = FIG5A_ISOTYPES),
    isotypeLabel = factor(ISOTYPE_LABELS[as.character(isotype)],
                          levels = ISOTYPE_LABELS[FIG5A_ISOTYPES])) %>%
  ## Donors ordered by group, then numerically within group (P635-2 before P635-10).
  dplyr::arrange(isotype, studyGroup,
                 as.numeric(sub("^P635-", "", project)))

## The four plotted isotypes are the only isotypes called in these data, so each
## donor's four percentages must sum to exactly 100. This guards the per-donor
## normalization: the submitted panel divided donor P635-1's counts by 51 rather
## than its own denominator of 53, so its four percentages summed to 103.9%.
donorSums <- fig5a %>%
  dplyr::group_by(studyGroup, project, donorId) %>%
  dplyr::summarise(pctSum = sum(pctIsotype), .groups = "drop")
if (any(abs(donorSums$pctSum - 100) > 1e-6)) {
  print(donorSums %>% dplyr::filter(abs(pctSum - 100) > 1e-6))
  stop("per-donor isotype percentages do not sum to 100")
}
cat(sprintf("per-donor percentages sum to 100 for all %d donors\n", nrow(donorSums)))

## ---------------------------------------------------------------------------
## Statistics
## ---------------------------------------------------------------------------

## Both tests are computed and reported. At n = 8 per group the Wilcoxon makes no
## distributional assumption; the Welch t-test is what the submitted panel used.
fig5aStats <- fig5a %>%
  dplyr::group_split(isotype) %>%
  lapply(function(d) {
    x <- d$pctIsotype[d$studyGroup == "TC"]
    y <- d$pctIsotype[d$studyGroup == "DS"]
    tibble::tibble(
      isotype      = as.character(d$isotype[1]),
      isotypeLabel = as.character(d$isotypeLabel[1]),
      nCtrl = length(x), nDS = length(y),
      meanCtrl = mean(x),   meanDS = mean(y),
      medianCtrl = median(x), medianDS = median(y),
      pWelchT   = t.test(x, y)$p.value,
      pWilcoxon = suppressWarnings(wilcox.test(x, y)$p.value))
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    ## BH across the four isotypes of this panel, for both tests, regardless of
    ## which one is drawn. The four percentages are compositional (they sum to
    ## 100% within a donor) so the tests are positively dependent, under which
    ## BH remains valid.
    padjWelchT   = p.adjust(pWelchT, method = "BH"),
    padjWilcoxon = p.adjust(pWilcoxon, method = "BH"),
    testReported = TEST,
    pReported    = if (TEST == "wilcoxon") pWilcoxon else pWelchT,
    padjMethod   = P_ADJUST,
    pAdjusted    = if (P_ADJUST == "none") NA_real_
                   else if (TEST == "wilcoxon") padjWilcoxon else padjWelchT,
    isotypeLabel = factor(isotypeLabel, levels = ISOTYPE_LABELS[FIG5A_ISOTYPES]))

print(fig5aStats %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))) %>%
        as.data.frame())

## ---------------------------------------------------------------------------
## Panel
## ---------------------------------------------------------------------------

centerFun <- if (CENTER_LINE == "median") median else mean

## Three significant figures, never scientific notation, no trailing zeros ---
## so 0.586, 0.000142, 0.714 rather than 0.586000 or 1.42e-04.
fmtP <- function(p) format(signif(p, 3), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
fig5aStats$panelLabel <- if (P_ADJUST == "none") {
  paste0("p=", fmtP(fig5aStats$pReported))
} else {
  paste0("q=", fmtP(fig5aStats$pAdjusted))
}

## Headroom for the annotation so it never sits on the data.
yTop <- max(fig5a$pctIsotype)

pFig5A <- ggplot(fig5a, aes(studyGroup, pctIsotype, colour = studyGroup)) +
  stat_summary(fun = centerFun, geom = "crossbar", width = 0.55, colour = "black",
               fatten = 1, linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.9) +
  geom_text(data = fig5aStats, inherit.aes = FALSE, vjust = 0, size = 3,
            aes(x = 1.5, y = yTop * 1.06, label = panelLabel)) +
  facet_wrap(~ isotypeLabel, nrow = 1, strip.position = "bottom") +
  scale_colour_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = GROUP_LABELS) +
  scale_y_continuous(breaks = seq(0, 80, by = 20),   # as published
                     limits = c(0, yTop * 1.2),
                     expand = ggplot2::expansion(mult = c(0.02, 0))) +
  labs(x = NULL, y = "% Isotype+ per donor",
       title = "Figure 5A: IgH isotype frequencies",
       subtitle = sprintf("%% of isotype-assignable cells per donor | line = %s | %s%s",
                          CENTER_LINE,
                          if (TEST == "wilcoxon") "Wilcoxon rank-sum" else "Welch t-test",
                          if (P_ADJUST == "none") "" else sprintf(", %s-adjusted across 4 isotypes", P_ADJUST))) +
  theme_classic(base_size = 11) +
  theme(strip.placement = "outside", strip.background = element_blank())

print(pFig5A)
savePanel(pFig5A, "Fig5A_isotype_frequencies", width = 7, height = 4.2)

## ---------------------------------------------------------------------------
## JCI Insight source data
## ---------------------------------------------------------------------------

saveSource(
  fig5a %>%
    dplyr::transmute(
      studyGroup = as.character(studyGroup),
      groupLabel = GROUP_LABELS[as.character(studyGroup)],
      project, donorId,
      isotype = as.character(isotype),
      isotypeLabel = as.character(isotypeLabel),
      nCells, nCellsWithIsotype, pctIsotype),
  "Fig5A_isotype_frequencies", "Percentage of each donor's isotype-assignable cells carrying each IgH isotype (the plotted points)",
  panel = "Fig 5A", role = "plotted")

saveSource(fig5aStats %>% dplyr::select(-panelLabel),
           "Fig5A_statistics",
           "Group summaries and two-sided tests per isotype (n = 8 donors per group)",
           panel = "Fig 5A", role = "statistics")
