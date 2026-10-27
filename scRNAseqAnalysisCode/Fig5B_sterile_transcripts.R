## Main text Figure 5B --- sterile transcripts at the IgM and IgG loci
##
## "Frequency of sterile transcripts at IgM and IgG loci in spike-specific IgM+
## B cells from participants with and without DS."
##
## Input:  inputs/sterile_transcript_annot_all_donors_Jan2026.csv
##         (per-cell sterile read counts; see inputs/README.md for provenance.
##          Produced upstream by P635_Sterile_tx_all_donors_Jan2026.Rmd running
##          sciCSR over ~80 GB of per-donor BAMs, which are not in this repo --
##          this cached table is the reproducible starting point, and is what the
##          prior analyst's own plotting code read.)
## Output: pseudobulk_clean_output/figures/Fig5B_sterile_transcripts.pdf/.png
##         pseudobulk_clean_output/figure_source_data/Fig5B_sterile_transcripts.csv
##         pseudobulk_clean_output/figure_source_data/Fig5B_donor_medians.csv
##         pseudobulk_clean_output/figure_source_data/Fig5B_statistics.csv
##
## Run:  Rscript figure_panels/Fig5B_sterile_transcripts.R

## Locate _shared.R whether this is run from the repo root or from figure_panels/.
local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})

## ---------------------------------------------------------------------------
## Parameters
## ---------------------------------------------------------------------------

## Cells plotted: the productive isotype of the cell. The panel is about sterile
## (germline) transcription in cells that have NOT yet switched, so it is drawn
## on IgM-expressing cells.
FIG5B_CELL_ISOTYPE <- "IGHM"

## Sterile loci shown, in the published left-to-right order. The quantification
## also covers IGHA1-S, IGHA2-S and IGHE-S; those are not in the panel.
FIG5B_LOCI   <- c("IGHM-S", "IGHG1-S", "IGHG2-S", "IGHG3-S", "IGHG4-S")
LOCUS_LABELS <- c(`IGHM-S` = "IgM", `IGHG1-S` = "IgG1", `IGHG2-S` = "IgG2",
                  `IGHG3-S` = "IgG3", `IGHG4-S` = "IgG4")

## The test is a donor-level Wilcoxon rank-sum on each donor's MEDIAN per-cell
## sterile read count (exact = FALSE), as published. Cells are not independent
## within a donor, so the donor is the unit of analysis even though every cell is
## plotted.
##
## UNTESTABLE_LABEL: at a locus where every donor's median is 0 in both groups
## the rank-sum statistic has zero variance and wilcox.test() returns NaN. The
## submitted panel prints that through as "p=NaN" (IgG3). "nan" reproduces the
## submitted panel; "explicit" prints "not testable" instead.
##   "nan"      prints "p=NaN", exactly as the submitted panel does
##   "na"       prints "p=NA" -- same meaning, without the R-internal spelling
##   "explicit" prints "not testable"
## The statistics table always records this as a missing value (NA), never NaN
## and never a fabricated 1, with a `testable` flag alongside it.
UNTESTABLE_LABEL <- "nan"

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------

sterile <- read_csv(file.path(REPO_ROOT, "inputs",
                              "sterile_transcript_annot_all_donors_Jan2026.csv"),
                    show_col_types = FALSE, na = "")

## `donorId` in this table is the P635-n project code, which the rest of the
## pipeline calls `project`; renamed here so the source tables agree.
fig5b <- sterile %>%
  dplyr::rename(project = donorId) %>%
  dplyr::filter(c_gene_family == FIG5B_CELL_ISOTYPE, sterile_type %in% FIG5B_LOCI) %>%
  dplyr::mutate(
    studyGroup = factor(studyGroup, levels = c("TC", "DS")),
    locus      = factor(sterile_type, levels = FIG5B_LOCI),
    locusLabel = factor(LOCUS_LABELS[as.character(sterile_type)],
                        levels = LOCUS_LABELS[FIG5B_LOCI])) %>%
  dplyr::arrange(locus, studyGroup, as.numeric(sub("^P635-", "", project)), libid)

## Every plotted cell must appear once per locus, and both groups must be present.
stopifnot(dplyr::n_distinct(fig5b$libid) * length(FIG5B_LOCI) == nrow(fig5b))
stopifnot(all(c("TC", "DS") %in% levels(droplevels(fig5b$studyGroup))))
cat(sprintf("%d %s+ cells from %d donors x %d loci = %d rows\n",
            dplyr::n_distinct(fig5b$libid), FIG5B_CELL_ISOTYPE,
            dplyr::n_distinct(fig5b$project), length(FIG5B_LOCI), nrow(fig5b)))

## ---------------------------------------------------------------------------
## Statistics --- donor-level Wilcoxon on per-donor medians
## ---------------------------------------------------------------------------

donorMedians <- fig5b %>%
  dplyr::group_by(locus, locusLabel, studyGroup, project) %>%
  dplyr::summarise(nCells = dplyr::n(),
                   donorMedian = median(sterile_reads),
                   donorMean   = mean(sterile_reads), .groups = "drop")

fig5bStats <- donorMedians %>%
  dplyr::group_split(locus) %>%
  lapply(function(d) {
    x <- d$donorMedian[d$studyGroup == "TC"]
    y <- d$donorMedian[d$studyGroup == "DS"]
    p <- suppressWarnings(tryCatch(
      wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_))
    tibble::tibble(
      locus = as.character(d$locus[1]), locusLabel = as.character(d$locusLabel[1]),
      nDonorCtrl = length(x), nDonorDS = length(y),
      medianOfDonorMediansCtrl = median(x), medianOfDonorMediansDS = median(y),
      pctCellsNonZeroCtrl = NA_real_, pctCellsNonZeroDS = NA_real_,
      pDonorWilcoxon = p,
      testable = is.finite(p))
  }) %>%
  dplyr::bind_rows()

## Share of cells with any sterile read, per group --- this is what makes a locus
## untestable at the donor level even when individual cells are non-zero.
nonZero <- fig5b %>%
  dplyr::group_by(locus, studyGroup) %>%
  dplyr::summarise(pct = 100 * mean(sterile_reads > 0), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = studyGroup, values_from = pct)
fig5bStats <- fig5bStats %>%
  dplyr::select(-pctCellsNonZeroCtrl, -pctCellsNonZeroDS) %>%
  dplyr::left_join(nonZero %>% dplyr::rename(pctCellsNonZeroCtrl = TC,
                                             pctCellsNonZeroDS   = DS),
                   by = "locus") %>%
  dplyr::mutate(locus = factor(locus, levels = FIG5B_LOCI),
                locusLabel = factor(locusLabel, levels = LOCUS_LABELS[FIG5B_LOCI])) %>%
  dplyr::arrange(locus)

print(fig5bStats %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))) %>%
        as.data.frame())

## ---------------------------------------------------------------------------
## Panel
## ---------------------------------------------------------------------------

fmtP <- function(p) format(signif(p, 3), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
untestableLabel <- switch(UNTESTABLE_LABEL,
                          nan = "p=NaN", na = "p=NA", explicit = "not testable",
                          stop("UNTESTABLE_LABEL must be \"nan\", \"na\" or \"explicit\""))
fig5bStats$panelLabel <- ifelse(fig5bStats$testable,
                                paste0("p=", fmtP(fig5bStats$pDonorWilcoxon)),
                                untestableLabel)

## A locus with no defined statistic is a missing value in the table, not NaN.
fig5bStats$pDonorWilcoxon[!fig5bStats$testable] <- NA_real_

## Annotation sits above the tallest cell at that locus; the scales are free, so
## the headroom is computed per locus.
annot <- fig5b %>%
  dplyr::group_by(locus, locusLabel) %>%
  dplyr::summarise(yTop = max(sterile_reads), .groups = "drop") %>%
  dplyr::left_join(fig5bStats %>% dplyr::select(locus, panelLabel), by = "locus")

pFig5B <- ggplot(fig5b, aes(studyGroup, sterile_reads, fill = studyGroup)) +
  geom_violin(trim = FALSE, width = 0.5, alpha = 0.7, linewidth = 0.3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, fatten = 0,
               colour = "white", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 0.7, colour = "black") +
  geom_text(data = annot, inherit.aes = FALSE, vjust = 0, size = 3,
            aes(x = 1.5, y = yTop * 1.35, label = panelLabel)) +
  facet_wrap(~ locusLabel, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = GROUP_LABELS) +
  ## log1p is heavily compressed at the top, so the default break algorithm
  ## stacks labels on top of each other. Anchor at 0 and take a few log-spaced
  ## breaks above it, per locus (the scales are free).
  scale_y_continuous(
    trans  = "log1p",
    breaks = function(lims) {
      hi <- max(lims, na.rm = TRUE)
      if (!is.finite(hi) || hi <= 1) return(c(0, 1))
      unique(c(0, scales::breaks_log(4)(c(1, hi))))
    },
    labels = scales::label_number(big.mark = ",", accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "Sterile reads per cell (log1p)",
       title = "Figure 5B: sterile transcription at the IgM and IgG loci",
       subtitle = sprintf("%s+ cells | violin and points = cells, white line = median | donor-level Wilcoxon on per-donor medians",
                          LOCUS_LABELS[paste0(FIG5B_CELL_ISOTYPE, "-S")])) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank())

print(pFig5B)
savePanel(pFig5B, "Fig5B_sterile_transcripts", width = 9, height = 4.2)

## ---------------------------------------------------------------------------
## JCI Insight source data
## ---------------------------------------------------------------------------

saveSource(
  fig5b %>% dplyr::transmute(
    studyGroup = as.character(studyGroup),
    groupLabel = GROUP_LABELS[as.character(studyGroup)],
    project, libid,
    cellIsotype = c_gene_family,
    locus = as.character(locus), locusLabel = as.character(locusLabel),
    sterileReads = sterile_reads),
  "Fig5B_sterile_transcripts",
  "Sterile reads per cell at each locus in IgM+ cells (the plotted points)",
  panel = "Fig 5B", role = "plotted")

saveSource(
  donorMedians %>% dplyr::transmute(
    studyGroup = as.character(studyGroup),
    groupLabel = GROUP_LABELS[as.character(studyGroup)],
    project, locus = as.character(locus), locusLabel = as.character(locusLabel),
    nCells, donorMedian, donorMean) %>%
    dplyr::arrange(locus, studyGroup, as.numeric(sub("^P635-", "", project))),
  "Fig5B_donor_medians",
  "Per-donor median and mean sterile reads per cell (the donor-level values the test is computed from)",
  panel = "Fig 5B", role = "supporting")

saveSource(fig5bStats %>% dplyr::select(-panelLabel) %>%
             dplyr::mutate(locus = as.character(locus), locusLabel = as.character(locusLabel)),
           "Fig5B_statistics",
           "Donor-level Wilcoxon rank-sum per locus (n = 8 donors per group)",
           panel = "Fig 5B", role = "statistics")
