## Main text Figure 5D --- somatic hypermutation at each Ig locus
##
## "Quantitating somatic hypermutation (mutations per cell) at each Ig locus in
## spike-specific IgG+ (above) and IgM+ (below) B cells."
##
## Input:  inputs/imgt_vquest/{IgH,IgK,IgL}_{IgG,IgM}.txt
##         IMGT/V-QUEST "8_V-REGION-nt-mutation-statistics.txt" exports, one per
##         locus x productive isotype. See inputs/README.md for provenance.
## Output: pseudobulk_clean_output/figures/Fig5D_somatic_hypermutation.pdf/.png
##         pseudobulk_clean_output/figure_source_data/Fig5D_somatic_hypermutation.csv
##         pseudobulk_clean_output/figure_source_data/Fig5D_donor_medians.csv
##         pseudobulk_clean_output/figure_source_data/Fig5D_statistics.csv
##
## Run:  Rscript figure_panels/Fig5D_somatic_hypermutation.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
suppressPackageStartupMessages({ library(stringr); library(purrr) })

## ---------------------------------------------------------------------------
## Parameters
## ---------------------------------------------------------------------------

FIG5D_LOCI     <- c("IgH", "IgK", "IgL")   # columns, left to right
FIG5D_ISOTYPES <- c("IgG", "IgM")          # rows: IgG+ above, IgM+ below
ISOTYPE_LABELS <- c(IgG = "IgG+ B cells", IgM = "IgM+ B cells")

## V-QUEST reports several mutation counts; the panel plots the total.
METRIC       <- "Total_mutations"
METRIC_LABEL <- "Total Mutations"

## The test is a two-sided Welch t-test on each donor's MEDIAN mutations per
## sequence. Donor, not sequence, is the unit of analysis: sequences from one
## donor are not independent.
##
## UNTESTABLE_LABEL: where every donor's median is 0 in both groups (IgM+ at the
## IgH locus -- IgM+ cells are largely unmutated, which is the expected biology)
## the t statistic is undefined. The submitted panel prints "p = NaN".
##   "nan" (as published) | "na" | "explicit"
UNTESTABLE_LABEL <- "nan"

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------

## V-QUEST writes counts like "12 (3.4%)"; take the leading integer, and treat a
## missing or unparseable entry as zero mutations, as the original analysis did.
extractCount <- function(x) tidyr::replace_na(as.numeric(str_extract(as.character(x), "^\\d+")), 0)

readLocus <- function(locus, isotype) {
  path <- file.path(REPO_ROOT, "inputs", "imgt_vquest", sprintf("%s_%s.txt", locus, isotype))
  if (!file.exists(path)) stop("missing V-QUEST export: ", path)
  read_delim(path, delim = "\t", show_col_types = FALSE, name_repair = "minimal") %>%
    dplyr::transmute(
      locus       = locus,
      cellIsotype = isotype,
      sequenceId  = `Sequence ID`,
      ## Study group and donor are encoded in the FASTA header V-QUEST echoes
      ## back, e.g. "..._P635-12_TC". Cross-checked against step 1 below.
      studyGroup  = dplyr::if_else(str_detect(`Sequence ID`, "TC$"), "TC", "DS"),
      project     = str_extract(`Sequence ID`, "P635-\\d+"),
      !!METRIC   := extractCount(.data[["V-REGION Nb of mutations"]]))
}

fig5d <- tidyr::expand_grid(locus = FIG5D_LOCI, isotype = FIG5D_ISOTYPES) %>%
  purrr::pmap_dfr(function(locus, isotype) readLocus(locus, isotype)) %>%
  dplyr::mutate(
    studyGroup  = factor(studyGroup, levels = c("TC", "DS")),
    locus       = factor(locus, levels = FIG5D_LOCI),
    cellIsotype = factor(cellIsotype, levels = FIG5D_ISOTYPES),
    isotypeLabel = factor(ISOTYPE_LABELS[as.character(cellIsotype)],
                          levels = ISOTYPE_LABELS[FIG5D_ISOTYPES])) %>%
  dplyr::arrange(cellIsotype, locus, studyGroup,
                 as.numeric(sub("^P635-", "", project)))

stopifnot(!any(is.na(fig5d$project)))

## The group label here is parsed from a sequence-ID suffix, so check it against
## the donor table step 1 wrote rather than trusting the string.
step1 <- read_csv(file.path(OUT_DIR, "step1_library_qc_metrics.csv"), show_col_types = FALSE) %>%
  dplyr::distinct(project, studyGroup) %>%
  dplyr::rename(studyGroupStep1 = studyGroup)
check <- fig5d %>% dplyr::distinct(project, studyGroup) %>%
  dplyr::left_join(step1, by = "project")
if (any(is.na(check$studyGroupStep1)) ||
    any(as.character(check$studyGroup) != check$studyGroupStep1)) {
  print(check %>% dplyr::filter(is.na(studyGroupStep1) |
                                as.character(studyGroup) != studyGroupStep1))
  stop("V-QUEST sequence-ID study groups disagree with step 1")
}
cat(sprintf("%d sequences across %d locus x isotype panels; group labels agree with step 1\n",
            nrow(fig5d), length(FIG5D_LOCI) * length(FIG5D_ISOTYPES)))

## ---------------------------------------------------------------------------
## Statistics --- Welch t-test on per-donor medians
## ---------------------------------------------------------------------------

donorMedians <- fig5d %>%
  dplyr::group_by(locus, cellIsotype, isotypeLabel, studyGroup, project) %>%
  dplyr::summarise(nSequences  = dplyr::n(),
                   donorMedian = median(.data[[METRIC]]),
                   donorMean   = mean(.data[[METRIC]]), .groups = "drop")

fig5dStats <- donorMedians %>%
  dplyr::group_split(cellIsotype, locus) %>%
  lapply(function(d) {
    x <- d$donorMedian[d$studyGroup == "TC"]
    y <- d$donorMedian[d$studyGroup == "DS"]
    p <- suppressWarnings(tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_))
    tibble::tibble(
      cellIsotype = as.character(d$cellIsotype[1]), locus = as.character(d$locus[1]),
      nDonorCtrl = length(x), nDonorDS = length(y),
      nSeqCtrl = sum(d$nSequences[d$studyGroup == "TC"]),
      nSeqDS   = sum(d$nSequences[d$studyGroup == "DS"]),
      meanOfDonorMediansCtrl = mean(x), meanOfDonorMediansDS = mean(y),
      pWelchTOnDonorMedians = p, testable = is.finite(p))
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(locus = factor(locus, levels = FIG5D_LOCI),
                cellIsotype = factor(cellIsotype, levels = FIG5D_ISOTYPES)) %>%
  dplyr::arrange(cellIsotype, locus)

print(fig5dStats %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))) %>%
        as.data.frame())

## ---------------------------------------------------------------------------
## Panel
## ---------------------------------------------------------------------------

fmtP <- function(p) format(signif(p, 3), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
untestable <- switch(UNTESTABLE_LABEL, nan = "p = NaN", na = "p = NA",
                     explicit = "not testable",
                     stop("UNTESTABLE_LABEL must be \"nan\", \"na\" or \"explicit\""))
fig5dStats$panelLabel <- ifelse(fig5dStats$testable,
                                paste0("p = ", fmtP(fig5dStats$pWelchTOnDonorMedians)),
                                untestable)
fig5dStats$pWelchTOnDonorMedians[!fig5dStats$testable] <- NA_real_

annot <- fig5d %>%
  dplyr::group_by(locus, cellIsotype, isotypeLabel) %>%
  dplyr::summarise(yTop = max(.data[[METRIC]]), .groups = "drop") %>%
  dplyr::left_join(fig5dStats %>% dplyr::select(locus, cellIsotype, panelLabel),
                   by = c("locus", "cellIsotype"))

pFig5D <- ggplot(fig5d, aes(studyGroup, .data[[METRIC]], colour = studyGroup)) +
  ## Points are sequences; the bar and whiskers are the mean and SD across them,
  ## as published. The reported test is the donor-level one above.
  geom_jitter(width = 0.25, height = 0, size = 0.8, alpha = 0.8) +
  stat_summary(fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.25, colour = "black", linewidth = 0.4) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, fatten = 0,
               colour = "black", linewidth = 0.6) +
  geom_text(data = annot, inherit.aes = FALSE, vjust = 0, size = 3,
            aes(x = 1.5, y = yTop * 1.08, label = panelLabel)) +
  facet_grid(isotypeLabel ~ locus, scales = "free_y", switch = "y") +
  scale_colour_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = GROUP_LABELS) +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = METRIC_LABEL,
       title = "Figure 5D: somatic hypermutation at each Ig locus",
       subtitle = "Points = sequences, bar = mean +/- SD | Welch t-test on per-donor medians") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.placement = "outside")

print(pFig5D)
savePanel(pFig5D, "Fig5D_somatic_hypermutation", width = 8, height = 6)

## ---------------------------------------------------------------------------
## JCI Insight source data
## ---------------------------------------------------------------------------

saveSource(
  fig5d %>% dplyr::transmute(
    studyGroup = as.character(studyGroup),
    groupLabel = GROUP_LABELS[as.character(studyGroup)],
    project, sequenceId,
    cellIsotype = as.character(cellIsotype), locus = as.character(locus),
    totalMutations = .data[[METRIC]]),
  "Fig5D_somatic_hypermutation",
  "Total V-region mutations per sequence at each locus (the plotted points)",
  panel = "Fig 5D", role = "plotted")

saveSource(
  donorMedians %>% dplyr::transmute(
    studyGroup = as.character(studyGroup),
    groupLabel = GROUP_LABELS[as.character(studyGroup)],
    project, cellIsotype = as.character(cellIsotype), locus = as.character(locus),
    nSequences, donorMedian, donorMean) %>%
    dplyr::arrange(cellIsotype, locus, studyGroup, as.numeric(sub("^P635-", "", project))),
  "Fig5D_donor_medians",
  "Per-donor median and mean mutations per sequence (the donor-level values the test is computed from)",
  panel = "Fig 5D", role = "supporting")

saveSource(
  fig5dStats %>% dplyr::select(-panelLabel) %>%
    dplyr::mutate(locus = as.character(locus), cellIsotype = as.character(cellIsotype)),
  "Fig5D_statistics",
  "Welch t-test on per-donor medians, per locus and productive isotype",
  panel = "Fig 5D", role = "statistics")
