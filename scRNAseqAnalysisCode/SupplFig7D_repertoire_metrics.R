## Supplemental Figure 7D --- donor-level repertoire metrics
##
## The six interpretable repertoire metrics compared between groups. These sit
## to the right of the clonal-expansion plot and share its panel letter: the
## legend reads "(D) Comparing multiple features of the repertoire to assess
## clonal expansion", so panel D is the whole row.
##
## Input:  pseudobulk_clean_output/step2_cells_chains.rds  ($repertoireMetrics)
## Output: figures/SupplFig7D_repertoire_metrics.pdf/.png
##         figure_source_data/SupplFig7D_repertoire_metrics.csv
##         figure_source_data/SupplFig7D_metric_statistics.csv
##
## Run:  Rscript figure_panels/SupplFig7D_repertoire_metrics.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
local({ for (p in c("figure_panels/_s7_common.R", "_s7_common.R")) if (file.exists(p)) return(source(p)) })

## Metrics reported, in the published left-to-right order. Depth-driven metrics
## (raw Shannon, Hill numbers, group counts) are computed in step 2 but not
## reported, because they correlate with cells sequenced at r ~ 1.
REPORT_METRICS <- c("rarefiedRichness", "rarefiedHill1", "evenness",
                    "pctCellsInExpandedGroups", "top1Share", "vGeneShannon")
METRIC_LABELS  <- c(rarefiedRichness         = "Richness\n(rarefied)",
                    rarefiedHill1            = "Hill\n(rarefied)",
                    evenness                 = "Evenness",
                    pctCellsInExpandedGroups = "%Cells in\nexpanded groups",
                    top1Share                = "top1Share",
                    vGeneShannon             = "Shannon")

s2 <- step2Object()
rep <- s2$repertoireMetrics %>%
  dplyr::mutate(studyGroup = factor(studyGroup, levels = c("TC", "DS")))
missing <- setdiff(REPORT_METRICS, names(rep))
if (length(missing)) stop("repertoireMetrics is missing: ", paste(missing, collapse = ", "))

s7e <- rep %>%
  dplyr::select(project, studyGroup, nCells, dplyr::all_of(REPORT_METRICS)) %>%
  tidyr::pivot_longer(dplyr::all_of(REPORT_METRICS), names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = REPORT_METRICS),
                metricLabel = factor(METRIC_LABELS[as.character(metric)],
                                     levels = METRIC_LABELS[REPORT_METRICS])) %>%
  dplyr::arrange(metric, studyGroup, as.numeric(sub("^P635-", "", project)))

## Donor-level Wilcoxon per metric, plus each metric's correlation with
## sequencing depth -- a metric that tracks depth is measuring library size.
s7eStats <- s7e %>%
  dplyr::group_split(metric) %>%
  lapply(function(d) {
    x <- d$value[d$studyGroup == "TC"]; y <- d$value[d$studyGroup == "DS"]
    tibble::tibble(metric = as.character(d$metric[1]),
                   metricLabel = gsub("\n", " ", as.character(d$metricLabel[1])),
                   nCtrl = length(x), nDS = length(y),
                   medianCtrl = median(x), medianDS = median(y),
                   pWilcoxon = suppressWarnings(wilcox.test(x, y)$p.value),
                   corWithDepth = cor(d$value, d$nCells, method = "spearman"))
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(padjBH = p.adjust(pWilcoxon, method = "BH"),
                metric = factor(metric, levels = REPORT_METRICS)) %>%
  dplyr::arrange(metric)

print(s7eStats %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
cat(sprintf("\nmetrics with BH-adjusted p < 0.05: %d of %d\n",
            sum(s7eStats$padjBH < 0.05), nrow(s7eStats)))

labs2 <- s7eStats %>%
  dplyr::mutate(metricLabel = factor(METRIC_LABELS[as.character(metric)],
                                     levels = METRIC_LABELS[REPORT_METRICS]),
                label = sprintf("p = %.2f", pWilcoxon))

pS7E <- ggplot(s7e, aes(studyGroup, value, colour = studyGroup)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black",
               fatten = 1, linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.85) +
  geom_text(data = labs2, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.4, size = 3) +
  facet_wrap(~ metricLabel, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = GROUP_LABELS) +
  labs(x = NULL, y = NULL,
       title = "Supplemental Figure 7D: repertoire metrics, DS vs Ctrl",
       subtitle = "line = median | donor-level Wilcoxon rank-sum") +
  theme_bw(base_size = 10)

print(pS7E)
savePanel(pS7E, "SupplFig7D_repertoire_metrics", width = 11, height = 3.4)

saveSource(
  s7e %>% dplyr::transmute(studyGroup = as.character(studyGroup),
                           groupLabel = GROUP_LABELS[as.character(studyGroup)],
                           project, nCellsSequenced = nCells,
                           metric = as.character(metric),
                           metricLabel = gsub("\n", " ", as.character(metricLabel)),
                           value),
  "SupplFig7D_repertoire_metrics",
  "Per-donor value of each reported repertoire metric (the plotted points)",
  panel = "Suppl 7D", role = "plotted")

saveSource(s7eStats %>% dplyr::mutate(metric = as.character(metric)),
           "SupplFig7D_metric_statistics",
           "Donor-level Wilcoxon per metric, with BH adjustment and correlation with sequencing depth",
           panel = "Suppl 7D", role = "statistics")
