## Supplemental Figure 7D --- clonal expansion in spike-specific B cells
##
## Per-donor proportion of cells called clonally expanded.
##
## PROVENANCE. This panel is the prior analyst's analysis, not the step 2 one,
## and it is NOT reproducible from this repository. Her definition calls a cell
## expanded when it shares an identical V gene, J gene, C gene and CDR3 with
## another cell from the same donor, and her April 2026 run recomputed
## nearest-neighbour distances on a cell set that was never saved; the distance
## tables that are saved are from an earlier (August 2025) run and give
## p = 0.666, not the published 0.829.
##
## So the per-donor values are carried as an input, exactly as the sciCSR,
## V-QUEST and MSigDB steps elsewhere in this repo are. They currently come from
## `inputs/clonal_expansion_per_donor_prior_analyst.csv`, which was recovered
## from the vector coordinates of her published plot
## (Clonal_Expanded_proprotions_stripplot_new_threshold_April2026.pdf) and
## verified: the recovered values reproduce the published Wilcoxon p = 0.829 and
## the published group medians to four decimal places. **Replace that file with
## her own per-donor table when she can supply it** -- values read back off a
## figure are a stopgap, not a provenance.
##
## Output: figures/SupplFig7D_clonal_expansion.pdf/.png
##         figure_source_data/SupplFig7D_expansion_per_donor.csv
##         figure_source_data/SupplFig7D_statistics.csv
##
## Run:  Rscript figure_panels/SupplFig7D_clonal_expansion.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})

AS_PERCENT <- TRUE   # the submitted panel's axis is "% Expanded Cells"

## Cross-check against the published panel. Both must still hold after any
## change to the input file, or the panel has drifted from the figure.
PUBLISHED_P          <- 0.829
PUBLISHED_MEDIAN_TC  <- 0.0371
PUBLISHED_MEDIAN_DS  <- 0.0418

src <- file.path(REPO_ROOT, "inputs", "clonal_expansion_per_donor_prior_analyst.csv")
if (!file.exists(src)) stop("missing ", src)

s7d <- read_csv(src, show_col_types = FALSE) %>%
  dplyr::mutate(studyGroup = factor(studyGroup, levels = c("TC", "DS")))
stopifnot(nrow(s7d) == 16, all(table(s7d$studyGroup) == 8))

yv   <- if (AS_PERCENT) "pctExpanded" else "propExpanded"
ylab <- if (AS_PERCENT) "% Expanded Cells" else "Proportion of expanded cells"

x <- s7d$propExpanded[s7d$studyGroup == "TC"]
y <- s7d$propExpanded[s7d$studyGroup == "DS"]
pW <- suppressWarnings(wilcox.test(x, y)$p.value)

s7dStats <- tibble::tibble(
  nDonorCtrl = length(x), nDonorDS = length(y),
  medianCtrl = median(x), medianDS = median(y),
  meanCtrl = mean(x), meanDS = mean(y),
  pWilcoxon = pW)
print(s7dStats %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))) %>% as.data.frame())

ok <- function(a, b, tol) isTRUE(abs(a - b) < tol)
if (!ok(pW, PUBLISHED_P, 5e-4) ||
    !ok(median(x), PUBLISHED_MEDIAN_TC, 5e-4) || !ok(median(y), PUBLISHED_MEDIAN_DS, 5e-4))
  stop(sprintf(paste("this input no longer reproduces the published panel:",
                     "p = %.4f (expected %.3f), medians %.4f / %.4f (expected %.4f / %.4f)"),
               pW, PUBLISHED_P, median(x), median(y), PUBLISHED_MEDIAN_TC, PUBLISHED_MEDIAN_DS))
cat(sprintf("reproduces the published panel: p = %.4f, medians %.4f / %.4f\n",
            pW, median(x), median(y)))

fmtP <- function(p) format(signif(p, 3), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)

pS7D <- ggplot(s7d, aes(studyGroup, .data[[yv]], colour = studyGroup)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, colour = "black",
               fatten = 1, linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, size = 2.4, alpha = 0.9) +
  annotate("text", x = 1.5, y = max(s7d[[yv]]) * 1.08, size = 3,
           label = paste0("Wilcoxon p = ", fmtP(pW))) +
  scale_colour_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = GROUP_LABELS) +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = ylab,
       title = "Supplemental Figure 7D: clonal expansion",
       subtitle = "line = median | donor-level Wilcoxon rank-sum") +
  theme_classic(base_size = 11)

print(pS7D)
savePanel(pS7D, "SupplFig7D_clonal_expansion", width = 4.6, height = 4.2)

saveSource(
  s7d %>% dplyr::transmute(studyGroup = as.character(studyGroup),
                           groupLabel = GROUP_LABELS[as.character(studyGroup)],
                           donorRank, propExpanded, pctExpanded),
  "SupplFig7D_expansion_per_donor",
  "Per-donor proportion of clonally expanded cells (the plotted points)",
  panel = "Suppl 7D", role = "plotted")

saveSource(s7dStats, "SupplFig7D_statistics",
           "Donor-level Wilcoxon rank-sum on the per-donor proportions",
           panel = "Suppl 7D", role = "statistics")
