## Main text Figure 5G --- directional concordance of DS effects between isotypes
##
## "For each isotype, directional concordance of each DEG's dysregulation in the
## other isotype was assessed (G)."
##
## Input:  inputs/mast_de/DEGs_IGH{G,M}_MAST.csv
##         Per-cell Seurat/MAST DE, DS vs control, run separately in IgG+ and
##         IgM+ cells. See inputs/README.md -- note this is the DEGs_*_MAST.csv
##         run, NOT the similarly named DEGs_*_MAST_timebinary_2026.csv.
## Output: pseudobulk_clean_output/figures/Fig5G_isotype_concordance.pdf/.png
##         pseudobulk_clean_output/figure_source_data/Fig5G_isotype_concordance.csv
##         pseudobulk_clean_output/figure_source_data/Fig5G_statistics.csv
##
## Run:  Rscript figure_panels/Fig5G_isotype_concordance.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
suppressPackageStartupMessages(library(ggrepel))

## ---------------------------------------------------------------------------
## Parameters
## ---------------------------------------------------------------------------

## A gene is shown if it is differentially expressed in IgG+ cells, IgM+ cells,
## or both; the marker says which. Same list that feeds the Figure 5H enrichment.
DEG_P_COLUMN    <- "p_val_adj"
DEG_P_THRESHOLD <- 0.05

## Genes labelled on the published panel.
LABEL_GENES <- c("IFNGR2", "SAMSN1")

GROUP_LEVELS <- c("IgG+", "IgM+", "Both")
GROUP_SHAPES <- c(`IgG+` = 25, `IgM+` = 22, `Both` = 24)   # down tri, square, up tri
GROUP_FILLS  <- c(`IgG+` = "black", `IgM+` = "red", `Both` = "forestgreen")

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------

## Gene symbols are in an unnamed first column, so read_csv names it "...1".
readDE <- function(isotype)
  read_csv(file.path(REPO_ROOT, "inputs", "mast_de", sprintf("DEGs_%s_MAST.csv", isotype)),
           show_col_types = FALSE, name_repair = "minimal") %>%
    dplyr::rename(gene = 1)

deG <- readDE("IGHG"); deM <- readDE("IGHM")
sig <- function(d) d$gene[!is.na(d[[DEG_P_COLUMN]]) & d[[DEG_P_COLUMN]] < DEG_P_THRESHOLD]

fig5g <- tibble::tibble(gene = sort(union(sig(deG), sig(deM)))) %>%
  dplyr::left_join(deG %>% dplyr::select(gene, log2FC_IgG = avg_log2FC,
                                         p_IgG = p_val, padj_IgG = p_val_adj), by = "gene") %>%
  dplyr::left_join(deM %>% dplyr::select(gene, log2FC_IgM = avg_log2FC,
                                         p_IgM = p_val, padj_IgM = p_val_adj), by = "gene") %>%
  dplyr::mutate(
    significantIgG = !is.na(padj_IgG) & padj_IgG < DEG_P_THRESHOLD,
    significantIgM = !is.na(padj_IgM) & padj_IgM < DEG_P_THRESHOLD,
    significantIn  = factor(dplyr::case_when(significantIgG & significantIgM ~ "Both",
                                             significantIgG ~ "IgG+", TRUE ~ "IgM+"),
                            levels = GROUP_LEVELS),
    ## A gene tested in only one isotype has no second fold change, so it cannot
    ## be placed on the plane. Kept in the source table, with the reason.
    plotted        = !is.na(log2FC_IgG) & !is.na(log2FC_IgM),
    concordantSign = dplyr::if_else(plotted, sign(log2FC_IgG) == sign(log2FC_IgM), NA))

plotted <- fig5g %>% dplyr::filter(plotted)
cat(sprintf("%d genes DE in either isotype; %d have a fold change in both and are plotted\n",
            nrow(fig5g), nrow(plotted)))
print(plotted %>% dplyr::count(significantIn, name = "nPlotted") %>% as.data.frame())
notPlotted <- fig5g %>% dplyr::filter(!plotted)
if (nrow(notPlotted))
  cat(sprintf("not plotted (tested in one isotype only): %s\n",
              paste(notPlotted$gene, collapse = ", ")))

## ---------------------------------------------------------------------------
## Statistics
## ---------------------------------------------------------------------------

fig5gStats <- tibble::tibble(
  nGenesDEEither   = nrow(fig5g),
  nPlotted         = nrow(plotted),
  nDEInIgGOnly     = sum(plotted$significantIn == "IgG+"),
  nDEInIgMOnly     = sum(plotted$significantIn == "IgM+"),
  nDEInBoth        = sum(plotted$significantIn == "Both"),
  nConcordantSign  = sum(plotted$concordantSign),
  propConcordant   = mean(plotted$concordantSign),
  ## Are the fold changes of these genes more often same-signed than chance?
  binomialP        = binom.test(sum(plotted$concordantSign), nrow(plotted), 0.5)$p.value,
  pearsonPlotted   = cor(plotted$log2FC_IgG, plotted$log2FC_IgM),
  spearmanPlotted  = cor(plotted$log2FC_IgG, plotted$log2FC_IgM, method = "spearman"))
print(fig5gStats %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))) %>%
        as.data.frame())

## ---------------------------------------------------------------------------
## Panel
## ---------------------------------------------------------------------------

lim <- max(abs(c(plotted$log2FC_IgG, plotted$log2FC_IgM))) * 1.08

pFig5G <- ggplot(plotted, aes(log2FC_IgM, log2FC_IgG)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(shape = significantIn, fill = significantIn),
             colour = "black", size = 2.6, stroke = 0.3) +
  ggrepel::geom_text_repel(
    data = plotted %>% dplyr::filter(gene %in% LABEL_GENES),
    aes(label = gene), size = 3, min.segment.length = 0, segment.size = 0.3,
    box.padding = 0.8, nudge_y = 1.1, seed = 1) +
  scale_shape_manual(values = GROUP_SHAPES, drop = FALSE, name = "Significantly different in:") +
  scale_fill_manual(values = GROUP_FILLS, drop = FALSE, name = "Significantly different in:") +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = expression(log[2] ~ "FC IgM (DS vs Ctrl)"),
       y = expression(log[2] ~ "FC IgG (DS vs Ctrl)"),
       title = "Figure 5G: directional concordance between isotypes",
       subtitle = sprintf("%d genes DE in IgG+ cells, IgM+ cells or both (MAST, adjusted p < %g)",
                          nrow(plotted), DEG_P_THRESHOLD)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

print(pFig5G)
savePanel(pFig5G, "Fig5G_isotype_concordance", width = 6, height = 6.4)

## ---------------------------------------------------------------------------
## JCI Insight source data
## ---------------------------------------------------------------------------

saveSource(
  fig5g %>% dplyr::transmute(
    gene, significantIn = as.character(significantIn), plotted, concordantSign,
    log2FC_IgG, p_IgG, padj_IgG, log2FC_IgM, p_IgM, padj_IgM),
  "Fig5G_isotype_concordance",
  "Genes differentially expressed in IgG+ cells, IgM+ cells or both, with the fold change in each (the plotted points)",
  panel = "Fig 5G", role = "plotted")

saveSource(fig5gStats, "Fig5G_statistics",
           "Counts by isotype and sign concordance across the plotted genes",
           panel = "Fig 5G", role = "statistics")
