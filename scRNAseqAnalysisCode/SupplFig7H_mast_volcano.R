## Supplemental Figure 7H --- DS effect, MAST with donor as a random effect
##
## Volcano of the pooled DS effect from the combined MAST model, which fits both
## isotypes together with donor as a random effect. Chromosome 21 genes are
## distinguished because trisomy 21 makes them the expected signal.
##
## Input:  pseudobulk_clean_output/figure_source_data/mast_combined/
##           Step4_mastComb_additive_DS_effect.csv  (written by
##           P635_step4_DE_mast_combined.Rmd; knit that first)
## Output: figures/SupplFig7H_mast_volcano.pdf/.png
##         figure_source_data/SupplFig7H_mast_volcano.csv
##         figure_source_data/SupplFig7H_called_genes.csv
##
## Run:  Rscript figure_panels/SupplFig7H_mast_volcano.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
suppressPackageStartupMessages(library(ggrepel))

## Calling thresholds, as in the step 4 model notebook.
FDR     <- 0.10
LFC_MIN <- 1

MAST_SRC <- file.path(SRC_DIR, "mast_combined", "Step4_mastComb_additive_DS_effect.csv")
if (!file.exists(MAST_SRC))
  stop("missing ", MAST_SRC, " -- knit P635_step4_DE_mast_combined.Rmd first")

addDS <- read_csv(MAST_SRC, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(padj), !is.na(log2FC)) %>%
  dplyr::mutate(significant = padj < FDR & abs(log2FC) > LFC_MIN,
                chr21 = !is.na(chr) & chr == "21")

cat(sprintf("%d genes on the volcano; %d called (%d of them chr21)\n",
            nrow(addDS), sum(addDS$significant), sum(addDS$significant & addDS$chr21)))

## One categorical variable drives both shape and fill, so the two scales merge
## into a single four-key legend rather than two the reader must combine.
LV <- c("DE, chr21", "DE, non-chr21", "not DE, chr21", "not DE, non-chr21")
addDS <- addDS %>% dplyr::mutate(grp = factor(dplyr::case_when(
   significant &  chr21 ~ LV[1],
   significant & !chr21 ~ LV[2],
  !significant &  chr21 ~ LV[3],
  TRUE                  ~ LV[4]), levels = LV))

## Called genes, plus the single smallest adjusted p-value among the uncalled.
## That gene (SON) misses the fold-change cutoff but is the highest point on the
## plot, so leaving it anonymous would mislead. Picked by rank, not by name, so
## it survives a re-run that shifts the numbers.
labV <- dplyr::bind_rows(
  addDS %>% dplyr::filter(significant),
  addDS %>% dplyr::filter(!significant) %>% dplyr::slice_min(padj, n = 1, with_ties = FALSE)
) %>% dplyr::distinct()
cat(sprintf("labelled: %s\n", paste(sort(labV$gene), collapse = ", ")))

xl <- max(abs(addDS$log2FC))
keyGuide <- guide_legend(override.aes = list(size = 3, colour = "white", stroke = 0.4))

pS7I <- ggplot(addDS, aes(log2FC, -log10(padj))) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_hline(yintercept = -log10(FDR), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(-LFC_MIN, LFC_MIN), linetype = "dashed", colour = "grey60") +
  geom_point(data = ~ dplyr::filter(.x, grp == LV[4]),
             aes(shape = grp, fill = grp), colour = "grey85", size = 1.2, stroke = 0, alpha = 0.55) +
  geom_point(data = ~ dplyr::filter(.x, grp == LV[3]),
             aes(shape = grp, fill = grp), colour = "white", size = 2.7, stroke = 0.35) +
  geom_point(data = ~ dplyr::filter(.x, significant),
             aes(shape = grp, fill = grp), colour = "white", size = 3.6, stroke = 0.5) +
  geom_text_repel(data = labV, aes(label = gene, fontface = ifelse(chr21, "bold", "plain")),
                  size = 2.9, colour = "grey15", max.overlaps = Inf, box.padding = 0.45,
                  segment.size = 0.2, segment.colour = "grey60", min.segment.length = 0) +
  scale_shape_manual(name = NULL, breaks = LV, limits = LV, drop = FALSE,
                     values = setNames(c(22, 21, 22, 21), LV), guide = keyGuide) +
  scale_fill_manual(name = NULL, breaks = LV, limits = LV, drop = FALSE,
                    values = setNames(c("#B2182B", "#B2182B", "#2166ac", "grey85"), LV),
                    guide = keyGuide) +
  coord_cartesian(xlim = c(-xl, xl)) +
  labs(x = expression(log[2] * " fold change, DS vs Ctrl"),
       y = expression(-log[10] * "(FDR)"),
       title = "Supplemental Figure 7H: DS effect, MAST with donor random effect",
       subtitle = sprintf("called at FDR < %.2f and |log2FC| > %g", FDR, LFC_MIN)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", legend.text = element_text(size = 8.5))

print(pS7I)
savePanel(pS7I, "SupplFig7H_mast_volcano", width = 7, height = 5.6)

saveSource(
  addDS %>% dplyr::transmute(gene, log2FC, ci.lo, ci.hi, p, padj,
                             chr, chr21, significant, category = as.character(grp)),
  "SupplFig7H_mast_volcano",
  "Pooled DS effect per gene from the combined MAST model (the plotted points)",
  panel = "Suppl 7H", role = "plotted")

saveSource(
  labV %>% dplyr::transmute(gene, log2FC, p, padj, chr, chr21, significant) %>%
    dplyr::arrange(padj),
  "SupplFig7H_called_genes",
  "The labelled genes: every called gene, plus the smallest adjusted p-value among the uncalled",
  panel = "Suppl 7H", role = "supporting")
