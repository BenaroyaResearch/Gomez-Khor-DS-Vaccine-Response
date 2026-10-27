## Supplemental Figure 7I --- Hallmark gene sets, cameraPR
##
## Competitive gene-set test on the pooled DS effect from the combined MAST
## model. cameraPR corrects for inter-gene correlation, which an
## over-representation test on a thresholded gene list does not.
##
## Input:  pseudobulk_clean_output/figure_source_data/mast_combined/
##           Step4_mastComb_GSEA_cameraPR.csv  (written by
##           P635_step4_DE_mast_combined.Rmd; knit that first)
## Output: figures/SupplFig7I_cameraPR_genesets.pdf/.png
##         figure_source_data/SupplFig7I_cameraPR_genesets.csv
##
## Run:  Rscript figure_panels/SupplFig7I_cameraPR_genesets.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})

GSEA_FDR    <- 0.10   # calling threshold, as in the step 4 model notebook
GSEA_N_SHOW <- 10     # how many sets to draw, by smallest cameraPR p

CAM_SRC <- file.path(SRC_DIR, "mast_combined", "Step4_mastComb_GSEA_cameraPR.csv")
if (!file.exists(CAM_SRC))
  stop("missing ", CAM_SRC, " -- knit P635_step4_DE_mast_combined.Rmd first")

cam <- read_csv(CAM_SRC, show_col_types = FALSE)
cat(sprintf("%d gene sets x %d effect(s): %s\n", dplyr::n_distinct(cam$pathwayShort),
            dplyr::n_distinct(cam$effect), paste(unique(cam$effect), collapse = "; ")))

camSig <- cam %>% dplyr::filter(cameraFDR < GSEA_FDR)
cat(sprintf("sets at cameraPR FDR < %.2f: %d%s\n", GSEA_FDR, nrow(camSig),
            if (nrow(camSig)) paste0(" (", paste(camSig$pathwayShort, collapse = ", "), ")") else ""))

## The drawn sets: the smallest cameraPR p-values, plus anything significant
## that would otherwise fall outside that cut.
highlight <- union(
  cam %>% dplyr::group_by(pathwayShort) %>%
    dplyr::summarise(bestP = min(cameraP), .groups = "drop") %>%
    dplyr::slice_min(bestP, n = GSEA_N_SHOW, with_ties = FALSE) %>% dplyr::pull(pathwayShort),
  camSig$pathwayShort) %>% unique()

## The axis is the FDR, not the raw p, so bar length and bar colour encode the
## same quantity: a bar crosses the dashed line exactly when it is coloured as
## significant.
s7j <- cam %>%
  dplyr::filter(pathwayShort %in% highlight) %>%
  dplyr::mutate(signedLogFDR = ifelse(direction == "Up", 1, -1) * -log10(cameraFDR),
                significant  = cameraFDR < GSEA_FDR) %>%
  dplyr::arrange(dplyr::desc(signedLogFDR))

print(s7j %>% dplyr::transmute(pathwayShort, nGenesInSet, direction,
                               cameraP = signif(cameraP, 3),
                               cameraFDR = signif(cameraFDR, 3),
                               signedLogFDR = round(signedLogFDR, 3)) %>% as.data.frame())

pS7J <- ggplot(s7j, aes(signedLogFDR, reorder(pathwayShort, signedLogFDR),
                        fill = significant)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_vline(xintercept = c(-1, 1) * -log10(GSEA_FDR), linetype = "dashed", colour = "grey55") +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey75"),
                    labels = c(`TRUE` = sprintf("cameraPR FDR < %.2f", GSEA_FDR),
                               `FALSE` = "n.s."), name = NULL) +
  labs(x = bquote("signed " * -log[10] * "(cameraPR FDR)   (positive: up in DS; dashed: FDR " * .(GSEA_FDR) * ")"),
       y = NULL,
       title = "Supplemental Figure 7I: Hallmark gene sets, cameraPR",
       subtitle = "pooled DS effect, combined MAST model with donor random effect") +
  theme_classic(base_size = 11) +
  theme(legend.position = if (nrow(camSig)) "right" else "none")

print(pS7J)
savePanel(pS7J, "SupplFig7I_cameraPR_genesets", width = 7.5, height = 5)

saveSource(
  s7j %>% dplyr::transmute(effect, pathwayShort, nGenesInSet, direction,
                           cameraP, cameraFDR, signedLogFDR, significant),
  "SupplFig7I_cameraPR_genesets",
  "cameraPR result for each plotted Hallmark gene set (the plotted bars)",
  panel = "Suppl 7I", role = "plotted")

saveSource(
  cam %>% dplyr::transmute(effect, pathwayShort, nGenesInSet, direction, cameraP, cameraFDR) %>%
    dplyr::arrange(cameraP),
  "SupplFig7I_cameraPR_all_sets",
  "cameraPR result for all tested Hallmark sets, not only the ten drawn",
  panel = "Suppl 7I", role = "supporting")
