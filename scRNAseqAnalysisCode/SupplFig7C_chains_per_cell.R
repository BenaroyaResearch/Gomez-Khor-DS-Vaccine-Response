## Supplemental Figure 7C --- heavy and light chains recovered per cell
##
## "Number of heavy and light chain sequences recovered per cell."
##
## Drawn upstream of the doublet cut on purpose: the cells carrying two or more
## heavy chains are the ones the doublet rule removes, and the point of the
## panel is to show that distribution.
##
## Input:  pseudobulk_clean_output/step2_cells_chains.rds  ($libAnno)
## Output: figures/SupplFig7C_chains_per_cell.pdf/.png
##         figure_source_data/SupplFig7C_chains_per_cell.csv
##         figure_source_data/SupplFig7C_chain_count_crosstab.csv
##
## Run:  Rscript figure_panels/SupplFig7C_chains_per_cell.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
local({ for (p in c("figure_panels/_s7_common.R", "_s7_common.R")) if (file.exists(p)) return(source(p)) })

LIBRARY_SET <- "qc_passing"

s2 <- step2Object()
libs <- if (LIBRARY_SET == "post_doublet") dplyr::filter(s2$libAnno, !is_doublet) else s2$libAnno

s7c <- libs %>%
  dplyr::transmute(libid, donorId, project,
                   studyGroup = factor(studyGroup, levels = c("TC", "DS")),
                   numheavy, numlight, is_doublet)

cat(sprintf("%d cells; %d flagged as doublets (>=2 heavy or >=3 light)\n",
            nrow(s7c), sum(s7c$is_doublet)))

## The jitter cannot convey how many cells sit on each integer pair, so the
## exact counts go out as a source table beside the plotted points.
crossTab <- s7c %>%
  dplyr::count(studyGroup, numheavy, numlight, name = "nCells") %>%
  dplyr::arrange(studyGroup, numheavy, numlight)
print(crossTab %>%
        tidyr::pivot_wider(names_from = numlight, values_from = nCells,
                           values_fill = 0, names_prefix = "light_") %>%
        as.data.frame())

pS7C <- ggplot(s7c, aes(numlight, numheavy, colour = studyGroup)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.6, size = 0.9) +
  facet_wrap(~ studyGroup, labeller = labeller(studyGroup = GROUP_LABELS)) +
  scale_colour_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_continuous(breaks = 0:max(s7c$numlight)) +
  scale_y_continuous(breaks = 0:max(s7c$numheavy)) +
  labs(x = "# Light Chains", y = "# Heavy Chains",
       title = "Supplemental Figure 7C: chains recovered per cell",
       subtitle = sprintf("%s (n = %d cells)", LIBRARY_SET, nrow(s7c))) +
  theme_bw(base_size = 11)

print(pS7C)
savePanel(pS7C, "SupplFig7C_chains_per_cell", width = 7, height = 3.8)

saveSource(
  s7c %>% dplyr::transmute(studyGroup = as.character(studyGroup),
                           groupLabel = GROUP_LABELS[as.character(studyGroup)],
                           project, donorId, libid, numheavy, numlight, is_doublet),
  "SupplFig7C_chains_per_cell",
  "Heavy and light chain counts for every cell (the plotted points)",
  panel = "Suppl 7C", role = "plotted")

saveSource(crossTab, "SupplFig7C_chain_count_crosstab",
           "Cells per heavy x light chain combination (the counts the jitter hides)",
           panel = "Suppl 7C", role = "supporting")
