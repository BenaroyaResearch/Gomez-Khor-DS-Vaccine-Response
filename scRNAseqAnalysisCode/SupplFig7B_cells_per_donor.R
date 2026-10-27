## Supplemental Figure 7B --- cells sequenced per donor
##
## "Number of cells successfully sequenced per donor."
##
## Input:  pseudobulk_clean_output/step2_cells_chains.rds  ($libAnno)
## Output: figures/SupplFig7B_cells_per_donor.pdf/.png
##         figure_source_data/SupplFig7B_cells_per_donor.csv
##
## Run:  Rscript figure_panels/SupplFig7B_cells_per_donor.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
local({ for (p in c("figure_panels/_s7_common.R", "_s7_common.R")) if (file.exists(p)) return(source(p)) })

## Which libraries the panel counts. "qc_passing" is upstream of the doublet
## cut, matching panel C, which needs the doublets visible.
LIBRARY_SET <- "qc_passing"   # "qc_passing" or "post_doublet"

s2 <- step2Object()
libs <- if (LIBRARY_SET == "post_doublet") dplyr::filter(s2$libAnno, !is_doublet) else s2$libAnno

s7b <- libs %>%
  dplyr::count(studyGroup, project, donorId, name = "nCells") %>%
  dplyr::mutate(studyGroup = factor(studyGroup, levels = c("TC", "DS")),
                project = factor(project, levels = donorLevelsFrom(.))) %>%
  dplyr::arrange(project)

cat(sprintf("%d donors, %d libraries (%s)\n", nrow(s7b), sum(s7b$nCells), LIBRARY_SET))

pS7B <- ggplot(s7b, aes(project, nCells, fill = studyGroup)) +
  geom_col() +
  scale_fill_manual(values = GROUP_COLORS, labels = GROUP_LABELS, name = NULL) +
  labs(x = "Donor", y = "# Cells",
       title = "Supplemental Figure 7B: cells sequenced per donor",
       subtitle = sprintf("%s (n = %d libraries)", LIBRARY_SET, sum(s7b$nCells))) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(pS7B)
savePanel(pS7B, "SupplFig7B_cells_per_donor", width = 6.5, height = 4)

saveSource(
  s7b %>% dplyr::transmute(studyGroup = as.character(studyGroup),
                           groupLabel = GROUP_LABELS[as.character(studyGroup)],
                           project = as.character(project), donorId, nCells),
  "SupplFig7B_cells_per_donor",
  "Number of QC-passing libraries (cells) sequenced per donor (the plotted bars)",
  panel = "Suppl 7B", role = "plotted")
