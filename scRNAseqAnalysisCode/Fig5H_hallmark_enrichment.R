## Main text Figure 5H --- Hallmark gene-set enrichment of the dysregulated genes
##
## "Given the overall concordance, DEGs in either spike-specific IgM+ or IgG+ B
## cells from people with versus without DS were combined and analyzed by
## geneset enrichment (H)."
##
## Input:  inputs/GSEA_Hallmark_geneset_DysregGenes_IgGIgM.tsv
##           MSigDB "Investigate Gene Sets > Compute Overlaps" result, Hallmark
##           (H) collection, run on the combined 17-gene list. The overlap test
##           was run on the MSigDB website, so it is NOT reproducible from a
##           checkout -- this export is the reproducible starting point.
##         inputs/mast_de/DEGs_IGH{G,M}_MAST.csv
##           the per-cell Seurat/MAST results the gene list came from, used here
##           to re-derive that list and check it gene-by-gene against the one
##           recorded inside the MSigDB export.
## Output: pseudobulk_clean_output/figures/Fig5H_hallmark_enrichment.pdf/.png
##         pseudobulk_clean_output/figure_source_data/Fig5H_hallmark_enrichment.csv
##         pseudobulk_clean_output/figure_source_data/Fig5H_input_gene_list.csv
##
## Run:  Rscript figure_panels/Fig5H_hallmark_enrichment.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
suppressPackageStartupMessages(library(stringr))

## ---------------------------------------------------------------------------
## Parameters
## ---------------------------------------------------------------------------

## The gene list submitted to MSigDB: genes differentially expressed in IgG+ or
## IgM+ cells in the per-cell Seurat/MAST analysis (NOT the donor-level
## pseudobulk, which is a separate supplemental figure).
##
## Note which MAST run this is. The same folder also holds
## DEGs_*_DS_vs_TC_MAST_timebinary_2026.csv, a later run; its DEGs are a
## different set and do not reproduce this panel. Figures 5G and 5H both come
## from the DEGs_*_DS_vs_TC_MAST.csv pair at adjusted p < 0.05.
DEG_P_THRESHOLD <- 0.05
DEG_P_COLUMN    <- "p_val_adj"

## Colour ramp for the FDR scale, low to high, as published.
FDR_COLOURS <- c("red", "purple", "blue")

## ---------------------------------------------------------------------------
## The MSigDB overlap export
## ---------------------------------------------------------------------------

tsvPath <- file.path(REPO_ROOT, "inputs", "GSEA_Hallmark_geneset_DysregGenes_IgGIgM.tsv")
raw <- readLines(tsvPath, warn = FALSE)

## The export opens with a metadata block ("Collection(s):", "# genes in
## comparison (n):", ...) before the table header. Find the header rather than
## assuming a fixed number of lines to skip.
headerLine <- grep("^Gene Set Name\t", raw)
if (length(headerLine) != 1) stop("could not find the table header in ", basename(tsvPath))

metaValue <- function(key) {
  hit <- grep(key, raw[seq_len(headerLine - 1)], fixed = TRUE, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(".*\t", "", hit[1]))
}
nGenesSubmitted <- suppressWarnings(as.integer(metaValue("# genes in comparison (n):")))
nUniverse       <- suppressWarnings(as.integer(metaValue("# genes in universe (N):")))
collection      <- metaValue("Collection(s):")
cat(sprintf("MSigDB export: collection %s | %s genes submitted | universe %s\n",
            collection, nGenesSubmitted, format(nUniverse, big.mark = ",")))

## The export continues past the results table with a "Gene/Gene Set Overlap
## Matrix" section, so read only as far as the first blank line after the header.
tableEnd <- headerLine + which(!nzchar(trimws(raw[(headerLine + 1):length(raw)])))[1] - 1
fig5h <- read_tsv(I(raw[headerLine:tableEnd]), show_col_types = FALSE) %>%
  dplyr::rename(geneSet      = `Gene Set Name`,
                nGenesInSet  = `# Genes in Gene Set (K)`,
                description  = Description,
                nOverlap     = `# Genes in Overlap (k)`,
                geneRatio    = `k/K`,
                pValue       = `p-value`,
                fdrQValue    = `FDR q-value`) %>%
  dplyr::mutate(
    ## "HALLMARK_INTERFERON_GAMMA_RESPONSE" -> "Interferon_Gamma_Response"
    label = str_remove(geneSet, "^HALLMARK_") %>%
            str_to_lower() %>%
            str_replace_all("(^|_)([a-z])", function(x) toupper(x)) %>%
            str_replace_all("Il6", "IL6") %>%
            str_replace_all("Jak", "JAK") %>%
            str_replace_all("Stat3", "STAT3") %>%
            str_replace_all("P53", "P53") %>%
            str_replace_all("Tgf", "TGF") %>%
            str_replace_all("Tnfa", "TNFA") %>%
            str_replace_all("Nfkb", "NFKB")) %>%
  dplyr::arrange(fdrQValue) %>%
  ## Lowest FDR at the top of the panel.
  dplyr::mutate(label = factor(label, levels = rev(label)))

print(fig5h %>% dplyr::select(geneSet, nGenesInSet, nOverlap, geneRatio, pValue, fdrQValue) %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())

## ---------------------------------------------------------------------------
## Re-derive the submitted gene list and check it against the export
## ---------------------------------------------------------------------------

readDE <- function(isotype) {
  ## Gene symbols sit in an unnamed first column, so read_csv names it "...1".
  read_csv(file.path(REPO_ROOT, "inputs", "mast_de",
                     sprintf("DEGs_%s_MAST.csv", isotype)),
           show_col_types = FALSE, name_repair = "minimal") %>%
    dplyr::rename(gene = 1) %>%
    dplyr::mutate(isotype = isotype)
}
deG <- readDE("IGHG"); deM <- readDE("IGHM")

sigGenes <- function(d) d$gene[!is.na(d[[DEG_P_COLUMN]]) & d[[DEG_P_COLUMN]] < DEG_P_THRESHOLD]
geneList <- sort(union(sigGenes(deG), sigGenes(deM)))

inputGenes <- tibble::tibble(gene = geneList) %>%
  dplyr::left_join(deG %>% dplyr::select(gene, log2FC_IgG = avg_log2FC,
                                         p_IgG = p_val, padj_IgG = p_val_adj), by = "gene") %>%
  dplyr::left_join(deM %>% dplyr::select(gene, log2FC_IgM = avg_log2FC,
                                         p_IgM = p_val, padj_IgM = p_val_adj), by = "gene") %>%
  dplyr::mutate(
    significantIgG = !is.na(p_IgG) & p_IgG < DEG_P_THRESHOLD,
    significantIgM = !is.na(p_IgM) & p_IgM < DEG_P_THRESHOLD,
    significantIn  = dplyr::case_when(significantIgG & significantIgM ~ "Both",
                                      significantIgG ~ "IgG+", TRUE ~ "IgM+"))

cat(sprintf("re-derived gene list: %d genes at %s < %g in either isotype\n",
            nrow(inputGenes), DEG_P_COLUMN, DEG_P_THRESHOLD))

## Check the list gene by gene, not just by count: the export's overlap matrix
## names every gene MSigDB actually received. MSigDB silently drops symbols it
## cannot map, so the submitted list is allowed to be a superset -- but nothing
## MSigDB used may be missing from it.
matrixHeader <- grep("^Entrez Gene Id\t", raw)
msigdbGenes <- if (length(matrixHeader) == 1) {
  rows <- raw[(matrixHeader + 1):length(raw)]
  rows <- rows[nzchar(trimws(rows))]
  sort(unique(vapply(strsplit(rows, "\t", fixed = TRUE), function(x) x[2], character(1))))
} else character(0)

if (length(msigdbGenes)) {
  missing <- setdiff(msigdbGenes, inputGenes$gene)
  dropped <- setdiff(inputGenes$gene, msigdbGenes)
  if (length(missing))
    stop("MSigDB used genes that the re-derived list does not contain: ",
         paste(missing, collapse = ", "),
         " -- the enrichment and the gene list have drifted apart.")
  cat(sprintf("all %d genes in the MSigDB export are present in the re-derived list\n",
              length(msigdbGenes)))
  if (length(dropped))
    cat(sprintf("  (%s not mapped by MSigDB; its header records n = %s)\n",
                paste(dropped, collapse = ", "), nGenesSubmitted))
  inputGenes <- inputGenes %>% dplyr::mutate(usedByMSigDB = gene %in% msigdbGenes)
} else {
  warning("no gene/gene-set overlap matrix found in the export; gene list not verified")
  inputGenes <- inputGenes %>% dplyr::mutate(usedByMSigDB = NA)
}

## ---------------------------------------------------------------------------
## Panel
## ---------------------------------------------------------------------------

pFig5H <- ggplot(fig5h, aes(x = geneRatio, y = label)) +
  geom_point(aes(size = nOverlap, colour = fdrQValue)) +
  scale_colour_gradientn(
    colours = FDR_COLOURS,
    values  = scales::rescale(c(min(fig5h$fdrQValue), median(fig5h$fdrQValue),
                                max(fig5h$fdrQValue))),
    name = expression(italic(P)[adj])) +
  scale_size_continuous(name = "Gene count") +
  labs(x = "Gene ratio", y = NULL,
       title = "Figure 5H: Hallmark gene set",
       subtitle = sprintf("MSigDB Hallmark overlap | %d genes dysregulated in IgG+ or IgM+ cells",
                          nrow(inputGenes))) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank())

print(pFig5H)
savePanel(pFig5H, "Fig5H_hallmark_enrichment", width = 7.5, height = 4.5)

## ---------------------------------------------------------------------------
## JCI Insight source data
## ---------------------------------------------------------------------------

saveSource(
  fig5h %>% dplyr::transmute(
    geneSet, label = as.character(label), description,
    nGenesInSet, nOverlap, geneRatio, pValue, fdrQValue),
  "Fig5H_hallmark_enrichment",
  "MSigDB Hallmark overlap results for the combined dysregulated-gene list (the plotted points)",
  panel = "Fig 5H", role = "plotted")

saveSource(
  inputGenes %>% dplyr::select(gene, significantIn, usedByMSigDB,
                               log2FC_IgG, p_IgG, padj_IgG,
                               log2FC_IgM, p_IgM, padj_IgM),
  "Fig5H_input_gene_list",
  sprintf("The %d dysregulated genes behind the enrichment (MAST %s < %g in IgG+ or IgM+ cells)",
          nrow(inputGenes), DEG_P_COLUMN, DEG_P_THRESHOLD),
  panel = "Fig 5H", role = "supporting")
