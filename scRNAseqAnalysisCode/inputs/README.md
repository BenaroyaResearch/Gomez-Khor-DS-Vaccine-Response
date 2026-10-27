# Pipeline inputs

Every file the P635 pseudobulk pipeline reads. Copies, not moves — the originals
are untouched in the Box locations recorded below. The notebooks reference this
directory through `INPUT_DIR <- "inputs"`, so the analysis runs from a checkout
without reaching into Box.

| file | size | MD5 | used by |
|---|---|---|---|
| `P635Data.Rdata` | 6.8 MB | `46f470f92134e9c15c7af456cb9e0532` | step 1 |
| `Mixcr_summary_updated_with_Nucleotidesequences_1stAug2025_PseudogeneFiltered.csv` | 786 KB | `aaa20efd87c01660cbed782d9161a8d4` | step 2 |
| `QC_libids_HCLC_tobeExcluded_Aug2025.csv` | 18 KB | `4da1f11f9d208ee13f7de67836efa428` | step 2 |
| `P635_draw_processed.xlsx` | 9.6 KB | `485e89f895e9c8a0093d2e8e020ea7a0` | step 3 |
| `imgt_vquest/{IgH,IgK,IgL}_{IgG,IgM}.txt` | 428 KB total | (per file, below) | Figure 5D |
| `GSEA_Hallmark_geneset_DysregGenes_IgGIgM.tsv` | 3.9 KB | `38feb70428778d9964205a78c8b3e49c` | Figure 5H |
| `mast_de/DEGs_IGHG_MAST.csv` | 267 KB | `a8aa52fea0af62f211d00c8ac0ec27c3` | Figures 5G, 5H |
| `mast_de/DEGs_IGHM_MAST.csv` | 251 KB | `a69b0be4ff4345d05db558776d5db93a` | Figures 5G, 5H |
| `clonal_expansion_per_donor_prior_analyst.csv` | 1 KB | `7f8c61493c35debec0e0313b87719d58` | Suppl Fig 7D |
| `sterile_transcript_annot_all_donors_Jan2026.csv` | 279 KB | `2abd6aee5c5fccaef8a84b4e21630b17` | Figure 5B |

## Provenance and contents

### `P635Data.Rdata`

Original: `Box/2023 KhorLab DS COVIDvax response/Bioinformatics_analysis_Harshini/data/`

The original `apird` pull for project P635, saved by the prior analyst for
offline work. Contains four objects:

- `counts` — **raw integer gene counts** from `getGeneCounts()`, 864 libraries ×
  58,302 Ensembl gene IDs, stored libraries-as-rows (the notebook transposes it).
  This is the only source of raw counts in the project and the correct input for
  count-based differential expression.
- `anno` — sample annotation per library (donor, study group, timepoint, sex,
  plate, sort).
- `metrics` — sequencing QC metrics per library (total reads, percent aligned,
  median CV of coverage, duplication).
- `bcrs` — **not usable**: only 8 rows. BCR calls come from the MiXCR summary
  below, not from here.

Note: this is plate-based SMART-Seq, so one library equals one sorted cell.

### `Mixcr_summary_updated_with_Nucleotidesequences_1stAug2025_PseudogeneFiltered.csv`

Original: `Box/HChandrashekar@benaroyaresearch.org - Harshini Chandrashekar's Files and Folders/P635_Khor_CovidVacc_DS/data/`

MiXCR BCR annotations, one row per chain: V/J/C gene calls, CDR3 amino acid
(`junction`) and nucleotide (`cdr3_nt_sequence`) sequences, and the full-length
nucleotide sequence. Pseudogene-filtered, as the filename says.

### `QC_libids_HCLC_tobeExcluded_Aug2025.csv`

Original: same folder as the MiXCR summary.

Manual exclusion list of BCR sequences flagged during the prior analyst's
August 2025 QC review, matched on `libid` + `full_nt_sequence`. Only the first
two columns carry data; the rest are empty spreadsheet columns, so the notebook
reads it with `name_repair = "minimal"` and selects those two.

### `P635_draw_processed.xlsx`

Original: same folder as the MiXCR summary.

Draw and processing timestamps per donor, 16 rows. The only source for
`TimeElapsed` (hours from draw to processing) and its binary split, which step 3
shows to be the largest non-depth axis of cell-level variation and which the
step 4 differential expression adjusts for.

### `sterile_transcript_annot_all_donors_Jan2026.csv`

Original: `Box/HChandrashekar@benaroyaresearch.org - Harshini Chandrashekar's
Files and Folders/P635_Khor_CovidVacc_DS/tables/sterile_transcript_annot_all_donors_Jan2026.rds`

Sterile (germline) transcript read counts, 6,408 rows = 801 cells x 8 sterile
loci (`IGHM-S`, `IGHG1-S`..`IGHG4-S`, `IGHA1-S`, `IGHA2-S`, `IGHE-S`), with each
cell's donor, study group and productive isotype. This is the input to **main
text Figure 5B**.

Produced upstream by `P635_Sterile_tx_all_donors_Jan2026.Rmd`, which runs
`sciCSR` over per-donor merged BAM files
(`data/BAM_files_for_sterile_tx_quantification/P635-*.cbub.merged.bam`). Those
BAMs are ~2--6 GB each, ~80 GB for the 16 donors, so they are **not** part of
this repository and that quantification is not rerunnable from a checkout; the
sequence data goes to GEO/SRA. The cached per-cell table above is the
reproducible starting point, exactly as the prior analyst's own plotting code
used it (it re-reads the cached `.rds` rather than recomputing).

Converted here from `.rds` to `.csv` with `na = ""`. The upstream `.tsv` export
of the same object is **not** byte-faithful: 8 rows whose `c_gene` is the
literal string `"NA"` (a paste artefact where one chain was missing) are read
back as missing values. `c_gene_family`, which is what Figure 5B filters on, is
unaffected either way.

### `imgt_vquest/` --- IMGT/V-QUEST mutation statistics

Original: `Box/HChandrashekar@benaroyaresearch.org - Harshini Chandrashekar's
Files and Folders/P635_Khor_CovidVacc_DS/IMGT_VQUEST_outputs/`, the
`8_V-REGION-nt-mutation-statistics.txt` export from each run, renamed here to
`<locus>_<productive isotype>.txt`:

| this repo | V-QUEST run directory |
|---|---|
| `IgH_IgG.txt` | `P635_allsubjects_IGHG` |
| `IgH_IgM.txt` | `P635_allsubjects_IgM` |
| `IgK_IgG.txt` | `P635_allsubjects_IGK_IgG` |
| `IgK_IgM.txt` | `P635_allsubjects_IGK_IgM` |
| `IgL_IgG.txt` | `P635_allsubjects_IGL_IgG` |
| `IgL_IgM.txt` | `P635_allsubjects_IGL_IgM` |

| file | MD5 |
|---|---|
| `imgt_vquest/IgH_IgG.txt` | `1e7e0ee7dbc16769df1230fdf6792cc1` |
| `imgt_vquest/IgH_IgM.txt` | `25cc29cff076227f322075c294ad55a4` |
| `imgt_vquest/IgK_IgG.txt` | `49641d0ab09ae70fcc23f6442aedcddb` |
| `imgt_vquest/IgK_IgM.txt` | `ced2ad8325bebfbb369bbad9d4964d68` |
| `imgt_vquest/IgL_IgG.txt` | `01b8e17b910cbddfe8f9d1d44a2d7a94` |
| `imgt_vquest/IgL_IgM.txt` | `03be2f54a914acdaa44339aa7ba6d097` |

One row per sequence, 1,060 sequences in all. The input to **main text Figure
5D**. Study group and donor are parsed from the FASTA header V-QUEST echoes
back in `Sequence ID` (`..._P635-12_TC`); the panel script cross-checks that
parse against the donor table from step 1 and stops if they disagree.

The V-QUEST runs themselves were submitted to the IMGT web service by the prior
analyst from per-locus FASTA files, so they are not reproducible from a
checkout; these exports are the reproducible starting point.

### `mast_de/DEGs_IGH{G,M}_MAST.csv`

Original: `Box/HChandrashekar@benaroyaresearch.org - Harshini Chandrashekar's
Files and Folders/P635_Khor_CovidVacc_DS/plots_for_figures_Sept2025_timeElapsedCorrected/DEGs_IGH{G,M}_DS_vs_TC_MAST.csv`

Per-cell Seurat `FindMarkers` results, MAST test, DS vs control, run separately
in IgG+ and IgM+ cells with `TimeElapsed_binary` as a latent variable. The
source for **Figures 5G and 5H**.

**Two traps here.**

First, gene symbols sit in an **unnamed first column**, so `read_csv` names it
`...1`. Code that refers to a `gene` column fails against these files as
written.

Second, and more consequential: the same folder also holds
`DEGs_IGH{G,M}_DS_vs_TC_MAST_timebinary_2026.csv`, a later run with a very
similar name. **That is not the run behind Figures 5G and 5H.** Taking the union
of adjusted p < 0.05 across isotypes gives 18 genes from the files above (17 of
which MSigDB mapped, matching the gene list recorded in the Figure 5H export),
against 13 different genes from the `_timebinary_2026` pair. The two runs also
disagree about individual genes: HBB is significant in IgM+ only in the run used
here, which is how the published Figure 5G draws it, but significant in both in
the later run.

Note these are the donor-unaware, per-cell results. The donor-level pseudobulk
analysis of the same comparison is a separate supplemental figure.

### `GSEA_Hallmark_geneset_DysregGenes_IgGIgM.tsv`

Original: same Box folder, `data/`.

MSigDB "Investigate Gene Sets > Compute Overlaps" result for the Hallmark (H)
collection, run on the 17 mappable genes from the list above. The input to
**main text Figure 5H**.

The overlap test was run on the MSigDB **website**, not in R, so it is not
reproducible from a checkout; this export is the reproducible starting point.
Keep the full export rather than the abridged
`GSEA_Hallmark_geneset_DysregGenes_IgGIgM_subsetinfo.txt`: it carries a header
block (collection, gene counts, universe size) and a trailing "Gene/Gene Set
Overlap Matrix" naming every gene MSigDB received. The panel script uses that
matrix to verify, gene by gene, that the enrichment and the DE results still
agree.

### `clonal_expansion_per_donor_prior_analyst.csv`

The 16 per-donor proportions of clonally expanded cells behind **Supplemental
Figure 7D**.

This panel is the prior analyst's analysis, from
`code/ClonalExpansion_analysis_forFigure.Rmd`. A cell is called expanded when it
shares an identical V gene, J gene, C gene and CDR3 with another cell from the
same donor. **It is not reproducible from this repository**: her April 2026 run
recomputed nearest-neighbour distances on a cell set that was never saved, and
the distance tables that were saved are from an earlier (August 2025) run, which
gives p = 0.666 rather than the published 0.829.

**The values here were recovered from the vector coordinates of her published
plot** (`plots_for_figures_FINAL/Clonal_Expanded_proprotions_stripplot_new_threshold_April2026.pdf`),
then verified: they reproduce the published Wilcoxon p = 0.829 and the published
group medians (0.0371 Ctrl, 0.0418 DS) to four decimal places, and the panel
script asserts all three on every run.

Reading numbers back off a figure is a stopgap, not a provenance. **Replace this
file with the prior analyst's own per-donor table when she can supply it**; the
panel script will keep passing its cross-check if the values agree. Donors are
identified only by rank within group, because the plot does not carry donor
labels -- another reason to prefer her table.

## Not included

The prior analyst's QC'd Seurat object
(`plots_after_allQC_26thAug2025P635_seuratobject_after_ALL_QCs_26thAug2025.rds`,
25 MB) is used only by the optional cross-check chunk at the end of step 1,
which verifies that our library QC reproduces her 801-cell set. It is a
comparison artefact rather than a pipeline input, so it is left in Box; the
chunk skips itself when the file is absent. Its `counts` layer is in any case
already library-size normalized and must not be used for differential
expression.
