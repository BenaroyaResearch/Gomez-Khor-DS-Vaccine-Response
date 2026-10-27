# Figure panels

One script per manuscript figure panel. Each is a thin layer over the
`P635_pseudobulk_step1..4` pipeline: the pipeline does the analysis and saves
intermediates under `pseudobulk_clean_output/`, and each script here reads one
intermediate and writes

- the panel, to `pseudobulk_clean_output/figures/`, and
- its source-data table(s), to `pseudobulk_clean_output/figure_source_data/`
  (JCI Insight requires per-panel source data).

Exactly one script owns each panel's output files.

## Running

Knit the pipeline steps a panel depends on, then run the panel script from the
repo root:

```
Rscript figure_panels/Fig5A_isotype_frequencies.R
```

| Panel | Script | Depends on |
|---|---|---|
| Fig 5A | `Fig5A_isotype_frequencies.R` | step 1, step 2 (`step2_cell_metadata.csv`) |
| Fig 5B | `Fig5B_sterile_transcripts.R` | `inputs/sterile_transcript_annot_all_donors_Jan2026.csv` (no pipeline step) |
| Fig 5D | `Fig5D_somatic_hypermutation.R` | `inputs/imgt_vquest/` + step 1 (donor table, for a group-label cross-check) |
| Fig 5G | `Fig5G_isotype_concordance.R` | `inputs/mast_de/` |
| Fig 5H | `Fig5H_hallmark_enrichment.R` | `inputs/GSEA_Hallmark_geneset_DysregGenes_IgGIgM.tsv` + `inputs/mast_de/` |
| Suppl 7B | `SupplFig7B_cells_per_donor.R` | step 2 (`step2_cells_chains.rds` `$libAnno`) |
| Suppl 7C | `SupplFig7C_chains_per_cell.R` | step 2 (`$libAnno`) |
| Suppl 7D | `SupplFig7D_clonal_expansion.R` | `inputs/clonal_expansion_per_donor_prior_analyst.csv` |
| Suppl 7D | `SupplFig7D_repertoire_metrics.R` | step 2 (`$repertoireMetrics`) |
| Suppl 7E | `SupplFig7E_shm_subtypes.R` | `inputs/imgt_vquest/IgH_IgG.txt` + step 1 |
| Suppl 7H | `SupplFig7H_mast_volcano.R` | `P635_step4_DE_mast_combined.Rmd` |
| Suppl 7I | `SupplFig7I_cameraPR_genesets.R` | `P635_step4_DE_mast_combined.Rmd` |

## Source data for submission

JCI Insight requires the values behind every graph and every reported mean, with
a **separate tab per figure panel**, in one Excel file kept separate from other
supplemental material (<https://www.jci.org/kiosks/publish/supporting-data>).

On this manuscript the corresponding author maintains that master workbook, so
this repo emits **one file per panel** for copying into it. Each file holds a
single worksheet already named as the tab should be named in the master
("Fig 5A"), so copying it across — right-click the sheet tab > *Move or Copy* >
choose the master > tick *Create a copy* — carries the tab name with it.

```
Rscript figure_panels/Fig5A_isotype_frequencies.R    # writes CSVs + manifest rows
Rscript figure_panels/Fig5B_sterile_transcripts.R
Rscript figure_panels/build_supporting_data_values.R # -> supporting_data_values/Fig5A.xlsx, ...
```

Set `OUTPUT_MODE` in the builder to `"combined"` or `"both"` to also emit a
single all-panel workbook, for a journal or coauthor that wants one file.

`saveSource()` registers every table it writes in
`figure_source_data/_manifest.csv` (panel, role, description, dimensions), and
the builder reads that manifest. **The manifest is generated, never edited by
hand** — add a panel script and its tabs appear in the workbook automatically.

Each table declares a `role`, which orders the blocks within a panel's tab:

| role | meaning |
|---|---|
| `plotted` | the values behind the graph itself — what JCI actually requires |
| `supporting` | intermediate values a reported test is computed from |
| `statistics` | the test results annotated on the panel |

One tab per panel, with each block labelled by its description. A statistic that
does not exist is written as an empty cell with a `testable = FALSE` flag beside
it, never as `NaN` and never as a fabricated value.

The Methods Data Availability section must carry the callout JCI requires:
*"Values for all data points in graphs are reported in the Supporting Data
Values file."*

## Conventions

- `_shared.R` holds the paths, the published group colours, and the
  `savePanel()` / `saveSource()` writers. Source it first.
- Study-group codes are `TC` and `DS`; the manuscript calls `TC` "Ctrl", and the
  source tables carry both (`studyGroup`, `groupLabel`).
- Source tables are tidy and carry the raw counts and denominators, not just the
  derived percentage, so every plotted value can be recomputed from the table.
- Where a panel's displayed statistic is still under discussion, the choice is a
  named parameter at the top of the script, and *all* candidate statistics are
  written to the statistics table regardless of which one is drawn.

## Known discrepancies with the submitted figures

### Figure 5A

`Fig5A_isotype_frequencies.R` reproduces the submitted panel's data exactly
(all 64 plotted points, verified against the figure's vector coordinates), with
two deliberate differences:

1. The submitted panel divided donor **P635-1**'s isotype counts by 51 rather
   than that donor's own denominator of 53, so its four percentages summed to
   103.9% instead of 100%. The script asserts the per-donor sum, so this cannot
   recur silently.
2. The submitted panel's printed **IgD (p=0.84)** and **IgA (p=0.593)** values do
   not correspond to the plotted data under any test. Its printed IgM and IgG
   values match a Welch t-test exactly.

**Decided: the panel now reports a Wilcoxon rank-sum test** (`TEST <- "wilcoxon"`),
which is what the Figure 5 legend already states and what the paper's other
panels use. The reported values become:

| isotype | published | Wilcoxon (now reported) | Welch t |
|---|---|---|---|
| IgD | 0.84 | 0.344 | 0.586 |
| IgM | 0.0001 | 0.000155 | 0.000142 |
| IgG | 0.001 | 0.00295 | 0.000882 |
| IgA | 0.593 | 0.411 | 0.714 |

IgM and IgG remain strongly significant and IgD and IgA remain non-significant,
so no conclusion changes. The statistics table carries both tests and the
BH-adjusted value for each, whichever is drawn.

`CENTER_LINE` is still `"mean"`, matching the submitted panel. A median line
would pair more naturally with a rank-sum test and with the paper's other
legends ("lines denote median values"); it is a one-word change, and it moves
the IgG Ctrl line from 62.0 to 65.6 and the IgA DS line from 3.4 to 0.

### Figure 5B

`Fig5B_sterile_transcripts.R` reproduces the submitted panel exactly: all five
donor-level Wilcoxon p-values match (IgM 0.875, IgG1 0.68, IgG2 0.945, IgG3 NaN,
IgG4 1), on the same 219 IgM+ cells. Two things to be aware of:

1. The submitted panel prints **"p=NaN"** for IgG3. That is a real test failure,
   not a formatting slip: every donor's median is 0 in both groups, so the
   rank-sum statistic has zero variance. `UNTESTABLE_LABEL <- "explicit"` prints
   "not testable" instead.
2. More importantly, **the published summary statistic has almost no power at the
   IgG loci.** Only 16--48% of cells carry any sterile read there, so the median
   cell is 0 and nearly every donor's median is 0 — the test compares 0 against 0.
   That is what produces p=1 and p=NaN, and those values are not evidence of
   absence. Summarising each donor by the **mean** instead of the median, on the
   identical data, gives IgG3 p=0.008 and IgG4 p=0.024 (nominal, unadjusted;
   BH across the five loci gives q≈0.04 and q≈0.06), both in the direction of
   *more* sterile transcription in DS. The statistics table carries `donorMean`
   alongside `donorMedian` so this is checkable; the per-donor non-zero
   percentages are in the statistics table for the same reason.

This is a change in statistical approach, so it is flagged rather than applied.

### Figure 5D

`Fig5D_somatic_hypermutation.R` reproduces all six submitted p-values exactly
(IgG+: IgH 0.381, IgK 0.876, IgL 0.862; IgM+: IgH NaN, IgK 0.351, IgL 0.861).

The test is a **Welch t-test on per-donor medians**, not the Wilcoxon that the
same upstream notebook also computed and saved alongside it
(`mutation_stats_Wilcox_*.csv` vs `mutation_stats_*.csv`). Picking the wrong one
of that pair gives 0.791 where the figure shows 0.381, so the distinction
matters; the statistics table records which test was reported.

Two things to be aware of:

1. **IgH / IgM+ is again untestable** — every donor's median is 0 in both groups,
   because IgM+ cells are largely unmutated. The submitted panel prints
   "p = NaN"; `UNTESTABLE_LABEL` offers `"na"` or `"explicit"` instead. This is
   the same pathology as Figure 5B's IgG3, and the same caveat applies: it is not
   evidence of absence.
2. **IgL / IgM+ has 7 control donors, not 8** — one control contributed no IgL
   sequences in IgM+ cells. The per-panel donor counts are in the statistics
   table rather than being assumed to be 8 and 8.

### Figure 5H

`Fig5H_hallmark_enrichment.R` reproduces the submitted panel exactly: the same
nine Hallmark sets in the same order, with the same gene ratios, gene counts and
FDRs.

Figures 5G and 5H do **not** come from the step1--step4 pseudobulk pipeline.
They come from the prior analyst's per-cell Seurat/MAST results
(`inputs/mast_de/`). `pseudobulk_clean_output/` happens to contain
`Fig5G_isotype_concordance.*` and `Step4_GSEA_hallmark.*` from our own parallel
pseudobulk analysis; those are a different analysis and are not these panels.

The enrichment itself was run on the MSigDB website, so the panel script cannot
recompute it. What it does instead is re-derive the gene list from the MAST
results and check it **gene by gene** against the overlap matrix inside the
export, failing if MSigDB used a gene the DE results no longer support. That
catches drift between the two without pretending the web step is reproducible.

An earlier attempt matched MSigDB's gene *count* (17) from the wrong MAST run at
a raw-p threshold while sharing only 10 of the 17 genes, which is why the check
compares identities rather than counts.

### Figure 5G

`Fig5G_isotype_concordance.R` reproduces the submitted panel point for point.
Every marker's position was read out of the figure's vector content stream and
matched to a gene: 11 IgG+ (black down-triangles), 3 IgM+ (red squares) and 1
Both (green up-triangle, SON), all 15 agreeing with the recomputed fold changes
to within the calibration error.

Same source as Figure 5H: `inputs/mast_de/DEGs_IGH{G,M}_MAST.csv`, adjusted
p < 0.05. Both panels therefore move together if that threshold or that DE run
is ever changed.

Three of the 18 differentially expressed genes (HLA-A, HLA-DQB2, YBX3) were
tested in one isotype only, so they have no second fold change and cannot be
placed on the plane. They are kept in the source table with `plotted = FALSE`
rather than silently dropped.

**Name collision, now resolved.** `P635_pseudobulk_step4_DE.Rmd` used to write a
panel and tables called `Fig5G_isotype_concordance` from our *pseudobulk*
concordance analysis — a different analysis that is not this figure. Those
outputs are now named `Step4_pseudobulk_isotype_concordance` so the two cannot
overwrite each other. The pseudobulk version belongs with the supplemental
pseudobulk figure.

## Supplemental Figure 7

Panels B, C, D, E, H and I. (A is serology; F and G are FLAIRR-seq, handled
elsewhere.)

**The figure was relettered on 2026-10-21** and now agrees with the legend and
the manuscript text. Panel D covers the whole second row -- the clonal-expansion
plot *and* the six repertoire metrics -- so two scripts write into it. The
mutation-subtype panel is E, FLAIRR-seq is F-G, and the MAST volcano and
cameraPR sets are H and I. An earlier draft lettered these F, G-H, I-J with no
E at all; any script or note still using that scheme is stale.

Verified against the relettered figure:

- **7D (clonal expansion)** reproduces the published Wilcoxon p = 0.829 and both
  group medians. See `inputs/README.md` -- the values are recovered from the
  prior analyst's plot and should be replaced by her own table.
- **7D (repertoire metrics)** reproduces exactly: 0.325/0.325/0.269/0.325/0.959/0.878
  against the published 0.33/0.33/0.27/0.33/0.96/0.88.
- **7E** reproduces exactly: 0.082/0.684/0.855/0.535 against the published
  0.08/0.68/0.86/0.53. The legend specifies a Welch t test for this panel, and
  the data are the IgH sequences from IgG+ cells -- the same sequences as the
  top-left panel of Figure 5D, broken down by mutation type.
- **7H** reproduces exactly: the same eight labelled genes and four-key legend.
- **7I** reproduces exactly: the same ten sets in the same order, with
  KRAS_SIGNALING_DN the only one past the FDR threshold.

The vector-reading trick used on Figure 5 does not work on this figure as a
whole: it carries its markers in ~15,700 form XObjects whose streams the PDF
does not expose. Gridlines, centre bars and text are still readable. Individual
upstream plots saved as their own PDFs (such as the clonal-expansion stripplot)
are plain vector and can be read point by point.
