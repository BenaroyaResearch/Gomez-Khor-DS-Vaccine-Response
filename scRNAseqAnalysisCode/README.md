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


