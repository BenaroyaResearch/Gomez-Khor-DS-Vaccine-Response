# IMPACD, Level 1, and modifier analysis

Source code for IMPACD gating, Level 1 statistics, and modifier analysis of CyTOF, flow T cell, and flow B cell data.

The original code files were copied without changes. This collection preserves the research notebooks and shared functions, including historical and exploratory sections within those files. Data, results, plots, and caches are excluded from this directory.

## Source

Original repository: [BenaroyaResearch/Khor_covidvax_response_IMPACD](https://github.com/BenaroyaResearch/Khor_covidvax_response_IMPACD)

Source checkout HEAD recorded during transfer: `0a35cd8006072e3af8a48971408da3d1d6c115ed`

The original project documentation is preserved in [docs/SOURCE_README.md](docs/SOURCE_README.md). It contains historical progress notes that may describe earlier implementations.

## Analysis file map

| Data | Analysis notebook | Shared analysis functions |
| --- | --- | --- |
| CyTOF | `analysis_CyTOF_RG_2025_05.Rmd` | `R/cytof_analysis_RG_2025_05.R` |
| Flow T cells | `analysis_FLOW_3_RG_2025_05.Rmd` | `R/flow_4_analysis_RG_2025_09.R` |
| Flow B cells | `analysis_FLOW_Bcell_RG_2025_05.Rmd` | `R/flow_Bcell_analysis_RG_2025_11.R` |

`R/common.R` contains the shared Level 1, multiple testing, and spreadsheet export functions.

Antibody titre analysis is preserved in `analysis_vax_response_RG_2025_06.Rmd` and `R/abTiter_analysis_RG_2025_06.R`. Its modifier analysis compares model slopes and differs from the DS vs TC distribution comparisons described below.

## Gating and object preparation

1. CyTOF gating is in `gating_CyTOF_RG_2025_10.Rmd` and `gating_CyTOF_RG_2025_10_TEMRA.Rmd`, supported by `R/cytof_gating_RG_2025_10.R`. The two notebooks retain different selections of roots for rerunning gating.
2. Flow T cell gating is in `gating_flow_T_cells_RG_2025_09.Rmd`. Additional naive CD4 AIM++ gating is in `gating_flow_T_cells_Tn_AIMpp_RG_2026_01.Rmd`.
3. Flow B cell gating is in `gating_flow_B_cells_RG_2025_08.Rmd` and the accompanying PBplus, SM, NSM, Tn, and DN gating notebooks.
4. `Build_full_Flow_Bcell_objects_2025_09.Rmd` assembles the B cell analysis objects from the gated counts.

The analysis inputs contain `pctParent`, `popCounts`, `parentCounts`, `nodes`, and `design`. These provide population frequencies, event counts, parent counts, gating relationships, and sample annotations.

## Level 1 analysis

The DS vs TC analysis functions fit models to eligible populations and calculate contrasts. Their default Gaussian model uses `pctParent ~ studyGroup`, weighted by parent counts.

`getLevel1()` in `R/common.R` selects `depth == 1` and recalculates BH adjusted p values within each comparison using that Level 1 subset.

The Level 1 table and the full statistics table therefore use different multiple testing sets. The DS vs TC modifier functions receive the full statistics table, rather than the separately adjusted Level 1 table.

## DS vs TC modifier analysis

The modifier functions first select populations using the primary statistics and count filters.

For each selected population, DS frequencies are centered by subtracting that population's TC median frequency. Parent and child distributions are then compared using `wilcox.test(..., exact = FALSE)`.

Round 1 selects children with the same terminal node at the next depth. For parents with modifiers, the existing modifier sequence must also match.

Round 2 tests surviving depth 3 populations against additional depth 2 candidates matched using the child's `modifier_2` and the candidate parent's `modifier_1`.

The code retains modifier hits and comparison records, including Round 2 kept and kicked results.

The Wilcoxon comparisons use raw p values against the threshold named `fdrCutoff`. These functions do not apply an additional BH correction to those Wilcoxon p values. This filtering step should be distinguished from the BH correction of the primary statistics and Level 1 tables.

## Before running

1. Configure input and output paths. Several files retain absolute paths or `setwd()` calls pointing to the original repository. `R/project_settings.R` also retains the original base directory.
2. Set the working directory to this `IMPACD` directory after updating those paths.
3. Supply the required FCS files, FlowJo workspaces, annotations, and processed count objects for the selected workflow.
4. Install the packages loaded by that workflow. Gating uses IMPACD, flowWorkspace, and CytoML; analysis also uses the statistical, plotting, and export packages named in the scripts.
5. Load one assay's analysis helper at a time. The helpers define overlapping function names.
6. Review notebook chunks and active root selections before running. Some selections represent targeted reruns rather than a complete analysis.
7. Legacy aging chunks in the three main analysis notebooks reference `cytof_functions.R`, which is not included in this directory. Those chunks require correction or exclusion before rendering the whole notebook.

File copies were checked against their sources using `cmp`. The scientific analyses were not rerun during this transfer.
