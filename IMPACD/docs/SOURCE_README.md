# Table of Contents
- [Table of Contents](#table-of-contents)
- [Covid vax response in DS](#covid-vax-response-in-ds)
  - [Organization](#organization)
    - [Folders that are not synced to GitHub](#folders-that-are-not-synced-to-github)
- [Workflow](#workflow)
  - [Gating](#gating)
  - [Analysis](#analysis)
    - [Data sets](#data-sets)
      - [Flow data](#flow-data)
      - [CyTOF data](#cytof-data)
    - [Modifier analysis](#modifier-analysis)
    - [DS vs TC](#ds-vs-tc)
      - [Analysis steps](#analysis-steps)
      - [Tables that are generated](#tables-that-are-generated)
      - [Plots that are generated](#plots-that-are-generated)
      - [Model](#model)
      - [Modifier analysis](#modifier-analysis-1)
    - [Aging and Immune age](#aging-and-immune-age)
      - [Analysis steps](#analysis-steps-1)
      - [Tables that are generated](#tables-that-are-generated-1)
      - [Plots that are generated](#plots-that-are-generated-1)
      - [Model](#model-1)
      - [Modifier analysis](#modifier-analysis-2)
    - [Ab Titer response](#ab-titer-response)
      - [Analysis steps](#analysis-steps-2)
      - [Tables that are generated](#tables-that-are-generated-2)
      - [Plots that are generated](#plots-that-are-generated-2)
      - [Model](#model-2)
      - [Modifier analysis](#modifier-analysis-3)
  - [What has been done and next steps](#what-has-been-done-and-next-steps)
    - [Angle of attack](#angle-of-attack)
    - [New work](#new-work)
    - [Other ideas](#other-ideas)



# Covid vax response in DS
This is the code repository for the Khor groups "Covid vax response" project.

The project supports `renv` for package management. To install the packages used in the project, run the following command in the R console:

```R
renv::restore()
```

If this project gets downloaded from Github, it might be necessary to restore some of the folders that are too large to be synced to github. These folders are mentioned below in the section "Folders that are not synced to GitHub". They can be found at `nfs/spribitzer/projects/Khor_covid_vax_response/` on the pure storage.

## Organization
The repository is roughly structured as follows:

- `data/`: Contains `Rdata` objects that are used in the analysis.
- `meetings and slides/`: Contains the slides and notes from our meetings.
- `R/`: Contains the shared R functions that are used in the analysis.
- `misc code/`: Contains miscellaneous, outdated code that was used at some point in the project. A good place to look up various flowworkspace functions that could come in handy.
- `results/`: Contains the results of the analysis such as plots and tables, organized by the type of analysis (e.g. Abtiter response, aging effects)
- `plots/` and `results_old/`: Was the old folder structure where we stored plots and results (tables) in separate folders. This turned out to be not optimal for sharing, so we are in the process of moving everything to the new structure `results/`. Until all the analysis have been rerun, it is probably wise to keep these folders around so that we can recover plots and tables if needed.
- `renv/`: Contains the `renv.lock` file that specifies the packages used in the project.
- `tables/`: This folder contains various tables and spreadsheets that were shared by collaborators, such as meta data, counts from FlowJo counting, etc.

The main `RMD` analysis files are kept in the root directory of the repository. The are somewhat split by what they do:
- Files starting with `gating` were used to IMPACD gate the data. The gated data was then saved as `Rdata` objects in the `data/counts` folder. These objects should contain all the data (meta data, gating path, subsets, counts,...) necessary to run down stream analysis.
- Files starting with `'analysis` are used to run the analysis on the gated data. They load the `Rdata` objects from the `data/counts` folder and (batch) run the analysis. The results are saved in the `results/` folder.

### Folders that are not synced to GitHub
There are a few folders that are not synced to the github repository due to them being too large. These are:
- `analysis_results/`: This folder should contain `Rdata` files that contain the worksspace of the "analysis" of each root path. These are being created as part of the functions that get run in the analysis files. The idea behind this is that it will be easier to regenerate specifc graphs or check for specific values if we have the workspace saved and don't need to rerun the whole analysis.
- `data/raw_counts/`: The rawcounts as they are returned after gating. These get processed and turned into the `Rdata` objects in the `data/counts` folder. The raw counts are not synced to github as they are too large.
- `fcs_files/`: The fcs files that were used as the basis for gating. They originally came from the box folder that was shared by Bernard/Amber. They can be accessed there or on the pure storage folder under: `nfs/spribitzer/projects/Khor_covid_vax_response/fcs_files`

# Workflow
The general workflow for the project is as follows:
1. Gating - generate the events data for the analysis and save it. This should all be done at this point
2. Analysis - run the different types of analysis on the gated data. This is currently still in process. More details on how to run the analysis can be found below in [Analysis](#Analysis)
   - **DS vs TC**: Looks for differences between the DS and TC groups
   - **Aging and Immune age**: Explore populations that change with age. The ones from TC are used to build an immune age model, and that is then evaluated against the DS data.
   - **Abtiter response**: Looks at the response to the vaccine in terms of the antibody titers - this is our vaccine response.
   - **Supervised signatures**: Bernard provided a set of supervised signatures. We compared these to our data to see if they also show effects on the immune response and aging in our data.
## Gating
At this point (2025-03-06) the gating steps should be more or less done. The files starting with `gating_` provide a good starting point for future work. 
> Note: The Tcell and Bcell flow data was a very early version of the gating and might not be "the most efficient way" to gate the data anymore. The most mature gating script is probably the `gating_CyTOF.Rmd` file.

In general there are two major steps to the gating work;
1. Validate the flowJo worspace data and verify that we can load the gates and batches correctly in R.
   - The collaborators will provide a flowJo workspace and a spreadsheet with event counts for a selected set of populations and samples (at least one per batch)
   - Extract gates and transformations from the flowJo workspace
   - Find the fcs files and apply the transformationts to them
   - Then gate the populations that are in the spreadsheet and compare the counts to the spreadsheet
   - If the counts match (perfect diagonal), we can be confident that we are gating correctly
   - Possible issues:
     - Logical gates: sometimes we need to create additional logical gates in R that are not present in the flowJo workspace, e.g. the CD45RO- gate needs to be created from the CD45RO+ gate
     - Transformation: sometimes the transformation is not correctly applied in R. This can be due to the transformation not being correctly extracted from the flowJo workspace or the transformation not being correctly applied in R. Also, batch specific transformations are a thing!
2. Gating the data
   - Once we are confident that we can gate the data correctly, we can gate all the data using IMPACD
   - Identify the rootpaths that we are looking for, gate and save the files.

## Analysis
### Data sets
We have access to multiple data sets that we can use for the analysis, the biggest difference between flow and CyTOF data are the available time points and that flow has an additional antigen specifc compartment (which can somewhat be treated like its own "rootpath").

Here is an overview of the time points we have
| Time Point | Description                           |
|------------|---------------------------------------|
| T1         | Baseline                              |
| T2         | 2 weeks after vaccination             |
| T3         | Approximately 2 months after vaccination |


#### Flow data
The flow data was generated for B cells and T cells, and is stored in the `data/counts/flow` folder. For the B cell data we only have one subset, while we have several for the T cells (TFH, non TFH, Th1....).

The biggest differences between the flow and the CyTOF data are 
- Flow data has an antigen specific and global compartment. This means that we need to address this in the analysis functions, either by writing separate functions for the flow and CyTOF data or else by adding an argument to the functions that specifies the compartment.
- The flow data has samples for a lot more time points (T1, T2, T3)

#### CyTOF data
The CyTOF data was generated almost exclusitely only for T1 samples. However, there are some additional samples that were generated for the T2 and T3 time points. Make sure to subset appropriately when running the analysis. The files are stored in the `data/counts` folder. There are a lot more markers in the CyTOF data, so the data is a lot more complex than the flow data.

### Modifier analysis
The modifier analysis is a necessary step to reduce our "false discovery" rate. The idea is to identify the main sources of variation in the data and then adjust for them in the analysis. For example, if there is a subset `A` that shows a strong difference between DS and TC, we would also expect to see that or a similar difference in any of the subsets that have modifiers added to `A`, like `B/A`, `C/A`, `B/D/A`. The goal of the modifier analysis is to identify where the drivers for differences are actually originating from and to remove any "false" differences that are driven by the modifiers from the results. For example in the above case, if the effects between DS and TC are the same in all `A`, `B/A`, `C/A`, `B/D/A`, we would only report the difference in `A` and not in the other subsets.

However, there a few things to keep in mind when running the modifier analysis:
- If the difference in `A` and `B/A` is the same, but the difference in `C/A` is different, we would still report the difference in `A` and `C/A`. This is because the difference in `C/A` is not driven by the same effect as in `A` and `B/A`.
- The modifier analysis is somewhat finicky and still not "fully fledged out"
- Depending on what is being modeled (DS vs TC, aging), the way the modifier analysis works is completely different. More details can be found in the discussion of the individual analysis types below.

### DS vs TC
The DS vs TC analysis is the most straight forward analysis. We are looking for differences between the DS and TC groups. The analysis is run on the gated data and the results should be stored in the `results/DS_vs_TC` folder. This analysis has been run for the flow data and CyTOF data. For the CyTOF data look at the file `analyis_CyTOF.Rmd`. This file contains the DS vs TC comparisons as well as the aging and immune age analysis. For the flow data, this was run a long time ago and should probably be rerun using shared functionality (models, etc.) like the CyTOF data.

#### Analysis steps
1. Load the data
2. Build the model
3. Test for differences
4. Run the modifier analysis
5. Generate plots and tables

#### Tables that are generated
- All stats, p-values, etc. for the DS vs TC comparison before the modifier analysis
- Level 1 stats, basically the data from above subset to level 1 depth and fdr corrected for this sample size
- Modifier results: The modifier analysis should generate several tables, one that shows what survived after correction, one that shows all the comparisons that have been made and what the outcomes were.

#### Plots that are generated
- Boxplot of all subsets at level 1 depth, indicate the significant ones with stars
- Heatmap of the significant subsets at level 1 depth
- Boxplots for the subsets that survived the modifier analysis (not yet implemented)
- Heatmap for the significant subsets that survived the modifier analysis
- Volcano plot of the modifier analysis results (not yet implemented)
- Radial plots that show what terminal nodes are the most important for the differences between DS and TC (not yet implemented)
  
#### Model
For now we are using a linear response model for the DS vs TC comparison. This is a simple model that is easy to interpret and can be used to generate the p-values for the differences between DS and TC. The model is defined as follows:
```R
fit = glm(pctParent ~ studyGroup*ageAtDraw, data = data, weights = data$weights)
```
For the reimplemenation of the flow data, we should probably use this same model, but loop over compartment and time points. Other models to consider would be a beta regression model or other ways to restrain the model to the 0-1 range.

#### Modifier analysis
For the comparison between DS and TC, the modifier analysis is relatively simple. We do the following:
- start at the significant hits from level 1
- pick a terminal node, e.g. `A`
- find all the level 2 subsets that contain `A` as their terminal node, e.g. `A`, `B/A`, `C/A`
- center the DS data around the mean of the TC data for each terminal node, this step is necessary so that we can compare "distributions" and not "absolute values". absolute values of pctParent are not comparable between different terminal nodes.
- test if the recentered data of `A` is different from the recentered data of `B/A`, `C/A`, right now we are doing a wilcox test `wilcox.test(background_corrected_x, background_corrected_y, exact = FALSE)$p.value`, but it might be worth looking into other tests or doing a model and then comparing the coefficients.
- report the level 2 results only if the test above is showing a significant difference
- repeat for level 2
  
**Problems:** I have noticed that the current implementation is not careful enough at handling situations where `A`, `B/A`,`C/B/A` are all significant and show the same effect. In this case we would only want to report `A` and not the other two, but instead of what we are currently getting is `A` and `C/B/A`. The problem is that we are losing track of of the level 2 population.

One approach to this could be to compare level 1 to level 2 AND level 3. Then go to level 2 and compare against level 3

#### CyTOF Modifier Analysis Updates (May 2025):

A. Structural and workflow changes:
Directory structure updated with suffix _RG_2025_05 to differentiate current (May 2025) and ongoing analyses from previous implementations, making subset-specific results easier to locate and avoiding confusion with past results.

Error-handling (tryCatch) added to the main analysis loop in analysis_CyTOF_RG_2025_05.Rmd, ensuring that errors with individual subsets don't halt execution.

B. Modifier analysis improvements (Round 1):
Original problem: The original code for selecting child subsets during modifier analysis matched subsets based solely on two criteria:

```R
children <- dplyr::filter(
  stats_,
  terminal_node == population_$terminal_node,
  depth == population_$depth + 1
)
```
This approach identified child subsets that shared the same terminal node and were exactly one gate deeper, but did not ensure that modifiers matched sequentially. As a result, parent subsets such as CD4Teff/CXCR5+/CCR4+ could incorrectly be compared against child subsets like CD4Teff/CCR6+/TIGIT+/CCR4+, even though modifier positions didn't align exactly.

This mismatch is biologically irrelevant because the sequential order of modifiers directly reflects the gating strategy used in CyTOF analyses, where each subsequent gate (modifier_1, modifier_2, etc.) indicates a biologically distinct subpopulation. Comparing populations without strictly enforcing this modifier sequence violates the hierarchical logic of gating and produces misleadingly significant results, diluting the true biological meaning of the identified differences.

Implemented solution: Re‑created the original selection: candidate children are one gate deeper and share the same terminal node. Only add extra tests when the parent already has at least one modifier (depth ≥ 2). Level‑1 parents behave exactly as before. Narrows the child set to those whose modifier_1 … modifier₍d–1₎ perfectly match the parent’s. This enforces sequential‑modifier logic: a level‑3 child must share its parent’s full modifier history.
Code lines (cytof_analysis_RG_2025_05.R):

```R
children <- stats_ %>%
  dplyr::filter(
    terminal_node == population_$terminal_node,
    depth == population_$depth + 1
  )
if (population_$depth >= 2) {
  for (m in seq_len(population_$depth - 1)) {
    mod_col <- paste0("modifier_", m)
    children <- children %>%
      dplyr::filter(.data[[mod_col]] == population_[[mod_col]])
  }
}
```
C. Modifier analysis second-round filtering:

Problem being solved: the 1st round of modifier analysis kicked all level 2 subsets that were not not significantly different from level 1, but non-biologically relevant level 3 subsets could still be called as significant hits based upon the idea that Root/M1/M2/T = Root/M2/M1/T and the 1st round code fix of modifier analysis, as implemented, did not analyze a Root/M1/M2/T child against both Root/M1/T and Root/M2/T parents.

Implemented solution:

Round-2 re-tests depth-3 subsets from Round 1, matching on child's modifier_2 against parent's modifier_1, irrespective of order.

Distinguishes subsets as "kicked" or "kept."

Code snippet (cytof_analysis_RG_2025_05.R):
```R
if (max_depth >= 3) {
  depth2_parents <- dplyr::filter(stats_, depth == 2)
  depth3_round1  <- dplyr::filter(report, depth == 3)

  kicked_round2 <- data.frame()
  kept_round2   <- data.frame()

  for (iChild in seq_len(nrow(depth3_round1))) {
    child <- depth3_round1[iChild, ]

    candidate_parents <- depth2_parents %>%
      dplyr::filter(
        terminal_node == child$terminal_node,
        modifier_1 == child$modifier_2
      )

    if (nrow(candidate_parents) == 0) {
      kept_round2 <- rbind(kept_round2, child)
      next
    }

    best_p <- 1
    best_parent <- NA_character_

    for (ip in seq_len(nrow(candidate_parents))) {
      parent <- candidate_parents[ip, ]
      x <- pctParent_study_corrected[parent$subset , weights_study[parent$subset , ] > 0]
      y <- pctParent_study_corrected[child$subset , weights_study[child$subset , ] > 0]

      pval <- wilcox.test(x, y, exact = FALSE)$p.value
      if (pval < best_p) {
        best_p <- pval
        best_parent <- parent$subset
      }
    }

    child$modifierAnalysisPvalue_round2 <- best_p
    child$modifierAnalysisAgainst_round2 <- best_parent

    if (best_p <= fdrCutoff) {
      kept_round2 <- rbind(kept_round2, child)
    } else {
      kicked_round2 <- rbind(kicked_round2, child)
    }
  }

  # Export kicked and kept results immediately
  library(openxlsx)
  exportDir <- file.path("results_RG_2025_05", "CyTOF", analysis, "kicked_kept")
  dir.create(exportDir, recursive = TRUE, showWarnings = FALSE)

  write.xlsx(kept_round2, file = file.path(exportDir, paste0(cellSubset, "_kept_round2.xlsx")), overwrite = TRUE)
  write.xlsx(kicked_round2, file = file.path(exportDir, paste0(cellSubset, "_kicked_round2.xlsx")), overwrite = TRUE)

  report <- dplyr::bind_rows(
    dplyr::filter(report, depth < 3),
    kept_round2
  )

  removed_from_report <- rbind(removed_from_report, kicked_round2)
}
```

D. Wilcoxon test error-handling:
Added error catching around Wilcoxon tests to prevent crashes due to unexpected or invalid subsets (e.g., duplicated HLA-DR+/HLA-DR+ subsets). Safely assigns NA values, ensuring analysis completion.
```R
pval <- tryCatch({
  wilcox.test(x, y, exact = FALSE)$p.value
}, warning = function(w) {
  message("Wilcoxon warning: ", w$message)
  NA_real_
}, error = function(e) {
  message("Wilcoxon error: ", e$message)
  NA_real_
})
```
E. Statistical model simplification (interaction removal):
Original problem: Prior statistical model unnecessarily included interaction terms (sex, ageAtDraw), despite matched cohorts.

Implemented solution:
Simplified model removes interactions, now testing differences based solely on studyGroup.

Original model:
```R
fit <- glm(pctParent ~ studyGroup * sex * ageAtDraw, data = data, weights = data$weights)
```
Updated simplified model:
```R
fit <- glm(pctParent ~ studyGroup, data = data, weights = data$weights)
```
Aligns analysis with cohort matching design, producing clearer, more interpretable results.

F. Immediate exporting of Round-2 results:
Added code to save and export kept_round2 and kicked_round2 subsets immediately after Round-2 analysis to spot check that 'kicked' and 'kept' analyses were implemented as intended.
```R
kept_round2 <- modifierAnalysis$kept_round2
kicked_round2 <- modifierAnalysis$kicked_round2

library(openxlsx)
exportDir <- file.path("results_RG_2025_05", "CyTOF", analysis, "kicked_kept")
dir.create(exportDir, recursive = TRUE, showWarnings = FALSE)

write.xlsx(kept_round2, file = file.path(exportDir, paste0(cellSubset, "_kept_round2.xlsx")), overwrite = TRUE)
write.xlsx(kicked_round2, file = file.path(exportDir, paste0(cellSubset, "_kicked_round2.xlsx")), overwrite = TRUE)
```
#### Flow Modifier Analysis Overhaul (Jun 2025)

* * **Flexible biological filters + looping**  
  * Added clear filters for **compartment**, **timePoint**, and **conditionLong** in both `nodes` and `design`.  
  * A triple loop now cycles through every combination of these factors (defaults: `compartment = c("global","antigen_specific")`, `timePoint = c("T1","T2","T3")`, `conditionLong = c("CEFX","DMSO","spike","MN")`).  
  * Each combo gets its own `analysis_suffix` (e.g. `antigen_specific_T2_MN`), its own result/plot folders, and a full modifier‑analysis run—no more hard‑coding to *global+T1+spike*.


* **Analysis‑suffix & directory hygiene**  
  * After loading data we now build `analysis_suffix` as  
    `<compartment>_<timePoint>_<conditionLong>` (e.g. `global_T1_spike`).  
  * All result/plot paths use this suffix, preventing accidental overwrites.

* **Full factorial looping**  
  * New triple loop over `compartment`, `timePoint`, `conditionLong` lets one call run the entire matrix of analyses—no more manual re‑filtering.

* **Robust plotting helpers**  
  * `makeBoxplots()` skips stat labels when no sig hits.  
  * `makeHeatmap()` exits cleanly (with a message) when no data remain and scales figure height to subset count for readability.

* **Quality‑of‑life tweaks**  
  * Directory creation moved *after* suffix definition.  
  * Added progress messages (`🚀`, `✅`, `⚠️`) for loop iterations.  
  * Modified Wilcoxon try‑catch already covered in CyTOF README applies here too.

* **New analysis script**  
  * **New analysis scripts**  
  * `abTitre_Flow_corr_analysis_RG_2025_06.R` – correlates **T3 antibody titres** with **T1 global & antigen-specific CD4RO⁺ flow frequencies** (DS & TC), auto‑generates stats tables and correlation plot .png files for significant correlations and writes them here: `results_RG_2025_05/vax_response/Tcell_corr_plots/`.
  * `abTitre_ALL_Ag_specific_flow_corr_Bcells_TC_and_DS.R` – same correlation framework but run on **B‑cell subsets** for both TC and DS cohorts; writes tables and correlation plots to this directory: `results_RG_2025_05/vax_response/Bcell_corr_plots/` folders.


> All code changes live in `analysis_FLOW_3_RG_2025_05.Rmd`, `analysis_FLOW_Bcell_RG_2025_05.Rmd` and (`Khor_covidvax_response_IMPACD/R`).

### Aging and Immune age
The next analysis is the aging and immune age analysis. This analysis has already been run for all datasets, but it too needs a rework. While the flow data analysis can be found in the old "gating_" files, the CyTOF data analysis can be found in the `analysis_CyTOF.Rmd` file. This file contains the DS vs TC comparisons as well as the aging and immune age analysis. This should all be reworked and put into the shared analysis files.

#### Analysis steps
1. Load the data
2. Build the model as a function of age
3. Test for differences
4. Run the modifier analysis
5. Generate plots and tables
6. Take the subsets that show aging in healthy
7. Build a model for immune age on those
8. Compare for differences in DS

#### Tables that are generated
- All stats, p-values, etc. for the aging effect before the modifier analysis
- Level 1 stats, basically the data from above subset to level 1 depth and fdr corrected for this sample size
- Modifier results: The modifier analysis should generate several tables, one that shows what survived after correction, one that shows all the comparisons that have been made and what the outcomes were.
- Immune age modeling, the coefficients and p-values for the immune age model
- Immune age comparison, the comparison of the immune age model between DS and TC

#### Plots that are generated
- Scatterplot of the aging effect, if there is a significant effect, draw the line of best fit and error bars
- Heatmap for the significant subsets that survived the modifier analysis
- Plots for the immune age comparison
- Volcano plot that show the subsets that were found the show different aging in DS and TC
  
#### Model
For now we are using a linear response model for the DS vs TC comparison. This is a simple model that is easy to interpret and can be used to generate the p-values for the differences between DS and TC. The model is defined as follows:
```R
fit = glm(pctParent ~ studyGroup*ageAtDraw, data = data, weights = data$weights)
```
For the reimplemenation of the flow data, we should probably use this same model, but loop over compartment and time points. Other models to consider would be a beta regression model or other ways to restrain the model to the 0-1 range.

#### Modifier analysis
For the comparison between DS and TC, the approach this time is slightly different to the one in DS vs TC. Instead of comparing the distributions between DS and TC, we now run a modifier analysis with in the age effects of DS and TC respectively
- start at level 1
- pick a terminal node, e.g. `A`
- find all the level 2 subsets that contain `A` as their terminal node, e.g. `A`, `B/A`, `C/A`
- loop over all children and build a model like `fit = glm(pctParent ~ ageAtDraw * subset, data = data_,  weights = weights)`
- contrast by comparting the subsets: `age_trends = emtrends(fit, ~ subset, var = "ageAtDraw")`
- if the trend is different between parent and child, report
  
**Problems:** We are definitely not filtering out enough subsets as we can see in the heatmaps. We are probably filtering out the level 2 ones, but then we are not filtering out the level 3 ones. This is a problem that needs to be addressed.

### Ab Titer response
In this analysis we want to model the Ab Titer response as a function of subsets. Unfortunately, we currently only have good ab Titer data for T3. T2 was truncated at 150, which makes it hard to model if all points have the same value. In the CyTOF data, where we have only T1, we can model subsets to find out which are going to predict a good response `abTiter_T3 vs pctParent_T1`. In the flow data, we have even more options: `abTiter_T3 vs pctParent_T1`, `abTiter_T3 vs pctParent_T2`, `abTiter_T3 vs pctParent_T3`, and then we also want to do this for our global and antigen specific compartments separately.This analysis is currently only partially implemented, see file `analysis_vax_response.Rmd`. 

> There is also a second file `analysis_vax_response.Rmd` that contains a lot of "general" analysis that is not specific to the IMPACD gating. It generates all kinds of plots that show the Ab Titer response and try to model that as a function of the frequency of antigen specifc cells, IgG+ cells, etc. 

`analysis_vax_response.Rmd` is currently in work and most of the functions are still in the file itself and have not been moved into `R/abTiter_analysis.R` yet. This is something that should be done in the future. The current status is that the vax response can be modeled for the CyTOF data, but the modifier analysis is not yet implemented. The flow data is not yet implemented at all.
#### Analysis steps
1. Load the data
2. Build the model as a function of pctParent
3. Test for differences
4. Run the modifier analysis
5. Generate plots and tables
6. Similar to immune age, we can now model the lack of response in DS as a function of the subsets that are predictive of a good response in TC. This is not yet fully implemented.


#### Tables that are generated
- All stats, p-values, etc. for the aging effect before the modifier analysis
- Level 1 stats, basically the data from above subset to level 1 depth and fdr corrected for this sample size
- Modifier results: The modifier analysis should generate several tables, one that shows what survived after correction, one that shows all the comparisons that have been made and what the outcomes were.
- ... AbTiter model and lack in DS (not yet implemented)


#### Plots that are generated
- Scatterplot of the ab titer effect, if there is a significant effect, draw the line of best fit and error bars
- Heatmap for the significant subsets that survived the modifier analysis
- Plots for the comparison in response between DS and TC
  
#### Model
For now we are using a linear response model for the DS vs TC comparison. This is a simple model that is easy to interpret and can be used to generate the p-values for the differences between DS and TC. The model is defined as follows:
```R
fit = glm(abTiter_T3 ~ studyGroup * pctParent, data = data, weights = data$weights)
```


#### Modifier analysis
Needs to be implemented - there is some code that was started in the analysis, but it is not yet complete. This is probably going to be very similar to how we do modifier analysis in the aging analysis.
  
## What has been done and next steps
So far, the only work that is truly complete (more or less) is the gating. The analysis is still in progress and probably needs a complete rework. The below table should indicate what kind of analysis has been done in the past and what is completely missing.

| Analysis | Flow Data | CyTOF Data | Files | Needs rework |
|----------|-----------|------------|------|------|
| Gating   | Done      | Done       |`gating_*.Rmd`| No
| DS vs TC | Done      | Done       |`analysis_CyTOF.Rmd`| Yes |
| Aging    | Done      | Done       |`analysis_CyTOF.Rmd`| Yes |
| Ab Titer | Not done  | Partially done |`analysis_vax_response.Rmd`| Yes |


Additionally, we have the following analyses that do no require the IMPACD gating data:
| Analysis | Status | File |
|----------|--------|------|
| Ab Titer response | Done | `analysis_vax_response_high_level.Rmd` | 
| Supervised Signatures | Done | `analysis_supervised_signatures.Rmd` |
| Custom requests | TBD | `custom_requests.Rmd` |

These analysis are more or less done for now, but might need to get updated based on requests from Bernard.

### Angle of attack
A good next step might be to start with rewriting and streamlining the analysis. Some initial steps have already been taken to homogenize analysis, for example by putting all the color into `R/project_settings.R` and highly shared functions (like fdr correction) into `R/common.R`. Any analysis file should benefit from that and could be structured like:
- load shared functions
- load functions specific to this analyis
- run on CyTOF data
- run on flow data

This would allow for a more streamlined analysis and easier comparison between the different data sets. 

The first step would be to start with the DS vs TC analysis and to completely build this out. Then move on to the aging and immune age analysis and the Ab Titer response analysis.

By keeping things as organized as possible, it should be relatively straight forward in the future to swap out smaller parts, like reworking the modifier analysis or changing the model.

### New work
There are some additional analysis that could be done:
- We should be receiving ab avidity data soon. This could be used to model the response to the vaccine in a different way.
- Comparing subsets that were found as significant in the aging and immune response modeling. But this could be tricky if we are not careful with the modifier analysis, as we might be comparing apples to oranges.
- Integrated analysis - combine subsets from the different rootpaths based on their significance. Pull them out of the individual `Rdata` counts. Then go looking for networks, build correlation matrices, etc. This could be a good way to get a more global view of the immune response to the vaccine.


### Other ideas
We came across a paper that does something very similar to our IMPACD approach: [terraFlow](https://www.cell.com/cell-reports/fulltext/S2211-1247(24)00955-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124724009550%3Fshowall%3Dtrue). I was wondering if we could use a network/tree like shown in Figure 1 B there, either to help us with the modifier analysis, or at least for visualization. It could be a cool "map" to build the entire network and then color the nodes by pvalues. Some initial network code can found in the `misc_code` folder.