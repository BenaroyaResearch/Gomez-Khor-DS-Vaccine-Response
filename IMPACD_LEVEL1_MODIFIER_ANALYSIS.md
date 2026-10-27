# IMPACD, Level 1, and modifier analysis workflow

This note summarizes the IMPACD-based analysis workflow adapted from the source repo:
`BenaroyaResearch/Khor_covidvax_response_IMPACD`.

It is intended to capture the core analysis logic needed for this repo's vaccine-response work, especially the concepts of:

- IMPACD gating and data structure
- Level 1 hits
- Modifier analysis
- Comparing parent vs child populations in a biologically consistent way

## 1) Repository context

The source IMPACD repo describes a generalized workflow:

1. Gating / IMPACD-based population extraction
2. Saving gated counts as R objects in `data/counts`
3. Running analysis on those subsets
4. Defining a primary analysis (e.g. DS vs TC, aging, Ab titer response)
5. Running a modifier analysis to reduce false discoveries

The project README states that the gated data is generated first, then the analysis is run on those `Rdata` objects. The main analysis files are root-level `.Rmd` files beginning with `analysis_` or `gating_`.

The relevant concepts from the source repo are summarized below.

## 2) Level 1 analysis

The source repo describes `Level 1 stats` as the data subset to level 1 depth and FDR-corrected for that sample size. It also calls out that level 1 hits are the entry point for the modifier analysis.

Conceptually:

- Fit a model for each subset/population at the relevant depth
- Define the subset of populations that are significant at the level-1 stage
- Adjust for multiple testing with FDR
- Keep the strongest or most interpretable signals for downstream modifier testing

A simple example from the source repo is:

```r
fit = glm(pctParent ~ studyGroup * ageAtDraw, data = data, weights = data$weights)
```

and for a simplified cohort-matched version:

```r
fit <- glm(pctParent ~ studyGroup, data = data, weights = data$weights)
```

This is a good template for a single-subset or repeated-subset pipeline when testing DS vs TC or age-related shifts.

## 3) Modifier analysis overview

The source repo notes that modifier analysis is necessary to reduce the false discovery rate. The core idea is to identify the main sources of variation and then test whether the observed shift is specific to a given population, rather than a parent-child artifact or gating-level artifact.

Important points from the source repo:

- If the difference in `A` and `B/A` is the same, but the difference in `C/A` is different, the relevant comparison is still `A` vs `C/A`.
- Modifier analysis is delicate and not fully standardized.
- The approach differs depending on the research question: DS vs TC, aging, or Ab titer response.

## 4) Parent/child matching logic

The source README specifically calls out that the original modifier logic was too loose. It matched child subsets only by the same terminal node and depth, without ensuring the modifier sequence matched the biological gating path.

The updated logic in the source repo filters children using the same terminal node and sequential modifier structure, e.g.:

```r
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

This matters because biological meaning depends on the gating order: for example, subsequent gates such as `modifier_1`, `modifier_2`, and `modifier_3` should line up with the parent path rather than being treated as interchangeable cell states.

## 5) Round 2 filtering for deeper modifier comparisons

The source repo also includes a second-pass filtering step for deeper nodes. Once level-2 results are known, level-3 subsets are retested against the relevant level-2 parent using the matching modifier relationship.

Example code from the source repo:

```r
if (max_depth >= 3) {
  depth2_parents <- dplyr::filter(stats_, depth == 2)
  depth3_round1  <- dplyr::filter(report, depth == 3)

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
}
```

The workflow is meant to distinguish subsets that genuinely survive after accounting for their parent or sibling relationships from those that simply reflect the same underlying signal at a different level.

## 6) Recommended working pattern for this repo

For this repo, a useful workflow is:

1. Start with level 1 hits from the relevant model
2. Keep only biologically interpretable terminal nodes / populations
3. For each level-1 hit, find deeper child subsets that share the same terminal node and modifier sequence
4. Recenter or compare parent vs child distributions (for example using a TC-centered background correction or a Wilcoxon test)
5. Retain only those deeper populations that remain significant after the modifier comparison
6. Export kept vs kicked subsets for manual review

## 7) Minimal R template

A practical template matching this logic is:

```r
library(dplyr)

perform_modifier_analysis <- function(stats_df, parent_population, fdr_cutoff = 0.05) {
  # Parent population is the level-1 or level-2 node we are testing against child subsets
  children <- stats_df %>%
    dplyr::filter(
      terminal_node == parent_population$terminal_node,
      depth == parent_population$depth + 1
    )

  if (parent_population$depth >= 2) {
    for (m in seq_len(parent_population$depth - 1)) {
      mod_col <- paste0("modifier_", m)
      children <- children %>%
        dplyr::filter(.data[[mod_col]] == parent_population[[mod_col]])
    }
  }

  out <- children %>%
    dplyr::mutate(
      p_value = NA_real_,
      modifier_vs_parent = parent_population$subset,
      kept = NA
    )

  for (i in seq_len(nrow(children))) {
    child <- children[i, ]
    x <- parent_population$value_vec
    y <- child$value_vec

    pval <- tryCatch(
      wilcox.test(x, y, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    out$p_value[i] <- pval
    out$kept[i] <- pval <= fdr_cutoff
  }

  out
}
```

This is not a production function; it is a conceptual adaptation of the source repo logic meant to keep the analysis traceable and biologically consistent.

## 8) Practical note

The source repo README emphasizes that the modifier analysis is still under active development and should be treated as a structured filtering step rather than a purely automated final answer. The main takeaways for this repo are:

- level 1 is the starting point
- modifier analysis should respect gating hierarchy and modifier order
- deeper nodes should be retested against the correct parent, not just any same-terminal-node subset
- the final interpretation should emphasize biological plausibility and consistency with the gating scheme

## 9) Source repository

Source:
`BenaroyaResearch/Khor_covidvax_response_IMPACD`

This note is a compact, study-specific adaptation of the README and analysis logic from that project.
