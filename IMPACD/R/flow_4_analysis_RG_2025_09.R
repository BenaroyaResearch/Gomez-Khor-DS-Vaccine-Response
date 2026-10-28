source("R/common.R")


runModel = function(pctParent, popCounts, weights, design,
                    family = "gaussian", includeSex = TRUE, cutOff = 10,
                    depth_lookup = NULL, level1_parent_cutoff = TRUE) {
  #'@title runModel
  #'@description Fit GLMs per subset after a median-count gate.
  #'  For depth==1 (level-1) and when level1_parent_cutoff=TRUE,
  #'  the gate uses the PARENT counts (weights) median by group (e.g., RBD+ / AIM+).
  #'  For depth>=2, the gate uses the SUBSET counts median by group (original behavior).
  #'@param pctParent   matrix [subset x sample]
  #'@param popCounts   matrix [subset x sample]
  #'@param weights     matrix [subset x sample] (parent counts)
  #'@param design      data.frame with sample metadata; must include studyGroup
  #'@param family      "gaussian" or "binomial" (beta not implemented)
  #'@param includeSex  kept for API compatibility
  #'@param cutOff      numeric; median threshold
  #'@param depth_lookup named integer vector: names=subset, values=depth
  #'@param level1_parent_cutoff logical; if TRUE, apply parent-median rule for level-1
  #'@return list of fitted model objects (one per kept subset)
  
  models <- list()
  zeroCounts <- character(0)
  
  message("Running linear model using ", family, " response for ", nrow(pctParent), " subsets")
  
  pb <- txtProgressBar(min = 0, max = nrow(pctParent), style = 3)
  i <- 0L
  
  for (subset_ in rownames(pctParent)) {
    i <- i + 1L
    setTxtProgressBar(pb, i)
    
    # Build per-sample row for this subset
    data <- data.frame(
      sampleName = colnames(pctParent),
      pctParent  = as.numeric(pctParent[subset_, ]),
      count      = as.numeric(popCounts[subset_, ]),
      weights    = as.numeric(weights[subset_, ]),
      stringsAsFactors = FALSE
    )
    
    # Join metadata
    data <- merge(data, design, by.x = "sampleName", by.y = "sampleName")
    
    # Determine depth for this subset (default FALSE if unknown)
    is_level1 <- FALSE
    if (!is.null(depth_lookup) && !is.null(depth_lookup[[subset_]])) {
      is_level1 <- (as.integer(depth_lookup[[subset_]]) == 1L)
    }
    
    # ---- Median gate ----
    if (level1_parent_cutoff && is_level1) {
      # Parent-level cutoff for level-1: use parent counts (weights)
      parent_median <- data %>%
        dplyr::group_by(studyGroup) %>%
        dplyr::summarise(medianParent = stats::median(weights, na.rm = TRUE), .groups = "drop")
      
      if (!any(is.finite(parent_median$medianParent) & parent_median$medianParent >= cutOff)) {
        zeroCounts <- c(zeroCounts, subset_)
        next
      }
    } else {
      # Original subset-level cutoff for non–level-1
      medianCount <- data %>%
        dplyr::group_by(studyGroup) %>%
        dplyr::summarise(medianCount = stats::median(count, na.rm = TRUE), .groups = "drop")
      
      if (!any(is.finite(medianCount$medianCount) & medianCount$medianCount >= cutOff)) {
        zeroCounts <- c(zeroCounts, subset_)
        next
      }
    }
    
    # Guard: skip if all pctParent are zero/NA after slice
    if (sum(data$pctParent, na.rm = TRUE) == 0) {
      zeroCounts <- c(zeroCounts, subset_)
      next
    }
    
    # ---- Fit model (same as your original) ----
    if (family == "gaussian") {
      data$weights[data$weights == 0] <- 1e-10
      fit <- glm(pctParent ~ studyGroup, data = data, weights = data$weights)
      
    } else if (family == "binomial") {
      data$weights <- as.integer(round(data$weights))
      data$pctParent <- data$pctParent / 100
      fit <- glm(pctParent ~ studyGroup*sex*ageAtDraw,
                 data = data, weights = data$weights,
                 family = binomial(link = "logit"))
      
    } else if (family == "beta") {
      next
    } else {
      next
    }
    
    models[[subset_]] <- fit
  }
  
  close(pb)
  return(models)
}

contrastDSvsTC = function(models, subsets, design, pctParent, popCounts, weights) {
  
  # Message
  message("Calculating contrasts for DS vs TC for ", length(models), " subsets")
  
  stats = c()
  pb <- txtProgressBar(min = 0, max = length(models), style = 3)
  i = 0
  emm_options("msg.interaction" = FALSE)
  
  for (subset_ in names(models)) {
    i = i + 1
    setTxtProgressBar(pb, i)
    
    emm <- emmeans(models[[subset_]], ~ studyGroup, cov.reduce = mean)
    coefs = as.data.frame(emmeans(models[[subset_]], ~ studyGroup, cov.reduce = mean))
    
    contrasts_ <- contrast(emm, method = "pairwise")
    contrasts_ = as.data.frame(contrasts_)
    
    contrasts_$foldChange = NA
    for (j in 1:nrow(contrasts_)) {
      contrasts_$foldChange[j] = dplyr::filter(coefs, studyGroup == "DS")$emmean / dplyr::filter(coefs, studyGroup == "TC")$emmean
    }
    
    contrasts_$subset = subset_
    
    stats = rbind(stats, contrasts_)
  }
  # close the progress bar
  close(pb)
  
  emm_options("msg.interaction" = TRUE)
  
  stats = merge(stats, subsets, by = "subset")
  
  stats$comparison = "DS vs TC"
  comparisons = unique(stats$comparison)
  
  
  message("Calculating summary statistics")
  summaryStats = c()
  
  pb <- txtProgressBar(min = 0, max = length(comparisons), style = 3)
  i = 0
  
  for (comparison_ in comparisons) {
    setTxtProgressBar(pb, i)
    control_group_ = "TC"
    study_group_ = "DS"
    
    stats_ = stats %>% 
      dplyr::filter(comparison == comparison_) %>%
      dplyr::select(comparison, subset)
    
    
    control_samples = dplyr::filter(design, studyGroup == control_group_)$sampleName
    study_samples = dplyr::filter(design, studyGroup == study_group_)$sampleName
    
    counts_control = popCounts[stats_$subset, control_samples]
    counts_study = popCounts[stats_$subset , study_samples]
    pctParent_control = pctParent[stats_$subset, control_samples]
    pctParent_study = pctParent[stats_$subset, study_samples]
    
    weights_control = weights[stats_$subset, control_samples]
    weights_study = weights[stats_$subset, study_samples]
    
    # calculate the median and sd of the counts
    # we have to make sure that we do this based on only the samples that had a 
    # nonzero count in the parent population (which is what we use the weights for)
    for (iSubset in 1:nrow(stats_) ) {
      controls_idx = weights_control[iSubset, ] > 0
      study_idx = weights_study[iSubset, ] > 0
      
      counts_control_subset = as.numeric(counts_control[iSubset, controls_idx])
      counts_study_subset = as.numeric(counts_study[iSubset, study_idx])
      
      pctParent_control_subset = as.numeric(pctParent_control[iSubset, controls_idx])
      pctParent_study_subset = as.numeric(pctParent_study[iSubset, study_idx])
      
      stats_$medianCountsControl[iSubset] = median(counts_control_subset)
      stats_$meanCountsControl[iSubset] =   mean(counts_control_subset)
      stats_$sdCountsControl[iSubset] =     sd(counts_control_subset)
      stats_$medianCountsStudy[iSubset] =   median(counts_study_subset)
      stats_$meanCountsStudy[iSubset] =     mean(counts_study_subset)
      stats_$sdCountsStudy[iSubset] =       sd(counts_study_subset)
      
      stats_$medianPctParentControl[iSubset] = median(pctParent_control_subset)
      stats_$meanPctParentControl[iSubset] =   mean(pctParent_control_subset)
      stats_$sdPctParentControl[iSubset] =     sd(pctParent_control_subset)
      stats_$medianPctParentStudy[iSubset] =   median(pctParent_study_subset)
      stats_$meanPctParentStudy[iSubset] =     mean(pctParent_study_subset)
      stats_$sdPctParentStudy[iSubset] =       sd(pctParent_study_subset)
      
    }
    
    #merge into stats
    summaryStats = rbind(summaryStats, stats_)
    
    i = i + 1
  }
  # close the progress bar
  setTxtProgressBar(pb, i)
  close(pb)
  
  stats = merge(stats, summaryStats, by = c("comparison", "subset"))
  
  stats = multipleTestCorrection(stats)
  
  return(stats)
}





modifierAnalysisDSvsTC = function(stats, design, pctParent, weights, fdrCutoff = 0.05, pCutoff = 0.05) {
  
  fullReport = c()
  removed_from_report = c()
  fullModifierComparisons = c()
  
  comparisons = dplyr::select(stats, comparison) %>% distinct()
  
  verbose = T
  
  control_group_ = "TC"
  study_group_ = "DS"
  
  for (comparison_ in comparisons$comparison) {
    
    
    stats_ = dplyr::filter(stats, 
                           comparison == comparison_, 
                           (medianCountsControl >= 10 |  medianCountsStudy >= 10),
                           fdr < fdrCutoff)
    
    if (nrow(stats_) == 0) {
      message(paste0("No significant differences in ", comparison_, ". Skipping."))
      next
    }
    
    message(paste0("Starting modifier analysis of ", control_group_, " vs ", study_group_, "."))
    
    # background correct counts
    control_samples = dplyr::filter(design, studyGroup == control_group_)$sampleName
    study_samples = dplyr::filter(design, studyGroup == study_group_)$sampleName
    
    pctParent_control = pctParent[stats_$subset, control_samples]
    pctParent_study = pctParent[stats_$subset , study_samples]
    
    # subtract the median count from stats_$medianCountsControl from the counts
    # pctParent_control_corrected = pctParent_control - stats_$medianPctParentControl # don't need the control group since we are correcting by the median of the controls
    pctParent_study_corrected = pctParent_study - stats_$medianPctParentControl
    
    # ------------------ START OF DEBUGGING CODE ------------------
    # 
    # # debug the CD4Teff/CCR6+ subset
    # if ("CD4Teff/CCR6+" %in% stats_$subset) {
    #   ccr6_row_idx <- which(stats_$subset == "CD4Teff/CCR6+")
    #   
    #   print("🔎 Debugging CD4Teff/CCR6+ subset explicitly:")
    #   
    #   print("👉 Raw Control (TC) values for CD4Teff/CCR6+:")
    #   print(pctParent_control[ccr6_row_idx, ])
    #   
    #   print("👉 Raw Study (DS) values for CD4Teff/CCR6+:")
    #   print(pctParent_study[ccr6_row_idx, ])
    #   
    #   print("👉 Median of control (TC) group used for correction (CD4Teff/CCR6+):")
    #   print(stats_$medianPctParentControl[ccr6_row_idx])
    #   
    #   print("👉 Corrected Study (DS) values for CD4Teff/CCR6+:")
    #   print(pctParent_study_corrected[ccr6_row_idx, ])
    # }
    # 
    # # ------------------ END OF DEBUGGING CODE --------------------
    # # --------- Additional debugging for corrected TC group  ---------
    # if ("CD4Teff/CCR6+" %in% stats_$subset) {
    #   
    #   ccr6_row_idx <- which(stats_$subset == "CD4Teff/CCR6+")
    #   
    #   print("👉 Corrected Control (TC) values for CD4Teff/CCR6+: (should be around zero)")
    #   corrected_TC = pctParent_control[ccr6_row_idx, ] - stats_$medianPctParentControl[ccr6_row_idx]
    #   print(corrected_TC)
    #   
    #   # Wilcoxon test between DS and TC corrected values:
    #   # ensure numeric vectors passed to Wilcoxon test
    #   wilcox_result = tryCatch({
    #     wilcox.test(
    #       as.numeric(corrected_TC),
    #       as.numeric(pctParent_study_corrected[ccr6_row_idx, ]),
    #       exact = FALSE
    #     )
    #   }, warning = function(w) {
    #     print(paste("⚠️ Warning in  DS vs TC Wilcoxon:", w$message))
    #     return(NULL)
    #   }, error = function(e) {
    #     print(paste("❌ Error in eDSvs TC Wilcoxon:", e$message))
    #     return(NULL)
    #   })
    #   
    #   
    #   if (!is.null(wilcox_result)) {
    #     print(paste("✅ Explicit DS vs TC Wilcoxon test p-value for CD4Teff/CCR6+:", wilcox_result$p.value))
    #   }
    # }
    # -------------------------------------------------------------------------
    weights_control = weights[stats_$subset, control_samples]
    weights_study = weights[stats_$subset, study_samples]
    
    max_depth = max(stats_$depth)
    
    subpopulations = list()
    for (depth_ in 1:max_depth) {
      subpopulations[[depth_]] = dplyr::filter(stats_, depth == depth_)
    }
    
    initialized = F
    modifierComparisons = c()
    
    for (depth_ in 1:(max_depth-1) ) {
      nPops = nrow(subpopulations[[depth_]])
      
      if (nPops == 0) {next}
      
      if (!initialized) {
        report = subpopulations[[depth_]]
        report$modifierAnalysisAgainst = NA
        report$modifierAnalysisPvalue = NA
        initialized = T
      }
      
      if (verbose) message(paste0("Modifier analysis of ", nrow(subpopulations[[depth_ + 1]]), " subpopulations with fdr corrected p-value < ", fdrCutoff, " at gating depth ", (depth_+1), ".") )
      iPop = 0
      
      if (verbose) pb = txtProgressBar(min = iPop, max = nPops, style = 3)
      
      for (iPop in 1:nPops) {
        if (verbose) setTxtProgressBar(pb, iPop)
        
        population_ = subpopulations[[depth_]][iPop, ]
        
        background_corrected_x = as.numeric(pctParent_study_corrected[population_$subset, weights_study[population_$subset, ] > 0] )
        
        #### DEBUGGING statment 2025 05 20
        # if (population_$subset == "CD4Teff/CCR6+") {
        #   print(paste("✅ Parent subset (X) being analyzed:", population_$subset))
        #   print("Background-corrected X values:")
        #   print(background_corrected_x)
        # }
        
        # find children
        # --- NEW child-selection logic -------------------------------
        children <- stats_ %>%
          dplyr::filter(
            terminal_node == population_$terminal_node,
            depth         == population_$depth + 1
          )
        
        # ensure modifiers 1...(depth-1) match parent exactly
        if (population_$depth >= 2) {
          for (m in seq_len(population_$depth - 1)) {
            mod_col <- paste0("modifier_", m)
            children <- children %>%
              dplyr::filter(.data[[mod_col]] == population_[[mod_col]])
          }
        }
        # --------------------------------------------------------------
        if (nrow(children) == 0) {next}
        
        children$modifierAnalysisAgainst = NA
        children$modifierAnalysisPvalue = NA
        
        # parents = c()
        for (iChild_ in 1:nrow(children)) {
          
          
          subset_ = children[iChild_, "subset"]
          background_corrected_y = as.numeric(pctParent_study_corrected[subset_, weights_study[subset_, ] > 0] )
          
          ############## DEBUGGING statement 20250520
          
          if (subset_ == "CD4Teff/CCR6+") {
            print(paste("✅ Child subset (Y) being analyzed:", subset_))
            print("Background-corrected Y values:")
            print(background_corrected_y)
          }
          
          # if (F) {
          # # diagnostics plots
          # data = data.frame(pctParent = background_corrected_x, node = "leaf", subset = population_$subset, type = "corrected", studyGroup = "exp") %>%
          #   rbind(data.frame(pctParent = background_corrected_y, node = "parent", subset = parent_$subset, type = "corrected", studyGroup = "exp")) %>%
          #   rbind(data.frame(pctParent = as.numeric(pctParent_study[parent_$subset, ]), node = "parent", subset = parent_$subset, type = "raw", studyGroup = "exp")) %>%
          #   rbind(data.frame(pctParent = as.numeric(pctParent_study[population_$subset, ]), node = "leaf", subset = population_$subset, type = "raw", studyGroup = "exp")) %>%  
          #   rbind(data.frame(pctParent = as.numeric(pctParent_control[parent_$subset, ]), node = "parent", subset = parent_$subset, type = "raw", studyGroup = "ref")) %>%
          #   rbind(data.frame(pctParent = as.numeric(pctParent_control[population_$subset, ]), node = "leaf", subset = population_$subset, type = "raw", studyGroup = "ref"))
          # 
          # diagnosticplot = ggplot(data, aes(x = node, y = pctParent, color = studyGroup)) +
          #   geom_boxplot() +
          #   facet_wrap(~type) +
          #   geom_jitter(position=position_jitterdodge(dodge.width = 0.75)) +
          #   labs(title = paste0("Node: ", population_$subset, " vs Parent: ", parent_$subset, "\nEach distribution is centered correct by median of ref group (0)"))
          # 
          # ggsave(diagnosticplot, file = paste0("plots/Tcells/DS_vs_TC/modifierAnalysis/", contrast_,  "/", str_replace_all(population_$subset, "/", "-"), "_vs_", str_replace_all(parent_$subset, "/", "-"), ".png"), width = 8, height = 6, units = "in", dpi = 300)
          # }
          
          # Robust Wilcoxon calculation that handles errors or warnings gracefully
          ##### DEBUGGING for CD4Teff/CCR6+ included
          p_value_child_vs_parent <- tryCatch({
            wilcox.test(background_corrected_x, background_corrected_y, exact = FALSE)$p.value
          }, warning = function(w) {
            if (subset_ == "CD4Teff/CCR6+" || population_$subset == "CD4Teff/CCR6+") {
              print(paste("⚠️ Warning in Wilcoxon test (CD4Teff/CCR6+):", w$message))
            }
            NA_real_
          }, error = function(e) {
            if (subset_ == "CD4Teff/CCR6+" || population_$subset == "CD4Teff/CCR6+") {
              print(paste("❌ Error in Wilcoxon test (CD4Teff/CCR6+):", e$message))
            }
            NA_real_
          })
          
          ######DEBUGGING print computed wilcoxon p-value for CD4Teff/CCR6+
          if (subset_ == "CD4Teff/CCR6+" || population_$subset == "CD4Teff/CCR6+") {
            print(paste("🎯 Wilcoxon p-value computed for Child subset", subset_, 
                        "vs Parent subset", population_$subset, ":", p_value_child_vs_parent))
          }
          
          children$modifierAnalysisAgainst[iChild_] <- population_$subset
          children$modifierAnalysisPvalue[iChild_] <- p_value_child_vs_parent
          
          # Safely check if p-value is non-NA before comparison
          if (!is.na(p_value_child_vs_parent) && p_value_child_vs_parent <= fdrCutoff) {
            report <- rbind(report, children[iChild_, ])
          }
          
        }
        
        modifierComparisons_ = children %>% 
          dplyr::select(child = subset, comparison, 
                        parent = modifierAnalysisAgainst, pValueModifierAnalysis = modifierAnalysisPvalue, 
                        logFCChild = foldChange, fdrChild = fdr) %>% 
          dplyr::mutate(
            logFCChild = log2(logFCChild),
            logFCParent = log2(population_$foldChange),
            fdrParent = population_$fdr)
        
        modifierComparisons = rbind(modifierComparisons, modifierComparisons_)
        
      }
      
      # next, we need to add all the depth+1 populations that were significant, but not tested during the modifier analysis (since the leaf node was not present in this depth)
      report_ = dplyr::filter(subpopulations[[depth_+1]], !(subset %in% modifierComparisons$child) )
      if (nrow(report_) > 0) {
        report_$modifierAnalysisAgainst = NA
        report_$modifierAnalysisPvalue = NA
        report = rbind(report, report_)
      }
      
      if (verbose) close(pb)
    }
    
    ########### insert code to test 2nd round of modifier analysis 12 may 2025 ##########
   
   # <-- end of first‑round depth loop
  
  # >>> paste SECOND‑ROUND snippet here <<<
    # ===== ROUND 2 ANALYSIS FOR DEPTH-3 CHILDREN =============================
    if (max_depth >= 3) {
      
      depth2_parents <- dplyr::filter(stats_, depth == 2)
      depth3_round1  <- dplyr::filter(report, depth == 3)
      
      kept_round2   <- data.frame()
      kicked_round2 <- data.frame()
      
      if (nrow(depth3_round1) > 0) {
        
        message("Starting Round-2 modifier analysis for depth-3 subsets")
        pb2 <- txtProgressBar(min = 0, max = nrow(depth3_round1), style = 3)
        
        for (iChild in seq_len(nrow(depth3_round1))) {
          setTxtProgressBar(pb2, iChild)
          
          child <- depth3_round1[iChild, ]
        
        # Initialize both columns 
        child$modifierAnalysisPvalue_round2  <- NA_real_
        child$modifierAnalysisAgainst_round2 <- NA_character_
        
        candidate_parents <- depth2_parents %>%
          dplyr::filter(
            terminal_node == child$terminal_node,
            modifier_1    == child$modifier_2
          )
        
        best_p <- Inf
        best_parent <- "No_valid_parent"
        
        if (nrow(candidate_parents) > 0) {
          for (ip in seq_len(nrow(candidate_parents))) {
            parent <- candidate_parents[ip, ]
            
            x <- as.numeric(pctParent_study_corrected[parent$subset,
                                                      weights_study[parent$subset, ] > 0])
            y <- as.numeric(pctParent_study_corrected[child$subset,
                                                      weights_study[child$subset, ] > 0])
            
            pval <- tryCatch({
              wilcox.test(x, y, exact = FALSE)$p.value
            }, warning = function(w) {
              print(paste("⚠️ Wilcoxon warning (NA returned): comparing child subset", child$subset, "with parent subset", parent$subset, "|", w$message))
              NA_real_
            }, error = function(e) {
              print(paste("❌ Wilcoxon error (NA returned): comparing child subset", child$subset, "with parent subset", parent$subset, "|", e$message))
              NA_real_
            })
            
            # # HLA-DR+/HLA-DR+ issue:  print subset comparisons that resulted in an NA:
            # if (is.na(pval)) {
            #   print(paste("🚨 NA generated in Wilcoxon test for Child subset", child$subset, "vs Parent subset", parent$subset))
            # }
            
            # Only update if pval is non-NA and smaller
            if (!is.na(pval) && pval < best_p) {
              best_p <- pval
              best_parent <- parent$subset
            }
            
          }
        } else {
          best_p <- NA_real_
          best_parent <- "No_valid_parent"
        }
        
        
        child$modifierAnalysisPvalue_round2 <- best_p
        child$modifierAnalysisAgainst_round2 <- best_parent
        
        if (!is.na(best_p) && best_p <= fdrCutoff) {
          kept_round2 <- rbind(kept_round2, child)
        } else {
          kicked_round2 <- rbind(kicked_round2, child)
        }
      }
      close(pb2)
      } else {
        message("⚠️ No subsets available for Round-2 modifier analysis at depth-3.")
      } # <-- THIS BRACE CLOSES THE if/else STATEMENT ABOVE
      
      
      
      
      report <- dplyr::bind_rows(
        dplyr::filter(report, depth < 3),
        kept_round2
      )
      
      removed_from_report <- rbind(removed_from_report, kicked_round2)
    }# <-- THIS final brace is necessary (closes your overall block/function/loop)
    # ===========================================================================
    
    
  
  
  
  
    report = distinct(report)
    
    fullReport = rbind(fullReport, report)
    
    fullModifierComparisons = rbind(fullModifierComparisons, modifierComparisons)
    
    
    if (!is.null(removed_from_report)) {
      removed_from_report = distinct(removed_from_report)
    }
    
  }
  
  condensedReport <- tryCatch({
    
    fullReport %>%
      dplyr::select(
        comparison, subset, terminal_node, modifier_1, modifier_2,
        foldChange, p.value, fdr, path,
        medianCountsControl, meanCountsControl, sdCountsControl,
        medianCountsStudy, meanCountsStudy, sdCountsStudy,
        medianPctParentControl, meanPctParentControl, sdPctParentControl,
        medianPctParentStudy, meanPctParentStudy, sdPctParentStudy,
        modifierAnalysisAgainst_round2,           
        modifierAnalysisPvalue_round2             
      ) %>%
      distinct() %>%
      arrange(fdr)
    
  }, error = function(e) {
    
    # Helpful message indicating what happened:
    message("⚠️ Round-2 modifier analysis columns missing: ", e$message)
    
    # Safely select without missing columns:
    fullReport %>%
      dplyr::select(
        comparison, subset, terminal_node, modifier_1, modifier_2,
        foldChange, p.value, fdr, path,
        medianCountsControl, meanCountsControl, sdCountsControl,
        medianCountsStudy, meanCountsStudy, sdCountsStudy,
        medianPctParentControl, meanPctParentControl, sdPctParentControl,
        medianPctParentStudy, meanPctParentStudy, sdPctParentStudy
        # Explicitly exclude round-2 columns as they're missing.
      ) %>%
      distinct() %>%
      arrange(fdr)
    
  })
  
  modifierHits <- condensedReport
  modifierResults <- fullReport
  modifierComparisons <- fullModifierComparisons
  
  # Return robustly packaged list
  return(list(
    modifierHits = modifierHits,
    modifierResults = modifierResults,
    modifierComparisons = modifierComparisons,
    kept_round2 = kept_round2,
    kicked_round2 = kicked_round2
  ))
  
} # END OF FUNCTION DEFINITION (Nothing goes below this curly bracket!)


run_DS_vs_TC = function(cellSubset) {
  
  # Defaults ( retained from your original function)
  fdrCutoff = 0.05
  pCutoff = 0.05 # for modifier analysis
  medianCellCutoff = 10
  analysis = "DS vs TC"
  
  # Set parameters  for loops
  compartments <- c("antigen_specific", "global")
  timePoints   <- c("T1", "T2", "T3")
  conditions   <- c("DMSO", "CEFX", "spike", "MN")
  
  # Outer loop: compartments
  for (compartment_label in compartments) {
    # Middle loop: timePoints
    for (timepoint_label in timePoints) {
      # Inner loop: conditions
      for (condition_label in conditions) {
        
        # Generate a clear, suffix for filenames/directories
        analysis_suffix <- paste(compartment_label, timepoint_label, condition_label, sep = "_")
        
        message("🚀 Starting analysis for: ", cellSubset, " | ", analysis_suffix)
        
        
        # Load and prep data ---------------------------------------------------
        filePath <- paste0("data/counts/flow/", cellSubset, ".Rdata")
        load(filePath)
        
        pctParent    <- matrices[["pctParent"]]
        popCounts    <- matrices[["popCounts"]]
        parentCounts <- matrices[["parentCounts"]]
        
        nodes        <- matrices[["nodes"]]
        design       <- matrices[["design"]]
        
        # Filter subsets by compartment (from nodes)
        compartment_subsets <- nodes$subset[nodes$compartment == compartment_label]
        
        pctParent    <- pctParent[rownames(pctParent) %in% compartment_subsets, ]
        popCounts    <- popCounts[rownames(popCounts) %in% compartment_subsets, ]
        parentCounts <- parentCounts[rownames(parentCounts) %in% compartment_subsets, ]
        nodes        <- nodes[nodes$subset %in% compartment_subsets, ]
        
        message("✅ Subsets after compartment filtering: ", nrow(nodes))
        
        # Prepare metadata 
        design$studyGroup <- factor(ifelse(design$patientType == "Control", "TC", "DS"), levels = c("TC", "DS"))
        design$sex <- factor(design$sex, levels = sort(unique(design$sex)))
        
        # Filter by your criteria
        design_filtered <- design %>%
          dplyr::filter(!control,
                        timePoint == timepoint_label,
                        conditionLong == condition_label)
        
        # Skip iteration  if no samples remain after filtering
        if (nrow(design_filtered) == 0) {
          message("⚠️ No samples after filtering for ", analysis_suffix, ". Skipping.")
          next
        }
        
        # Ensure matrices align with filtered design samples
        pctParent    <- pctParent[, design_filtered$sampleName]
        popCounts    <- popCounts[, design_filtered$sampleName]
        parentCounts <- parentCounts[, design_filtered$sampleName]
        
        # *** INSERT THIS: keep parentCounts keyed exactly like pctParent ***
        if (is.null(rownames(parentCounts)) || !identical(rownames(parentCounts), rownames(pctParent))) {
          rownames(parentCounts) <- rownames(pctParent)
        }
        
        weights <- parentCounts
        
        # Set up directories  ----------------------------------------
        tablesDir <- file.path("results_RG_2025_05", "Flow", analysis, cellSubset, analysis_suffix)
        
        statsDir      <- file.path(tablesDir, "stats")
        level1Dir     <- file.path(tablesDir, "level1")
        modifierDir   <- file.path(tablesDir, "modifierAnalysis")
        
        plotDir   <- file.path("plots_RG_2025_05", "Flow", analysis, cellSubset, analysis_suffix)
        level1Plots   <- file.path(plotDir, "level1")
        modifierPlots <- file.path(plotDir, "modifierAnalysis")
        
        dirs_to_create <- list(statsDir, level1Dir, modifierDir, level1Plots, modifierPlots)
        lapply(dirs_to_create, dir.create, recursive = TRUE, showWarnings = FALSE)
        
        # ---------- Run linear model (with parent-median gate for level-1) ----------
        #subsets <- nodes  # used downstream by contrastDSvsTC()
        
        # Make a subsets table keyed to the matrix rownames so merge() won't drop rows
        subsets <- dplyr::left_join(
          data.frame(subset = rownames(pctParent), stringsAsFactors = FALSE),
          nodes,
          by = "subset"
        )
        
        
        # *** REPLACE your existing depth_lookup line with THIS ***
        depth_lookup <- data.frame(subset_mat = rownames(pctParent), stringsAsFactors = FALSE) |>
          dplyr::left_join(nodes |> dplyr::select(subset, depth), by = c("subset_mat" = "subset")) |>
          {\(.) setNames(.$depth, .$subset_mat)}()
        
        
        # --- DIAGNOSTIC: end ---------------------------------------------------------------------
        # ---- Guard: nothing to model ----
        if (is.null(pctParent) || nrow(pctParent) == 0) {
          message("⚠️ 0 subsets to model for ", analysis_suffix, " — skipping.")
          next
        }
        
        models <- runModel(
          pctParent, popCounts, weights, design_filtered,
          family = "gaussian",
          cutOff = medianCellCutoff,
          depth_lookup = depth_lookup,        # uses matrix rownames as keys
          level1_parent_cutoff = TRUE         # parent-median gate for level-1
        )
        
        if (length(models) == 0) {
          message("⚠️ No significant differences for ", analysis_suffix)
          next
        }
        # ---------------------------------------------------------------------------
        

        # --- POST-MODEL DIAGNOSTIC: which level-1 subsets got dropped? ---------------------------------------------------------------------
        lvl1_expected <- nodes$subset[nodes$depth == 1]
        lvl1_kept     <- intersect(lvl1_expected, names(models))
        lvl1_missing  <- setdiff(lvl1_expected, names(models))
        
        message(sprintf("Level-1 expected: %d | kept: %d | missing: %d",
                        length(lvl1_expected), length(lvl1_kept), length(lvl1_missing)))
        
        if (length(lvl1_missing) > 0) {
          # Compute the two gates that can drop a subset inside runModel()
          tc_samples <- design_filtered$sampleName[design_filtered$studyGroup == "TC"]
          ds_samples <- design_filtered$sampleName[design_filtered$studyGroup == "DS"]
          
          med_by_grp <- function(row_obj, samp) {
            if (is.matrix(row_obj)) {
              vals <- suppressWarnings(as.numeric(row_obj[, samp, drop = TRUE]))
            } else {
              idx  <- match(samp, names(row_obj))
              vals <- suppressWarnings(as.numeric(row_obj[idx]))
            }
            stats::median(vals, na.rm = TRUE)
          }
          
          # parentCounts must key exactly like pctParent for the parent-median gate to work
          if (is.null(rownames(parentCounts)) || !identical(rownames(parentCounts), rownames(pctParent))) {
            rownames(parentCounts) <- rownames(pctParent)
          }
          
          diag_tbl <- data.frame(subset = lvl1_missing, stringsAsFactors = FALSE)
          diag_tbl$parent_med_TC <- vapply(diag_tbl$subset,
                                           function(s) med_by_grp(parentCounts[s, , drop = FALSE], tc_samples),
                                           numeric(1))
          diag_tbl$parent_med_DS <- vapply(diag_tbl$subset,
                                           function(s) med_by_grp(parentCounts[s, , drop = FALSE], ds_samples),
                                           numeric(1))
          diag_tbl$parent_rule_pass <- with(diag_tbl,
                                            (is.finite(parent_med_TC) & parent_med_TC >= medianCellCutoff) |
                                              (is.finite(parent_med_DS) & parent_med_DS >= medianCellCutoff))
          
          diag_tbl$pct_sum <- vapply(diag_tbl$subset,
                                     function(s) sum(as.numeric(pctParent[s, , drop = TRUE]), na.rm = TRUE),
                                     numeric(1))
          
          message("Level-1 missing details:")
          print(diag_tbl, row.names = FALSE)
        }
        # --- END POST-MODEL DIAGNOSTIC ------------------------------------------------------------------
        
       
        # Create contrasts ------------------------------------------
        stats = contrastDSvsTC(models, subsets, design_filtered, pctParent, popCounts, weights)
        
        level1Results = getLevel1(stats)
        
        # Save results  ----------------------------------------------
        saveToSpreadsheet(stats, file.path(statsDir, paste0(cellSubset, "_", analysis_suffix, "_stats.xlsx")))
        saveToSpreadsheet(level1Results, file.path(level1Dir, paste0(cellSubset, "_", analysis_suffix, "_level1.xlsx")))
        
        # Modifier analysis  -----------------------------------------
        if (nrow(dplyr::filter(stats, fdr < fdrCutoff)) == 0) {
          message("⚠️ No significant modifier results for ", analysis_suffix)
          next
        }
        
        modifierAnalysis <- modifierAnalysisDSvsTC(stats, design_filtered, pctParent, weights, 
                                                   fdrCutoff = fdrCutoff, pCutoff = pCutoff)
        
        modifierHits <- modifierAnalysis$modifierHits
        modifierResults <- modifierAnalysis$modifierResults
        modifierComparisons <- modifierAnalysis$modifierComparisons
        
        kept_round2 <- modifierAnalysis$kept_round2
        kicked_round2 <- modifierAnalysis$kicked_round2
        
        #include cellSubset to avoid redundant exports:
        exportDir <- file.path("results_RG_2025_05", "Flow", analysis, cellSubset, analysis_suffix, "kicked_kept")
        dir.create(exportDir, recursive = TRUE, showWarnings = FALSE)
        
        write.xlsx(kept_round2, file.path(exportDir, paste0(cellSubset, "_", analysis_suffix, "_kept_round2.xlsx")), overwrite = TRUE)
        write.xlsx(kicked_round2, file.path(exportDir, paste0(cellSubset, "_", analysis_suffix, "_kicked_round2.xlsx")), overwrite = TRUE)
        
        
        saveModifierResultsToSpreadsheet(modifierHits, modifierResults, modifierComparisons, 
                                         file.path(modifierDir, paste0(cellSubset, "_", analysis_suffix, "_modifierAnalysis.xlsx")))
        
        # Save results  ----------------------------------------------
        results <- list(models = models, stats = stats, modifierHits = modifierHits, 
                        modifierResults = modifierResults, modifierComparisons = modifierComparisons)
        
        # Correct Rdata export
        dir.create(file.path("analysis_results", analysis, cellSubset, analysis_suffix), recursive = TRUE, showWarnings = FALSE)
        
        save(results, 
             file = file.path("analysis_results", analysis, cellSubset, analysis_suffix, paste0(cellSubset, "_", analysis_suffix, ".Rdata")))
        
        
        
        # Make heatmaps and boxplots  --------------------------------
        hmtitle <- paste0("Sig. level 1 of ", cellSubset, " (", analysis_suffix, ")")
        filePath <- file.path(level1Plots, paste0("hm_", cellSubset, "_", analysis_suffix, ".png"))
        makeHeatmap(level1Results, pctParent, design_filtered, colors, hmtitle, filePath = filePath)
        
        # Boxplots 
        boxPath <- file.path(level1Plots, paste0("box_", cellSubset, "_", analysis_suffix, ".png"))
        makeBoxplots(level1Results, pctParent, weights, design_filtered, colors, filePath = boxPath, title = hmtitle, add_stats = FALSE)
        
        boxSigPath <- file.path(level1Plots, paste0("box_", cellSubset, "_", analysis_suffix, "_with_significance.png"))
        makeBoxplots(level1Results, pctParent, weights, design_filtered, colors, filePath = boxSigPath, title = paste0(hmtitle, " with significance indicator"), add_stats = TRUE)
        
        # Modifier analysis heatmap ---------------------------------
        hmtitle_mod <- paste0("Modifier analysis results of ", cellSubset, " (", analysis_suffix, ")")
        filePath_mod <- file.path(modifierPlots, paste0("hm_", cellSubset, "_", analysis_suffix, ".png"))
        makeHeatmap(modifierResults, pctParent, design_filtered, colors, hmtitle_mod, filePath = filePath_mod)
        
        message("✅ Finished analysis for: ", cellSubset, " | ", analysis_suffix)
        
        
      } # conditionLong loop
    } # timePoint loop
  } # compartment loop
  
  message("🎉 All analyses completed successfully!")
}

  
  
  # Return ---------------------------------------------------------------------
  #return(NULL)
  


makeHeatmap = function(resultsTable, pctParent, design, colors, title = "heatmap", filePath = NULL) {
  
  # Prep data for heatmap
  top.table = dplyr::filter(resultsTable, fdr <= 0.05)
  
  if (nrow(top.table) == 0) {
    message("No significant results found")
    return(NULL)
  }
  
  study_group = "DS"
  control_group =  "TC"
  
  top.table = dplyr::select(top.table, subset, foldChange, p.value, fdr) %>%
    distinct()
  
  top.table$logFC = log2(top.table$foldChange)
  
  # Extract meta data and reorder
  design_ = dplyr::filter(design, studyGroup %in% c(study_group, control_group))
  
  reorder = order(design_$studyGroup, design_$ageAtDraw, design_$sex)
  design_ = design_[reorder, ]
  
  # Split heatmap by studyGroup
  colSplit = design_$studyGroup
  
  # Subset and reorder the data matrix
  pctParent_ = pctParent[match(top.table$subset, rownames(pctParent)), 
                         match(design_$sampleName, colnames(pctParent))]
  
  # Get the labels
  labels = row.names(pctParent_)
  
  
  # Scale data and prep for plotting
  if (nrow(pctParent_) == 1) {
    tp = scale(pctParent_)
  } else {
    tp = t(scale(t(pctParent_)))
    # remove any rows with all NA
  }
  
  # Remove rows with all NA, Using drop = FALSE ensures tp remains a matrix even with only one row.
  tp <- tp[!apply(is.na(tp), 1, all), , drop = FALSE]
  
  # Ensure tp isn't empty or NA-filled
  if (is.null(tp) || all(is.na(tp)) || length(tp) == 0 || is.na(nrow(tp)) || nrow(tp) == 0) {
    message("⚠️ No usable data available for heatmap: ", title)
    return(NULL)
  }
  
  
  
  # Decide on wether to add labels and row height
  showRowNames = length(labels) < 30
  if (nrow(top.table) < 100) {
    heightPerSubset = 0.2
  } else {
    heightPerSubset = 0.02
  }
  
  # Generate annotation
  ha <- HeatmapAnnotation(df = list(patientType = design_$patientType,
                                    sex = design_$sex,
                                    age = design_$ageAtDraw), 
                          col = list(patientType = colors$studyGroup,
                                     sex = colors$sex,
                                     age = colors$ageAtDraw)
  )
  
  # Build heatmap
  hm = Heatmap(tp, 
               show_column_names = FALSE, 
               show_row_names = showRowNames, 
               cluster_columns = FALSE, 
               name = "z scaled pct parent", # legend 
               top_annotation = ha, 
               col = colors$heatmap, 
               cluster_rows = TRUE, 
               column_title = title, # title
               show_row_dend = FALSE,
               na_col = "grey",
               column_split = colSplit)
  
  if (is.null(filePath)) {
    plot(hm)
  } else {
    subset_count <- nrow(top.table)
    height_ <- max(5, 1 + 0.25 * subset_count)
    
    
    png(filePath, height= height_, width = 9, units = "in", res = 600)
    print(hm)
    invisible(dev.off())
  }
  
  return(hm)
}



makeBoxplots = function(stats, pctParent, weights, design, colors, filePath = NULL, title = "Box plot", significant_only = F, add_stats = F) {
  bar_width <- 0.5
  
  # subset to only significant if requested 
  if (significant_only) {
    stats = dplyr::filter(stats, fdr < 0.05)
  }
  
  # subset the data
  data = pctParent[match(unique(stats$subset), rownames(pctParent)), ]
  
  # convert data to long format, use column names as sample names, and row names as subset names
  data = as.data.frame(data)
  data$subset = rownames(data)
  data = tidyr::gather(data, key = "sampleName", value = "pctParent", -subset)
  
  # merge the data with the design matrix
  data = merge(data, design, by = "sampleName")
  
  
  # get the weights for each row from the weights matrix (rows are subset and columns are sampleName)
  data$weights = NA
  for (subset_ in unique(data$subset)) {
    data$weights[data$subset == subset_] = unlist(weights[subset_, match(data$sampleName[data$subset == subset_], colnames(weights))])
  }
  
  # hist(data$weights)
  # keep only data with weights > 0
  data = dplyr::filter(data, weights > 0)
  
  # Define the dodge position for the scatter
  dodge = position_dodge(width = 0.75)
  
  ggboxplot = ggplot(data) + 
    # most bottom layer is boxplot
    geom_boxplot(aes(fill = studyGroup, x = studyGroup, y = pctParent), color = "black", outlier.shape = NA, alpha = 0.5, position = dodge) + 
    # add data points
    geom_jitter(aes(fill = studyGroup, x = studyGroup, y = pctParent), shape = 21, color = "black", size = 2, position=position_jitterdodge(dodge.width = 0.75), alpha = 0.8) + 
    # additional settings
    facet_wrap(subset~. , scales = "free") +
    theme_bw() +
    scale_fill_manual(values = c(colors$studyGroupShort, colors$studyGroupFit))  +
    labs(x = "studyGroup", y = "percent parent", fill = "study group", color = "fit", title = title) +
    theme(strip.text.x = element_text(size = 7)) 
  
  if (add_stats) {
    dataSummary <- data %>%
      dplyr::filter(weights > 0) %>%
      dplyr::group_by(subset) %>%
      dplyr::summarize(group_median_pct_parent = median(pctParent), 
                       group_sd_pct_parent = sd(pctParent))
    
    dataSummary <- merge(stats, dataSummary, by = "subset")
    
    dataSummary$signif_position <- dataSummary$group_median_pct_parent + 2 * dataSummary$group_sd_pct_parent
    dataSummary$signif_symbol <- dataSummary$signif_position + dataSummary$group_sd_pct_parent / 2
    dataSummary$label <- if_else(dataSummary$fdr < 0.005, "**", if_else(dataSummary$fdr < 0.05, "*", ""))
    
    dataSummary <- dplyr::filter(dataSummary, label != "")
    
    # ONLY difference: Safely handle no significant hits
    if (nrow(dataSummary) > 0) {
      dataSummary$studyGroup <- 1.5
      
      # add stats to plot
      ggboxplot <- ggboxplot + 
        geom_text(data = dataSummary, aes(x = studyGroup, y = signif_symbol, label = label), size = 6, vjust = 0.5) +
        geom_segment(data = dataSummary, aes(x = studyGroup - bar_width,
                                             xend = studyGroup + bar_width, y = signif_position, yend = signif_position), size = 0.8) +
        scale_y_continuous(expand = expansion(mult = c(0.1, .2)))
    } else {
      message("⚠️ No significant results to annotate on boxplot for ", title)
    }
  }
  
  
  
  if (is.null(filePath)) {
    print(ggboxplot)
  } else {
    ggsave(filePath, plot = ggboxplot, width = 14, height = 11, units = "in", dpi = 300)
  }
  
  return(ggboxplot)
}


ageEffects = function(models, subsets, design, pctParent, popCounts, weights) {
  
  # Message
  message("Testing for age dependencies in ", length(models), " subsets")
  
  stats = c()
  pb <- txtProgressBar(min = 0, max = length(models), style = 3)
  i = 0
  emm_options("msg.interaction" = FALSE)
  AICs = c()
  for (subset_ in names(models)) {
    i = i + 1
    setTxtProgressBar(pb, i)
    
    age_trends <- emtrends(models[[subset_]], 
                           ~ studyGroup, 
                           var = "ageAtDraw")
    
    age_trends = test(age_trends)
    
    age_trends$subset = subset_
    
    stats = rbind(stats, age_trends)
    
    AICs = c(AICs, AIC(models[[subset_]]))
    
    ### DONT FORGET TO ADD CORRELATION HERE
    #     suppressMessages(
    # tryCatch({
    #   cor_results_ = data_ %>%
    #         group_by(studyGroup, condition) %>%
    #         summarise(correlation = cor(pctToModel, abTiterT2, use = "complete.obs"), 
    #                   corr.p.value = cor.test(pctToModel, abTiterT2)$p.value)
    # }, error = function(e) {
    #   cor_results_ = data_ %>%
    #         group_by(studyGroup, condition) %>%
    #         summarise(correlation = NA, p_value = NA)
    # })
    # )
    # cor_results_$subset = subset_
    # cor_results_$TFHandRoot = TFHandRoot_
    # 
    # cor_results = rbind(cor_results, cor_results_)
    
  }
  # close the progress bar
  close(pb)
  
  emm_options("msg.interaction" = TRUE)
  
  names(AICs) = names(models)
  stats = merge(stats, subsets, by = "subset")
  
  
  stats$comparison = paste0(stats$studyGroup, " age effect")
  comparisons = unique(stats$comparison)
  
  
  message("Calculating summary statistics")
  summaryStats = c()
  
  pb <- txtProgressBar(min = 0, max = length(comparisons), style = 3)
  i = 0
  
  summaryStats = c()
  for (comparison_ in comparisons) {
    stats_ = stats %>% 
      dplyr::filter(comparison == comparison_) %>%
      dplyr::select(comparison, subset, studyGroup)
    
    studyGroup_ = unique(stats_$studyGroup)
    
    study_samples = dplyr::filter(design, studyGroup == studyGroup_)$sampleName
    
    counts_study = popCounts[stats_$subset , study_samples]
    pctParent_study = pctParent[stats_$subset, study_samples]
    weights_study = weights[stats_$subset, study_samples]
    
    # calculate the median and sd of the counts
    # we have to make sure that we do this based on only the samples that had a 
    # nonzero count in the parent population (which is what we use the weights for)
    for (iSubset in 1:nrow(stats_) ) {
      study_idx = weights_study[iSubset, ] > 0
      
      counts_study_subset = as.numeric(counts_study[iSubset, study_idx])
      pctParent_study_subset = as.numeric(pctParent_study[iSubset, study_idx])
      
      stats_$medianCounts[iSubset] =   median(counts_study_subset)
      stats_$meanCounts[iSubset] =     mean(counts_study_subset)
      stats_$sdCounts[iSubset] =       sd(counts_study_subset)
      
      stats_$medianPctParent[iSubset] = median(pctParent_study_subset)
      stats_$meanPctParent[iSubset] =   mean(pctParent_study_subset)
      stats_$sdPctParent[iSubset] =     sd(pctParent_study_subset)
      
    }
    
    #merge into stats
    summaryStats = rbind(summaryStats, stats_)
    i = i + 1
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  summaryStats = dplyr::select(summaryStats, -studyGroup)
  
  stats = merge(stats, summaryStats, by = c("comparison", "subset"))
  
  stats = multipleTestCorrection(stats)
  
  return(stats)
}

modifierAnalysisAgeEffect = function(stats, design, pctParent, weights, fdrCutoff = 0.05, pCutoff = 0.05, medianCountCutoff = 10) {
  
  # Setup
  fullReport = c()
  removed_from_report = c()
  fullModifierComparisons = c()
  
  verbose = T
  
  comparisons = dplyr::select(stats, comparison, studyGroup) %>% distinct()
  
  # Modifier analysis for each comparison --------------------------------
  for (comparison_ in comparisons$comparison) {
    
    studyGroup_ = comparisons$studyGroup[comparisons$comparison == comparison_]
    
    stats_ = dplyr::filter(stats, 
                           comparison == comparison_, 
                           (medianCounts >= medianCountCutoff),
                           fdr < fdrCutoff)
    
    if (nrow(stats_) == 0) {
      message(paste0("No significant differences in ", comparison_, ". Skipping."))
      next
    }
    
    message(paste0("Starting modifier analysis of age trends in ", studyGroup_, "."))
    
    design_ = dplyr::filter(design, 
                            studyGroup == studyGroup_)
    
    max_depth = max(stats_$depth)
    
    subpopulations = list()
    for (depth_ in 1:max_depth) {
      subpopulations[[depth_]] = dplyr::filter(stats_, depth == depth_)
    }
    
    initialized = F
    modifierComparisons = c()
    
    ## Step through different depths -------------------------------------------
    for (depth_ in 1:(max_depth-1) ) {
      nPops = nrow(subpopulations[[depth_]])
      
      if (nPops == 0) {next}
      
      if (!initialized) {
        # print(depth_)
        report = subpopulations[[depth_]]
        
        report$modifierAnalysisAgainst = NA
        report$modifierAnalysisDeltaEstimate = NA
        report$modifierAnalysisPvalue = NA
        report$modifierAnalysisSE = NA
        report$modifierAnalysisDf = NA
        report$modifierAnalysisT.ratio = NA
        
        initialized = T
      }
      
      if (verbose) message(paste0("Modifier analysis of ", nrow(subpopulations[[depth_ + 1]]), " subpopulations with fdr corrected p-value < ", fdrCutoff, " at gating depth ", (depth_+1), ".") )
      iPop = 0
      
      if (verbose) pb = txtProgressBar(min = iPop, max = nPops, style = 3)
      
      # Step through all pops at this level ------------------------------------
      for (iPop in 1:nPops) {
        if (verbose) setTxtProgressBar(pb, iPop)
        
        population_ = subpopulations[[depth_]][iPop, ]
        
        
        data_x = data.frame(sampleName = colnames(pctParent),
                            pctParent = as.numeric(pctParent[population_$subset,]),
                            weights = as.numeric(weights[population_$subset,])
        )
        data_x = dplyr::filter(data_x, sampleName %in% design_$sampleName)
        # data_x$weights[data_x$weights == 0] = 1e-10
        data_x = dplyr::filter(data_x, weights > 0)
        data_x$subset = "parent"
        
        # find children
        children = dplyr::filter(stats_, terminal_node == population_$terminal_node, depth == depth_ + 1)
        if (depth_ == 2) {
          children = dplyr::filter(children, (modifier_1 == population_$modifier_1 | modifier_2 == population_$modifier_1) ) 
        }
        if (nrow(children) == 0) {
          next
        }
        children$modifierAnalysisAgainst = NA
        children$modifierAnalysisDeltaEstimate = NA
        children$modifierAnalysisPvalue = NA
        children$modifierAnalysisSE = NA
        children$modifierAnalysisDf = NA
        children$modifierAnalysisT.ratio = NA
        
        
        if (sum(data_x$pctParent) == 0) {
          next
        }
        
        # Compare parent to children pops ------------------------------------
        for (iChild_ in 1:nrow(children)) {
          subset_ = children[iChild_, "subset"]
          
          data_y = data.frame(sampleName = colnames(pctParent),
                              pctParent = as.numeric(pctParent[subset_,]),
                              weights = as.numeric(weights[subset_,])
          )
          data_y = dplyr::filter(data_y, sampleName %in% design_$sampleName)
          data_y$weights[data_y$weights == 0] = 1e-10
          data_y$subset = "child"
          
          data_ = rbind(data_x, data_y)
          
          data_ = merge(data_, design_, by.x = "sampleName", by.y = "sampleName")
          
          fit = glm(pctParent ~ ageAtDraw * subset, data = data_,  weights = weights)
          
          age_trends = emtrends(fit, ~ subset, var = "ageAtDraw")
          
          contrast_results = contrast(age_trends, method = "pairwise")
          contrast_results = as.data.frame(contrast_results)
          
          p_value_child_vs_parent = contrast_results$p.value[1]
          
          children$modifierAnalysisAgainst[iChild_] = population_$subset
          children$modifierAnalysisDeltaEstimate[iChild_] = contrast_results$estimate[1]
          children$modifierAnalysisPvalue[iChild_] = p_value_child_vs_parent
          children$modifierAnalysisSE[iChild_] = contrast_results$SE[1]
          children$modifierAnalysisDf[iChild_] = contrast_results$df[1]
          children$modifierAnalysisT.ratio[iChild_] = contrast_results$t.ratio[1]
          
          if (p_value_child_vs_parent <= fdrCutoff) {
            report = rbind(report, children[iChild_, ])
          }
          
        }
        
        modifierComparisons_ = children %>% 
          dplyr::select(child = subset, comparison = comparison, parent = modifierAnalysisAgainst, 
                        pValueModifierAnalysis = modifierAnalysisPvalue, modifierAnalysisDeltaEstimate, 
                        modifierAnalysisSE, modifierAnalysisDf, modifierAnalysisT.ratio,
                        ageAtDraw.trend.child = ageAtDraw.trend, fdrChild = fdr)%>% 
          dplyr::mutate(
            ageAtDraw.trend.parent = population_$ageAtDraw.trend,
            fdrParent = population_$fdr
          )
        
        
        modifierComparisons = rbind(modifierComparisons, modifierComparisons_)
        
      }
      
      # Save to report  --------------------------------------------------------
      # next, we need to add all those populations that were significant, but not tested during the modifier analysis, since the leaf node was not present in this depth
      report_ = dplyr::filter(subpopulations[[depth_+1]], !(subset %in% modifierComparisons$child) )
      if (nrow(report_) > 0) {
        report_$modifierAnalysisAgainst = NA
        report_$modifierAnalysisDeltaEstimate = NA
        report_$modifierAnalysisPvalue = NA
        report_$modifierAnalysisSE = NA
        report_$modifierAnalysisDf = NA
        report_$modifierAnalysisT.ratio = NA
        
        report = rbind(report, report_)
      }
      
      if (verbose) close(pb) 
      
    }
    
    
    report <- distinct(report)
    
    fullReport <- rbind(fullReport, report)
    fullModifierComparisons <- rbind(fullModifierComparisons, modifierComparisons)
    
    if (!is.null(removed_from_report)) {
      removed_from_report <- distinct(removed_from_report)
    }
    
    # ensure these columns always exist:
    if (!"modifierAnalysisAgainst_round2" %in% colnames(fullReport)) {
      fullReport$modifierAnalysisAgainst_round2 <- NA_character_
    }
    
    if (!"modifierAnalysisPvalue_round2" %in% colnames(fullReport)) {
      fullReport$modifierAnalysisPvalue_round2 <- NA_real_
    }
    
  } # end of function here 
  
  # Cleanup reports and return -------------------------------------------------
  condensedReport = fullReport %>% 
    dplyr::select(comparison = comparison, subset, studyGroup,
                  terminal_node, modifier_1, modifier_2, 
                  ageAtDraw.trend, p.value, fdr, path, 
                  medianCounts,	meanCounts,	sdCounts,
                  medianPctParent,	meanPctParent,	sdPctParent
    ) %>% 
    dplyr::distinct() %>%
    dplyr::arrange(fdr)
  
  modifierHits = condensedReport
  modifierResults = fullReport
  modifierComparisons = fullModifierComparisons
  
  # return(list(modifierResults = modifierResults, modifierComparisons = modifierComparisons, removedFromReport = removed_from_report))
  return(list(modifierHits = modifierHits, modifierResults = modifierResults, modifierComparisons = modifierComparisons))
}


plotAgeTrends = function(data, design, pctParent, parentCounts, models, colors, fdrCutoff = 0.05, filePath = NULL, plotTitle = "Age trends plot", significant_only = T, add_trends = T) {
  
  # Set up for the plot --------------------------------------------------------  
  ageAtDraw = range(design$ageAtDraw)
  ageAtDraw_vec = seq(ageAtDraw[1], ageAtDraw[2], length.out = 50)
  
  if (significant_only) {
    data = dplyr::filter(data, fdr <= fdrCutoff)
  }
  
  # calculate plot dimensions
  nCols = ceiling(sqrt(length(unique(data$subset))))
  height_ = 2*nCols + 0.5
  width_ = 3*nCols
  
  # extract the data
  pctParent_ = pctParent[match(unique(data$subset), rownames(pctParent)), ]
  pctParent_ = as.data.frame(pctParent_)
  pctParent_$subset = rownames(pctParent_)
  pctParent_ = tidyr::gather(pctParent_, key = "sampleName", value = "pctParent", -subset)
  pctParent_ = merge(pctParent_, design, by.x = "sampleName", by.y = "sampleName")
  
  pctParent_$weights = NA
  
  # get the weights for each row from the weights matrix (rows are subset and columns are sampleName)
  for (subset_ in unique(pctParent_$subset)) {
    pctParent_$weights[pctParent_$subset == subset_] = unlist(parentCounts[subset_, match(pctParent_$sampleName[pctParent_$subset == subset_], colnames(parentCounts))])
  }
  
  pctParent_ = dplyr::filter(pctParent_, !is.na(weights), weights > 0)
  
  # Start the plot -------------------------------------------------------------
  ggscatterplot = ggplot(pctParent_) 
  
  ## Calculate and add trends and ribbons if requested -------------------------
  if (add_trends) {
    predictions = c()
    emm_options("msg.interaction" = FALSE)
    
    for (comparison_ in unique(data$comparison)) {
      
      data_ = dplyr::filter(data, comparison == comparison_)
      studyGroup_ = unique(data_$studyGroup)
      
      for (subset_ in unique(data_$subset)) {
        
        if (data_$fdr[data_$subset == subset_] > fdrCutoff) { next }
        
        predictions_ = emmeans(models[[subset_]], 
                               ~ studyGroup*ageAtDraw, cov.reduce = mean,
                               at = list(ageAtDraw = ageAtDraw_vec, 
                                         studyGroup = studyGroup_)
        )
        
        predictions_ = as.data.frame(predictions_)
        predictions_$subset = subset_
        predictions_ = dplyr::rename(predictions_, pctParent = emmean)
        predictions = rbind(predictions, predictions_)
      }
    }
    emm_options("msg.interaction" = TRUE)
    
    ggscatterplot = ggscatterplot + 
      # add the prediction data
      geom_ribbon(data = predictions, aes(x = ageAtDraw, ymin = lower.CL, ymax = upper.CL, fill = studyGroup), alpha = 0.2) +
      geom_line(data = predictions, aes(x = ageAtDraw, y = pctParent, color = studyGroup), linewidth = 1, linetype = 4)
  }
  
  # Finish and export plot ----------------------------------------------------
  ggscatterplot = ggscatterplot + 
    geom_point(aes(fill = studyGroup, x = ageAtDraw, y = pctParent), color = "black", shape = 21, alpha = 0.5) + 
    # additional settings
    facet_wrap(~ subset, scales = "free", ncol = nCols) +
    theme_bw() +
    scale_fill_manual(values = colors$studyGroupShort) + 
    # scale_fill_manual(values = c(colors$studyGroupShort, colors$studyGroupFit))  +
    scale_color_manual(values = colors$studyGroupShort) +
    labs(x = "age", y = "percent parent", fill = "study group",
         title = plotTitle ) +
    guides(color = "none") +
    theme(strip.text.x = element_text(size = 7))
  
  if (is.null(filePath)) {
    print(ggscatterplot)
  } else {
    ggsave(filePath, plot = ggscatterplot, width = width_, height = height_, units = "in", dpi = 300, limitsize = FALSE)
  }
  
  return(ggscatterplot)
}

run_age_analysis = function(cellSubset) {
  message("Starting analysis of ", str_replace_all(cellSubset, "_", "/"))
  
  # Defaults -------------------------------------------------------------------
  fdrCutoff = 0.05
  pCutoff = 0.05 # for modifier analysis
  medianCellCutoff = 10
  analysis = "ageEffects"
  
  # Set up directories ---------------------------------------------------------
  tablesDir = file.path("results", "CyTOF", analysis)
  statsDir = file.path(tablesDir, "stats")
  level1Dir = file.path(tablesDir, "level1")
  modifierDir = file.path(tablesDir, "modifierAnalysis")
  agingDir = file.path(tablesDir, "immuneAge")
  
  plotDir = file.path("plots", "CyTOF", analysis)
  level1Plots = file.path(plotDir, "level1")
  modifierPlots = file.path(plotDir, "modifierAnalysis")
  agingPlots = file.path(plotDir, "immuneAge")
  
  
  dir.create(statsDir, showWarnings = FALSE, recursive = TRUE)
  dir.create(level1Dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(modifierDir, showWarnings = FALSE, recursive = TRUE)
  dir.create(agingDir, showWarnings = FALSE, recursive = TRUE)
  
  dir.create(level1Plots, showWarnings = FALSE, recursive = TRUE)
  dir.create(modifierPlots, showWarnings = FALSE, recursive = TRUE)
  dir.create(agingPlots, showWarnings = FALSE, recursive = TRUE)
  
  # Load and prep data ---------------------------------------------------------
  filePath = paste0("data/counts/", cellSubset, ".Rdata")
  
  
  load(filePath)
  
  pctParent = matrices[["pctParent"]]
  popCounts = matrices[["popCounts"]]
  parentCounts = matrices[["parentCounts"]]
  
  nodes = matrices[["nodes"]]
  design = matrices[["design"]]
  
  design$studyGroup = "DS"
  design$studyGroup[design$patientType == "Control"] = "TC"
  design$studyGroup = factor(design$studyGroup, levels = c("TC", "DS"))
  design$sex = factor(design$sex, levels = sort(unique(design$sex)))
  
  design = dplyr::filter(design, !control)
  
  # matrices already had IC removed, but now I am just going to make sure they are in the same order
  pctParent = pctParent[, design$sampleName]
  popCounts = popCounts[, design$sampleName]
  parentCounts = parentCounts[, design$sampleName]
  
  weights = parentCounts
  
  table(design$timePoint, design$patientType)
  
  # remove all the T2 and T3
  design = dplyr::filter(design, timePoint == "T1")
  
  table(design$ageAtDraw)
  
  pctParent = pctParent[, design$sampleName]
  popCounts = popCounts[, design$sampleName]
  parentCounts = parentCounts[, design$sampleName]
  
  weights = parentCounts
  
  
  # Run linear model -----------------------------------------------------------
  models = runModel(pctParent, popCounts, weights, design, family = "gaussian", cutOff = medianCellCutoff, includeSex = F)
  
  if (length(models) == 0) {
    message("No significant differences between DS and TC")
    return(NULL)
  }
  
  subsets = nodes
  
  # Create contrasts -----------------------------------------------------------
  stats = ageEffects(models, subsets, design, pctParent, popCounts, weights)
  
  stats = dplyr::filter(stats, depth > 0)
  stats = dplyr::filter(stats, fdr != 0)
  
  level1Results = getLevel1(stats)
  
  # Save results to spreadsheets -----------------------------------------------
  saveToSpreadsheet(stats, file.path(statsDir, paste0(cellSubset, "_stats.xlsx")))
  saveToSpreadsheet(level1Results, file.path(level1Dir, paste0(cellSubset, "_level1.xlsx")))
  
  # Run modifier analysis ------------------------------------------------------
  if (nrow(dplyr::filter(stats, fdr < fdrCutoff)) == 0) {
    message("No significant results for ", cellSubset)
    return(NULL)
  }
  
  modifierAnalysis = modifierAnalysisAgeEffect(stats, design, pctParent, weights, fdrCutoff = fdrCutoff, pCutoff = pCutoff)
  
  modifierHits = modifierAnalysis$modifierHits
  modifierResults = modifierAnalysis$modifierResults
  modifierComparisons = modifierAnalysis$modifierComparisons
  
  # remove lines where mean and sd are 0
  modifierHits = dplyr::filter(modifierHits, !(meanPctParent == 100 & sdPctParent == 0))
  modifierHits = dplyr::filter(modifierHits, !(meanPctParent == 0 & sdPctParent == 0))
  # remove hits were slope is too small
  modifierHits = dplyr::filter(modifierHits, abs(ageAtDraw.trend) > 0.01)
  
  # remove lines where mean and sd are 0
  modifierResults = dplyr::filter(modifierResults, !(meanPctParent == 100 & sdPctParent == 0))
  modifierResults = dplyr::filter(modifierResults, !(meanPctParent == 0 & sdPctParent == 0))
  # remove hits were slope is too small
  modifierResults = dplyr::filter(modifierResults, abs(ageAtDraw.trend) > 0.01)
  
  
  # Save modifier results to spreadsheet ---------------------------------------
  saveModifierResultsToSpreadsheet(modifierHits, modifierResults, modifierComparisons, file.path(modifierDir, paste0(cellSubset, "_modifierAnalysis.xlsx")))
  
  results = list(models = models, stats = stats, modifierHits = modifierHits, modifierResults = modifierResults, modifierComparisons = modifierComparisons, design = design, pctParent = pctParent, popCounts = popCounts, parentCounts = weights)
  dir.create(file.path("analysis_results", analysis), showWarnings = FALSE, recursive = TRUE)
  save(results, file = file.path("analysis_results", analysis, paste0(cellSubset, ".Rdata")))
  
  # Level 1 scatter plots ------------------------------------------------------
  plotTitle = paste0("Age effects in level 1 of ", str_replace_all(cellSubset, "_", "/"))
  filePath = file.path(level1Plots, paste0("scatter_", cellSubset, ".png"))
  scatterPlot = plotAgeTrends(level1Results, design, pctParent, parentCounts, models, colors, plotTitle = plotTitle, filePath = filePath, add_trends = T, significant_only = F)
  
  if (nrow(dplyr::filter(modifierHits, fdr < fdrCutoff)) == 0) {
    message("No significant results for ", cellSubset)
    return(NULL)
  }
  
  # Modifier analysis scatter plots --------------------------------------------
  plotTitle = paste0("Significant age effects after modifier analysis of ", str_replace_all(cellSubset, "_", "/"))
  filePath = file.path(modifierPlots, paste0("scatter_", cellSubset, ".png"))
  scatterPlot = plotAgeTrends(modifierHits, design, pctParent, parentCounts, models, colors, plotTitle = plotTitle, filePath = filePath, add_trends = T, significant_only = T)
  
  
  # Explore immune aging -------------------------------------------------------
  immuneAging = modelImmuneAge(modifierHits, design, pctParent, parentCounts, colors, plotDir, cellSubset = cellSubset, targetFolder = agingPlots, refGroup = "TC", fdrCutoff = 0.05, pCutoff = 0.05)
  
  immuneAgeStats = immuneAging$immuneAgeStats
  
  immungeAgeStats = dplyr::select(immuneAgeStats, -pValue_DS_predicted_vs_actual_age)
  
  saveToSpreadsheet(immuneAgeStats, file.path(agingDir, paste0(cellSubset, "_immuneAging.xlsx")))
  # Return ---------------------------------------------------------------------
  return(NULL)
  #   
}

modelImmuneAge = function(ageTrends, design, pctParent, parentCounts, colors, plotDir, cellSubset = "Subset", refGroup = "TC", fdrCutoff = 0.05, pCutoff = 0.05) {
  # Find the subsets that have shown a trend between age and pctParent
  ageTrends = dplyr::filter(ageTrends, studyGroup == refGroup, fdr <= fdrCutoff)
  
  # Model age as function of pctParent -----------------------------------------
  message("Modeling immune age")
  
  ageModels = list()
  immuneAgeStats = c()
  modelStats = c()
  predictions = c()
  plotData = c()
  
  j = 0
  pb = txtProgressBar(min = j, max = length(unique(ageTrends$subset)), style = 3)
  for (subset_ in unique(ageTrends$subset)) {
    j = j + 1
    setTxtProgressBar(pb, j)
    
    ## Prep data and run model -------------------------------------------------
    # Get and prepare data for current subset
    data_ = data.frame(sampleName = colnames(pctParent),
                       pctParent = as.numeric(pctParent[subset_,]),
                       weights = as.numeric(parentCounts[subset_,]),
                       subset = subset_
    )
    data_ = dplyr::filter(data_, weights > 0)
    if (length(unique(data_$pctParent)) == 1) {next}
    
    # Combine with design matrix
    data_ = dplyr::filter(data_, sampleName %in% design$sampleName)
    data_ = merge(data_, design, by.x = "sampleName", by.y = "sampleName")
    
    # linear modeling
    fit = glm(ageAtDraw ~ pctParent, data = data_[data_$studyGroup == refGroup, ], weights = data_$weights[data_$studyGroup == refGroup])
    
    ageModels[[subset_]] = fit
    
    # Generate statistics and evaluate model ----------------------------------
    stats_ <- emtrends(fit, var = "pctParent")
    stats_ = test(stats_)
    stats_$subset = subset_
    stats_$studyGroup = refGroup
    
    modelStats = rbind(modelStats, stats_)
    
    # predict age for all samples
    data_$predictedAge = predict(fit, newdata = data_)
    data_$delta = data_$predictedAge - data_$ageAtDraw
    
    plotData = rbind(plotData, data_)
    
    # generate data for trend and error bars
    predictions_ <- emmeans(fit, ~ pctParent,
                            at = list(pctParent = seq(min(data_$pctParent), max(data_$pctParent), length.out = 50))
    )
    predictions_ = as.data.frame(predictions_)
    predictions_ <- predictions_ %>%
      rename(ageAtDraw = emmean) %>%
      mutate(lower.CL = ageAtDraw - 1.96 * SE,
             upper.CL = ageAtDraw + 1.96 * SE)
    predictions_$studyGroup = "TC"
    predictions_$subset = subset_
    predictions = rbind(predictions, predictions_)
    
    
    # Calculate stats for difference between DS and TC ------------------------
    residuals_summary <- summary(data_$delta[data_$studyGroup == "DS"])
    
    rmse <- sqrt(mean(data_$delta[data_$studyGroup == "DS"]^2))
    
    # Perform t-test between predicted and actual age
    test_result_DS_predicted_vs_actual_age = t.test(data_$predictedAge[data_$studyGroup == "DS"],
                                                    data_$ageAtDraw[data_$studyGroup == "DS"], paired = TRUE)
    
    # Perform test beween the deltas (delta = predicted - actual age) of DS vs TC
    # We are asking if the difference between predicted and actual age is different between DS and TC
    # which is to say is the model predicting age different in DS vs TC
    test_result_DS_vs_TC = wilcox.test(data_$delta[data_$studyGroup == "DS"],
                                       data_$delta[data_$studyGroup == "TC"], paired = FALSE, exact = FALSE)
    
    # Store stats
    immuneAgeStats = rbind(immuneAgeStats, 
                           data.frame(subset = subset_,
                                      studyGroup = "DS",
                                      refGroup = refGroup,
                                      pValue_DS_predicted_vs_actual_age = test_result_DS_predicted_vs_actual_age$p.value, 
                                      pValue_DS_vs_TC = test_result_DS_vs_TC$p.value,
                                      formula = "predicted - actual",
                                      meanAgeDifference = mean(data_$delta[data_$studyGroup == "DS"]),
                                      sdAgeDifference = sd(data_$delta[data_$studyGroup == "DS"]),
                                      rmse = rmse
                           )
    )
    
  }
  
  close(pb)
  
  # Compute fdrs
  immuneAgeStats$fdr = p.adjust(immuneAgeStats$pValue_DS_vs_TC, method = "fdr")
  immuneAgeStats = dplyr::relocate(immuneAgeStats, fdr, .after = pValue_DS_vs_TC)
  
  immuneAgeStats$immuneAgeInDS = ifelse(immuneAgeStats$meanAgeDifference > 0 & immuneAgeStats$fdr <= fdrCutoff, "up",
                                        ifelse(immuneAgeStats$meanAgeDifference < 0 & immuneAgeStats$fdr <= fdrCutoff, "down",
                                               "ns"))
  
  # Generate plots -------------------------------------------------------------
  message("Generating plots")
  
  nCols = ceiling(sqrt(length(unique(plotData$subset))))
  
  height_ = 2*nCols + 0.5
  
  width_ = 3*nCols
  
  targetFolder = file.path(plotDir, "immuneAge")
  plotTitle = paste0(cellSubset, " immune age trends and scatter")
  
  ## Scatter and trends plot --------------------------------------------------- 
  message("  Scatter and trends plot")
  ggscatterplot = ggplot(plotData) +
    # add the prediction data
    geom_ribbon(data = predictions, aes(x = pctParent, ymin = lower.CL, ymax = upper.CL, fill = studyGroup), alpha = 0.2) +
    geom_line(data = predictions, aes(x = pctParent, y = ageAtDraw, color = studyGroup), linewidth = 1, linetype = 4) +
    geom_point(aes(fill = studyGroup, x = pctParent, y = ageAtDraw), color = "black", shape = 21, alpha = 0.5) + 
    # additional settings
    facet_wrap(~ subset, scales = "free", ncol = nCols) +
    theme_bw() +
    scale_fill_manual(values = colors$studyGroupShort) + 
    # scale_fill_manual(values = c(colors$studyGroupShort, colors$studyGroupFit))  +
    scale_color_manual(values = colors$studyGroupShort) +
    labs(x = "percent parent", y = "age", fill = "study group",
         title = plotTitle ) +
    guides(color = "none") +
    theme(strip.text.x = element_text(size = 7))
  
  if (is.null(targetFolder)) {
    print(ggscatterplot)
  } else {
    filePath = file.path(targetFolder, paste0(cellSubset, "_trends_and_scatter.png"))
    ggsave(filePath, plot = ggscatterplot, width = width_, height = height_, units = "in", dpi = 300, limitsize = FALSE)
  }
  
  
  ## Boxplot of difference between predicted and actual ages ------------------- 
  plotTitle = paste0(cellSubset, " immune age distribution DS vs TC")
  
  sigSubsets = dplyr::filter(immuneAgeStats, pValue_DS_vs_TC <= pCutoff)$subset
  
  myTest = function(x, y) {
    wilcox.test(x, y, paired = FALSE, exact = FALSE)
  }
  
  printSignificance= function(x) {
    if (x < 0.01) {
      return("**")
    } else if (x < 0.05) {
      return("*")
    } else {
      return("ns")
    }
    # return(signif(x, digits = 3))
  }
  
  message("  Box plots DS vs TC")
  ggboxplot = ggplot(dplyr::filter(plotData, weights > 0), aes(x = studyGroup, y = delta, fill = studyGroup)) +
    # boxes and jitter
    geom_boxplot(color = "black", alpha = 0.5, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 2, color = "black", stroke = 0.5, position = position_jitterdodge(dodge.width = 0.75))  +
    # significance
    geom_signif(data = dplyr::filter(plotData, subset %in% sigSubsets), 
                comparisons = list(c("DS", "TC")), 
                map_signif_level = printSignificance, 
                test = myTest) +
    # additional settings
    facet_wrap(~ subset, scales = "free", ncol = nCols) +
    theme_bw() +
    scale_fill_manual(values = colors$studyGroupShort) + 
    labs(x = "age", y = "percent parent", fill = "study group",
         title = plotTitle ) +
    theme(strip.text.x = element_text(size = 7)) +
    # expand y axis so that the significance markers are not cut off
    scale_y_continuous(expand = expansion(mult = c(0.1, .2)))
  
  if (is.null(targetFolder)) {
    print(ggboxplot)
  } else {
    filePath = file.path(targetFolder, paste0(cellSubset, "_boxplots.png"))
    ggsave(filePath, plot = ggboxplot, width = width_, height = height_, units = "in", dpi = 300, limitsize = FALSE)
  }
  
  
  ## Volcano plot --------------------------------------------------------------
  message("  Volcano plot")
  # Figure out how many subsets we are going to label in the plot
  nSigs = sum(immuneAgeStats$fdr <= fdrCutoff)
  nDisplay = min(5, nSigs)
  
  plotTitle = paste0(str_replace_all(cellSubset, "_", "/"), " immune system aging in DS")
  xlims = max(abs(immuneAgeStats$meanAgeDifference))
  
  ggvolcano = ggplot(immuneAgeStats, aes(x = meanAgeDifference, y = -log10(fdr))) +
    geom_point(aes(fill = immuneAgeInDS), size = 3, shape = 21, alpha = 0.8) +
    # guide lines
    geom_hline(yintercept = -log10(fdrCutoff), linetype = 2, color = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50")
  
  # add labels to the top 5 (or less) most significant subsets
  if (nSigs > 0) {
    ggvolcano = ggvolcano + 
      geom_text_repel(data = top_n(immuneAgeStats, -5, fdr), aes(label = subset), 
                      box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50', size = 2) 
  } 
  
  # Additional settings
  ggvolcano = ggvolcano + 
    theme_bw() +
    scale_fill_manual(values = colors$volcano) +
    labs(x = "Mean immune system age difference of DS against TC", y = "-log10(fdr)", fill = "Immune system\nage in DS",
         title = plotTitle) +
    xlim(-xlims, xlims)
  
  if (is.null(targetFolder)) {
    print(ggvolcano)
  } else {
    filePath = file.path(targetFolder, paste0(cellSubset, "_volcano.png"))
    ggsave(filePath, plot = ggvolcano, width = 8, height = 5, units = "in", dpi = 300, limitsize = FALSE)
  }
  
  figures = list(ggscatterplot, ggboxplot, ggvolcano)
  ageing = list(ageModels = ageModels, modelStats = modelStats)
  
  # Return the results ---------------------------------------------------------
  return(list(immuneAgeStats = immuneAgeStats, ageModeling = ageing, figures = figures))
  
}