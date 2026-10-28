source("R/common.R")

#' runAbTiterModel = function(pctParent, popCounts, weights, design, family = "gaussian", includeSex = T, cutOff = 10) {
#'   #'@title runModel
#'   #'@description This function runs a linear model on the CyTOF data. It uses the pctParent, popCounts, weights and design data to fit a linear model using age, sex, and studyGroup as covariates. The function returns a list of models. Fits are only returned for subsets where the median count of at least on study group is greater than the cutOff.
#'   #'@param pctParent Matrix of percent parent data, with rows corresponding to subsets and columns corresponding to samples.
#'   #'@param popCounts Matrix of population counts, with rows corresponding to subsets and columns corresponding to samples.
#'   #'@param weights Matrix of weights (parent counts), with rows corresponding to subsets and columns corresponding to samples.
#'   #'@param design Table with meta data
#'   #'@param family Modeling family
#'   #'@param cutOff Optional (default = 10), minimum median count for a subset to be included in the analysis
#'   #'@return A list of fits.
#'   
#'   models = list()
#'   
#'   message("Running linear model using ", family, " response for AbTiter at T3 for ", nrow(pctParent), " subsets")
#'   
#'   pb <- txtProgressBar(min = 0, max = nrow(pctParent), style = 3)
#'   i = 0
#'   
#'   zeroCounts = c()
#'   
#'   for (subset_ in row.names(pctParent)) {
#'     # update the progress bar
#'     i = i + 1
#'     setTxtProgressBar(pb, i)
#'     
#'     data = data.frame(ALId = colnames(pctParent),
#'                       pctParent = as.numeric(pctParent[subset_,]),
#'                       count = as.numeric(popCounts[subset_,]),
#'                       weights = as.numeric(weights[subset_,])
#'     )
#'     
#'     data = merge(data, design, by.x = "ALId", by.y = "ALId")
#'     # summarize data, get median for count by studygroup
#'     medianCount = data %>%
#'       group_by(studyGroup) %>%
#'       summarize(medianCount = median(count))
#'     
#'     if (!any(medianCount$medianCount >= cutOff)) {
#'       zeroCounts = c(zeroCounts, subset_)
#'       next
#'     }
#'     
#'     if (length( unique(data$pctParent[!is.na(data$pctParent)]) ) <= 1) {
#'       zeroCounts = c(zeroCounts, subset_)
#'       next
#'     }
#'     
#'     if (family == "gaussian") {
#'       data$weights[data$weights == 0] = 1e-10
#'       
#'       # if (includeSex) {
#'       #   fit = glm(pctParent ~ studyGroup*sex*ageAtDraw, data = data, weights = data$weights)
#'       # } else {
#'       fit = glm(abTiter_T3 ~ studyGroup * pctParent, data = data, weights = data$weights)
#'       # }
#'       
#'     } else if (family == "binomial") {
#'       
#'       # # have to think about the weights - with a binomial model, they have to be integers (correspond to the "number of trials")
#'       # data$weights = round(data$weights)
#'       # data$weights = as.integer(data$weights)
#'       # 
#'       # data$pctParent = data$pctParent/100
#'       # fit = glm(pctParent ~ studyGroup*sex*ageAtDraw, data = data, weights = data$weights, family = binomial(link = "logit"))
#'     } else if (family == "beta") {
#'       # need to implement beta
#'       # this might help with including that + random effect even maybe?
#'       # https://www.reddit.com/r/rstats/comments/16azcem/family_choice_for_glm_with_01_dependent_variable/
#'     }
#'     
#'     models[[subset_]] = fit
#'     
#'   }
#'   # close the progress bar
#'   close(pb)
#'   
#'   return(models)
#' }
#' 
#' testVaxResponse = function(models, subsets, design, pctParent, popCounts, weights) {
#'   
#'   # Message
#'   message("Testing for effects of vax response in ", length(models), " subsets")
#'   
#'   stats = c()
#'   pb <- txtProgressBar(min = 0, max = length(models), style = 3)
#'   i = 0
#'   emm_options("msg.interaction" = FALSE)
#'   AICs = c()
#'   for (subset_ in names(models)) {
#'     i = i + 1
#'     setTxtProgressBar(pb, i)
#'     
#'     # try to run model, if it fails, add an empty line for stats
#'     tryCatch({
#'       trends <- emtrends(models[[subset_]], 
#'                          ~ studyGroup, 
#'                          var = "pctParent")
#'       
#'       trends = test(trends)
#'       
#'       trends$subset = subset_
#'       
#'     }, error = function(e) {
#'       print(subset_)
#'       trends = data.frame(subset = subset_)
#'     })
#'     
#'     stats = rbind(stats, trends)
#'     
#'     AICs = c(AICs, AIC(models[[subset_]]))
#'     
#'     
#'   }
#'   # close the progress bar
#'   close(pb)
#'   
#'   emm_options("msg.interaction" = TRUE)
#'   
#'   names(AICs) = names(models)
#'   stats = merge(stats, subsets, by = "subset")
#'   
#'   
#'   stats$comparison = paste0(stats$studyGroup, " ab Titer effect")
#'   comparisons = unique(stats$comparison)
#'   
#'   
#'   message("Calculating summary statistics")
#'   summaryStats = c()
#'   
#'   pb <- txtProgressBar(min = 0, max = length(comparisons), style = 3)
#'   i = 0
#'   
#'   summaryStats = c()
#'   for (comparison_ in comparisons) {
#'     stats_ = stats %>% 
#'       dplyr::filter(comparison == comparison_) %>%
#'       dplyr::select(comparison, subset, studyGroup)
#'     
#'     studyGroup_ = unique(stats_$studyGroup)
#'     
#'     study_samples = dplyr::filter(design, studyGroup == studyGroup_)$ALId
#'     
#'     counts_study = popCounts[stats_$subset , study_samples]
#'     pctParent_study = pctParent[stats_$subset, study_samples]
#'     weights_study = weights[stats_$subset, study_samples]
#'     
#'     # calculate the median and sd of the counts
#'     # we have to make sure that we do this based on only the samples that had a 
#'     # nonzero count in the parent population (which is what we use the weights for)
#'     for (iSubset in 1:nrow(stats_) ) {
#'       study_idx = weights_study[iSubset, ] > 0
#'       
#'       counts_study_subset = as.numeric(counts_study[iSubset, study_idx])
#'       pctParent_study_subset = as.numeric(pctParent_study[iSubset, study_idx])
#'       
#'       stats_$medianCounts[iSubset] =   median(counts_study_subset)
#'       stats_$meanCounts[iSubset] =     mean(counts_study_subset)
#'       stats_$sdCounts[iSubset] =       sd(counts_study_subset)
#'       
#'       stats_$medianPctParent[iSubset] = median(pctParent_study_subset)
#'       stats_$meanPctParent[iSubset] =   mean(pctParent_study_subset)
#'       stats_$sdPctParent[iSubset] =     sd(pctParent_study_subset)
#'       
#'     }
#'     
#'     #merge into stats
#'     summaryStats = rbind(summaryStats, stats_)
#'     i = i + 1
#'     setTxtProgressBar(pb, i)
#'   }
#'   
#'   close(pb)
#'   
#'   summaryStats = dplyr::select(summaryStats, -studyGroup)
#'   
#'   stats = merge(stats, summaryStats, by = c("comparison", "subset"))
#'   
#'   stats = multipleTestCorrection(stats)
#'   
#'   # add AICs to stats
#'   stats$AIC = AICs[match(stats$subset, names(AICs))]
#'   
#'   
#'   return(stats)
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 

modifierAnalysisAbTiter = function(stats, design, pctParent, weights, fdrCutoff = 0.05, pCutoff = 0.05, medianCountCutoff = 10) {
  
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
                           !is.na(fdr))
    
    stats_ = dplyr::filter(stats, fdr < fdrCutoff)
    
    if (nrow(stats_) == 0) {
      message(paste0("No significant differences in ", comparison_, ". Skipping."))
      next
    }
    
    message(paste0("Starting modifier analysis of Ab Titer effects in ", studyGroup_, "."))
    
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
        
        
        data_x = data.frame(ALId = colnames(pctParent),
                            pctParent = as.numeric(pctParent[population_$subset,]),
                            weights = as.numeric(weights[population_$subset,])
        )
        data_x = dplyr::filter(data_x, ALId %in% design_$ALId)
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
          
          data_y = data.frame(ALId = colnames(pctParent),
                              pctParent = as.numeric(pctParent[subset_,]),
                              weights = as.numeric(weights[subset_,])
          )
          data_y = dplyr::filter(data_y, ALId %in% design_$ALId)
          data_y$weights[data_y$weights == 0] = 1e-10
          data_y$subset = "child"
          
          data_ = rbind(data_x, data_y)
          
          data_ = merge(data_, design_, by.x = "ALId", by.y = "ALId")
          
          fit = glm(abTiter_T3 ~ pctParent * subset, data = data_,  weights = weights)
          
          trends = emtrends(fit, ~ subset, var = "pctParent")
          
          contrast_results = contrast(trends, method = "pairwise")
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
                        pctParent.trend.child = pctParent.trend, fdrChild = fdr)%>% 
          dplyr::mutate(
            pctParent.trend.parent = population_$pctParent.trend,
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
    
    
    report = distinct(report)
    
    fullReport = rbind(fullReport, report)
    
    fullModifierComparisons = rbind(fullModifierComparisons, modifierComparisons)
    
    if (!is.null(removed_from_report)) {
      removed_from_report = distinct(removed_from_report)
    }
  }
  {
  # Cleanup reports and return -------------------------------------------------
  }
  condensedReport = fullReport %>% 
    dplyr::select(comparison = comparison, subset, studyGroup,
                  terminal_node, modifier_1, modifier_2, 
                  pctParent.trend, p.value, fdr, path, 
                  medianCounts,	meanCounts,	sdCounts,
                  medianPctParent,	meanPctParent,	sdPctParent
    ) %>% 
    dplyr::distinct() %>%
    dplyr::arrange(fdr)
  
  # only keep the single hit per terminal_node with the largest absolute effect
  modifierHits = condensedReport %>%
    dplyr::group_by(terminal_node) %>%
    dplyr::slice_max(order_by = abs(pctParent.trend), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  modifierResults = fullReport
  modifierComparisons = fullModifierComparisons
  
  # return(list(modifierResults = modifierResults, modifierComparisons = modifierComparisons, removedFromReport = removed_from_report))
  return(list(modifierHits = modifierHits, modifierResults = modifierResults, modifierComparisons = modifierComparisons))
}



makeAbBoxplots = function(stats, pctParent, weights, design, models, colors, filePath = NULL, title = "Box plot", significant_only = T, add_stats = F) {
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
  data = tidyr::gather(data, key = "ALId", value = "pctParent", -subset)
  
  # merge the data with the design matrix
  data = merge(data, design, by = "ALId")
  
  
  # get the weights for each row from the weights matrix (rows are subset and columns are ALId)
  data$weights = NA
  for (subset_ in unique(data$subset)) {
    data$weights[data$subset == subset_] = unlist(weights[subset_, match(data$ALId[data$subset == subset_], colnames(weights))])
  }
  
  # hist(data$weights)
  # keep only data with weights > 0
  data = dplyr::filter(data, weights > 0)
  
  predictions = c()
  for (iRow in 1:nrow(stats)) {
    if (stats$fdr[iRow] > 0.05) {
      next
    }
    subset_ = stats$subset[iRow]
    studyGroup_ = stats$studyGroup[iRow]
    # generate data for trend and error bars
    pctParent_range = dplyr::filter(data, subset == subset_, studyGroup == studyGroup_)$pctParent
    pctParent_range = seq(min(pctParent_range), max(pctParent_range), length.out = 50)
    
    predictions_ <- emmeans(models[[subset_]], ~ pctParent * studyGroup,
                            at = list(pctParent = pctParent_range,
                                      studyGroup = studyGroup_)
    )
    
    predictions_ = as.data.frame(predictions_)
    predictions_ <- predictions_ %>%
      rename(abTiter_T3 = emmean) %>%
      mutate(lower.CL = abTiter_T3 - 1.96 * SE,
             upper.CL = abTiter_T3 + 1.96 * SE) %>%
      mutate(abTiter_T3 = abTiter_T3,
             lower.CL = lower.CL,
             upper.CL = upper.CL)
    predictions_$subset = subset_
    predictions = rbind(predictions, predictions_)
  }
  
  nCols = ceiling(sqrt(length(unique(data$subset))))
  
  height_ = 2*nCols + 0.5
  
  width_ = 3*nCols
  
  # targetFolder = file.path(plotDir, "immuneAge")
  plotTitle = paste0(cellSubset, " ab titer response trends and scatter")
  
  ggscatterplot = ggplot(data) +
    # add the prediction data
    geom_ribbon(data = predictions, aes(x = pctParent, ymin = lower.CL, ymax = upper.CL, fill = studyGroup), alpha = 0.2) +
    geom_line(data = predictions, aes(x = pctParent, y = abTiter_T3, color = studyGroup), linewidth = 1, linetype = 4) +
    geom_point(aes(fill = studyGroup, x = pctParent, y = abTiter_T3), color = "black", shape = 21, alpha = 0.5) + 
    # additional settings
    facet_wrap(~ subset, scales = "free", ncol = nCols) +
    theme_bw() +
    scale_fill_manual(values = colors$studyGroup) + 
    scale_color_manual(values = colors$studyGroup) +
    labs(x = "percent parent", y = "AbTiter at T3", fill = "study group",
         title = plotTitle ) +
    guides(color = "none") +
    ylim(0, 151)
  theme(strip.text.x = element_text(size = 7))
  
  ggscatterplot
  # if (is.null(filePath)) {
  #   print(ggboxplot)
  # } else {
  #   ggsave(filePath, plot = ggboxplot, width = 14, height = 11, units = "in", dpi = 300)
  # }
  
  return(ggscatterplot)
}