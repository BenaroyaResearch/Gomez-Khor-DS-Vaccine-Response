getLevel1 = function(stats) {
  level1 = c()
  
  for (comparison_ in unique(stats$comparison)) {
    level1_ = dplyr::filter(stats, comparison == comparison_, depth == 1)
    level1_$fdr = p.adjust(level1_$p.value, method = "BH")
    level1 = rbind(level1, level1_)
  }
  
  level1 = dplyr::select(level1, -modifier_1, -modifier_2)
  
  return(level1)
}

multipleTestCorrection = function(stats) {
  # for each unique group run multiple test correction on the p-values
  for (comparison_ in unique(stats$comparison)) {
    stats$fdr[stats$comparison == comparison_] = p.adjust(stats$p.value[stats$comparison == comparison_], method = "BH")
  }
  stats = dplyr::relocate(stats, fdr, .after = p.value)
}

saveModifierResultsToSpreadsheet = function(modifierHits, modifierResults, modifierComparisons, filename) {
  wb <- createWorkbook()
  
  for (comparison_ in unique(modifierHits$comparison)) {
    sheetName_ = paste0(comparison_, " hits")
    sheetName_ = substr(sheetName_, 1, 31)
    addWorksheet(wb, sheetName = sheetName_)
    writeData(wb, sheet = sheetName_, 
              x = dplyr::filter(modifierHits, comparison == comparison_) %>%
                dplyr::arrange(fdr)
    )
  }
  
  for (comparison_ in unique(modifierHits$comparison)) {
    sheetName_ = paste0(comparison_, " + comparisons")
    # truncate at 31 characters
    sheetName_ = substr(sheetName_, 1, 31)
    addWorksheet(wb, sheetName = sheetName_)
    writeData(wb, sheet = sheetName_, 
              x = dplyr::filter(modifierResults, comparison == comparison_)
    )
  }
  
  sheetName_ = "All modifer comparisons"
  sheetName_ = substr(sheetName_, 1, 31)
  addWorksheet(wb, sheetName = sheetName_)
  writeData(wb, sheet = sheetName_, x = modifierComparisons)
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

saveToSpreadsheet = function(stats, filename) {
  wb <- createWorkbook()
  
  if ("comparison" %in% colnames(stats)) {
    mode = "contrasts"
  } else {
    mode = "immune aging"
  }
  
  if (mode == "contrasts") {
    for (comparison_ in unique(stats$comparison)) {
      stats_ = dplyr::filter(stats, comparison == comparison_) %>%
        dplyr::arrange(fdr)
      
      # Write results to Excel file
      addWorksheet(wb, sheetName = comparison_)
      writeData(wb, sheet = comparison_, x = stats_)
    }
  } else {
    stats = stats %>% dplyr::arrange(fdr)
    addWorksheet(wb, sheetName = "Immune system aging")
    writeData(wb, sheet = "Immune system aging", x = stats)
  }
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}