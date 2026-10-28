getGatingSet = function(design) {
  gs = list()
  print("Setting up gating sets for each batch...")
  for (batch_ in unique(design$batch)) {
    samples_ = dplyr::filter(design, batch == batch_)
    print(paste0("batch: ", batch_, " samples: ", nrow(samples_)))
    flowSet = read.flowSet(samples_$transformedFcs, truncate_max_range = FALSE)
    gs[[batch_]] = GatingSet(flowSet)
  }
  return(gs)
}

# add_new_node_workspace = function(gs, node, gate, parent) {
#   nodeToAdd <- if (parent == "root") paste0("/", node) else paste0(parent, "/", node)
#   existingNodes = gs_get_pop_paths(gs, path = "full")
#   if (nodeToAdd %in% existingNodes) return(nodeToAdd)
#   gs_pop_add(gs, gate, parent = parent)
#   return(nodeToAdd)
# }

add_new_node_workspace <- function(gs, node, gate, parent) {
  nodeToAdd <- if (identical(parent, "root")) paste0("/", node) else paste0(parent, "/", node)
  existingNodes <- flowWorkspace::gs_get_pop_paths(gs, path = "full")
  if (nodeToAdd %in% existingNodes) return(nodeToAdd)
  
  # ---- DEBUG + normalize gate list names to match the GatingSet ----
  samps <- flowWorkspace::sampleNames(gs)
  cat("[add_node] parent=", parent, " node=", node, " n_samps=", length(samps), "\n", sep="")
  cat("[add_node] first names(gate): ", paste(head(names(gate)), collapse=", "), "\n", sep="")
  
  if (!is.list(gate)) gate <- rep(list(gate), length(samps))
  names(gate) <- samps
  gate <- gate[samps]
  
  # Parent must exist
  if (!(parent %in% existingNodes)) {
    stop(sprintf("Parent '%s' not found in GatingSet. First nodes: %s",
                 parent, paste(utils::head(existingNodes, 10), collapse=" | ")))
  }
  
  flowWorkspace::gs_pop_add(gs, gate, parent = parent)
  nodeToAdd
}


# ---------- NEW: RA/cellType/parent-aware resolver ----------
.pick_node_idx <- function(nm, parent_path, gateList, prev_idx = NA_integer_) {
  cand <- which(gateList$node == nm)
  if (length(cand) == 0) return(NA_integer_)
  if (length(cand) == 1) return(cand)
  
  cellType_ <- ifelse(str_detect(parent_path, "CD4"), "CD4",
                      ifelse(str_detect(parent_path, "CD8"), "CD8", NA))
  RA_ <- ifelse(str_detect(parent_path, "RA\\+-"), "RA+-",
                ifelse(str_detect(parent_path, "RA\\+(?!-)"), "RA+", NA))
  
  if (!is.na(cellType_)) cand <- cand[gateList$cellType[cand] == cellType_]
  if (length(cand) == 1) return(cand)
  
  if (!is.na(RA_))       cand <- cand[gateList$RA[cand] == RA_]
  if (length(cand) == 1) return(cand)
  
  if (!is.na(prev_idx)) {
    dep_hits <- vapply(cand, function(i) {
      di <- gateList$dependsOn[[i]]
      if (is.null(di)) return(FALSE)
      if (is.character(di)) return(FALSE)
      prev_idx %in% unlist(di, use.names = FALSE)
    }, logical(1))
    if (sum(dep_hits) >= 1) return(cand[which(dep_hits)][1])
  }
  
  ph <- which(gateList$parent[cand] == parent_path)
  if (length(ph) >= 1) return(cand[ph[1]])
  
  cand[1]
}
# ----------------------------------------------------------- for supervised signatures, only

# --- helper: make a per-sample gate list named like the gs (GatingSet) ---
.fix_gate_names <- function(gs, gate_obj) {
  samps <- flowWorkspace::sampleNames(gs)
  if (!is.list(gate_obj)) {
    gate_obj <- rep(list(gate_obj), length(samps))
  }
  names(gate_obj) <- samps
  gate_obj <- gate_obj[samps]
  gate_obj
}
# -----------------------------------------------------------for supervised signatures, only
# Build the concrete root path in gs/gates, using resolver to avoid RA mixups
getRootPathTcells = function(pathToGate, gateList, gates, gs = NULL) {
  pathToGate = str_split(pathToGate, "/")[[1]]
  parent = "root"
  prev_idx <- NA_integer_
  
  for (i in seq_along(pathToGate)) {
    nodeToAdd = pathToGate[i]
    if (is.null(gs)) { message(nodeToAdd) }
    
    iGate <- .pick_node_idx(nodeToAdd, parent, gateList, prev_idx)
    if (is.na(iGate)) {
      stop(paste0("Could not resolve gate for node '", nodeToAdd, "' under parent '", parent, "'"))
    }
    
    # add declared dependencies first
    additionalGates = unlist(gateList$dependsOn[iGate])
    if (length(additionalGates) > 1 | !(any(is.na(additionalGates)))) {
      for (k in seq_along(additionalGates)) {
        if (is.character(additionalGates[k])) {
          idx <- .pick_node_idx(additionalGates[k], parent, gateList, prev_idx)
        } else {
          idx <- gateList$idx[additionalGates[k]]
        }
        if (is.null(gs)) {
          print(paste0("Extra node: ", gateList$node[idx]))
        } else {
          #add_new_node_workspace(gs, gateList$node[idx], gates[[idx]], parent)
          # ----------------------------------------------------------- for supervised signatures, only
          gate_obj <- .fix_gate_names(gs, gates[[idx]])
          cat("[dep] parent=", parent, " node=", gateList$node[idx],
              " n_samps=", length(flowWorkspace::sampleNames(gs)),
              " gate_len=", length(gate_obj), "\n", sep="")
          add_new_node_workspace(gs, gateList$node[idx], gate_obj, parent)
          # -----------------------------------------------------------for supervised signatures, only
        }
      }
    }
    
    idx <- gateList$idx[iGate]
    if (is.null(gs)) {
      parent <- if (parent == "root") paste0("/", nodeToAdd) else paste0(parent, "/", nodeToAdd)
      print(parent)
    } else {
      #parent <- add_new_node_workspace(gs, nodeToAdd, gates[[idx]], parent)
      # ----------------------------------------------------------- for supervised signatures, only
      gate_obj <- .fix_gate_names(gs, gates[[idx]])
      cat("[add] parent=", parent, " node=", nodeToAdd,
          " n_samps=", length(flowWorkspace::sampleNames(gs)),
          " gate_len=", length(gate_obj), "\n", sep="")
      parent <- add_new_node_workspace(gs, nodeToAdd, gate_obj, parent)
      # -----------------------------------------------------------for supervised signatures, only
    }
    prev_idx <- idx
  }
  return(parent)
}

getRawCounts = function(gs, design) {
  population_counts = lapply(unique(design$batch), function(batch_) {
    gs_pop_get_count_fast(gs[[batch_]], format = "long")
  })
  counts = bind_rows(population_counts)
  counts$PercentParent = counts$Count / counts$ParentCount * 100
  counts$PercentParent[is.na(counts$PercentParent)] = NA
  return(counts)
}

processGatingPaths = function(summary, max_depth = 1) {
  columns = colnames(summary)
  if (!"path" %in% columns) summary$path = summary$Population
  summary$depth = str_count(summary$path, "/")
  root_path_length = max(summary$depth) - max_depth
  pathing = str_split(str_remove(summary$path, '^/'), "/", simplify = TRUE)
  summary$terminal_node = pathing[cbind(1:nrow(pathing), summary$depth)]
  summary$depth = summary$depth - root_path_length
  if (max_depth == 1) return(summary)
  max_modifiers = max_depth - 1
  for (i_modifier in 1:max_modifiers) {
    column_name = paste0("modifier_", i_modifier)
    summary[[column_name]] = ""
    idx = which(summary$depth > i_modifier)
    summary[idx, column_name] = pathing[cbind(idx, root_path_length + i_modifier)]
  }
  return(summary)
}

assignRoot = function(nodes, rootPaths) {
  escaped_rootPaths <- sapply(rootPaths, function(path) {
    path <- str_replace_all(path, "\\+", "\\\\+")
    path <- str_replace_all(path, "\\-", "\\\\-")
    return(path)
  })
  nodes <- nodes %>%
    rowwise() %>%
    mutate(
      rootPath = {
        match <- sapply(names(escaped_rootPaths), function(x) {
          str_detect(path, paste0(escaped_rootPaths[[x]], "(/|$)"))
        })
        if (any(match)) names(escaped_rootPaths)[which(match)[1]] else NA
      },
      pathInCompartment = ifelse(!is.na(rootPath),
                                 str_remove(path, paste0(escaped_rootPaths[rootPath], "(/|$)")),
                                 NA)
    ) %>% ungroup()
  return(nodes)
}

convertToMatrix = function(counts, depths = c(1, 2, 3)) {
  counts = dplyr::filter(counts, depth %in% depths)
  pctParent    = spread(dplyr::select(counts, sampleName, PercentParent, subset), key = sampleName, value = PercentParent)
  parentCounts = spread(dplyr::select(counts, sampleName, ParentCount,   subset), key = sampleName, value = ParentCount)
  popCounts    = spread(dplyr::select(counts, sampleName, Count,         subset), key = sampleName, value = Count)
  rownames(pctParent)    = pctParent$subset;    pctParent    = dplyr::select(pctParent, -subset)
  rownames(parentCounts) = parentCounts$subset; parentCounts = dplyr::select(parentCounts, -subset)
  rownames(popCounts)    = popCounts$subset;    popCounts    = dplyr::select(popCounts, -subset)
  return(list(pctParent = pctParent, parentCounts = parentCounts, popCounts = popCounts))
}

gateAndSave = function(path_) {
  print(paste0("IMPACD gating all batches for: ", rootPaths$name[path_]) )
  
  gs = getGatingSet(design)
  
  for (batch_ in unique(design$batch)) {
    print(paste0("Batch ", batch_,": root"))
    x_ = getRootPathTcells(rootPaths$path[path_], gateListAll, gatesAll[[batch_]], gs[[batch_]])
    
    print(paste0("Batch ", batch_,": IMPACD"))
    rootPath = str_split(rootPaths$path[path_], "/")[[1]]
    
    # assemble modifiers for this batch
    modifierNodes = dplyr::filter(gateListAll, type == "modifier")
    modifierNodes = dplyr::filter(modifierNodes, !(node %in% rootPath))
    if (any(str_detect(rootPaths$path[path_], "CD27"))) {
      modifierNodes = dplyr::filter(modifierNodes, node != "CD27+")
    }
    
    modifierGates = gatesAll[[batch_]][modifierNodes$idx]
    names(modifierGates) = modifierNodes$node
    
    gatemetadata = data.frame(
      name = modifierNodes$node,
      gate_type = "2D",
      marker1 = NA,
      marker2 = NA,
      original_gate_position = NA,
      logical_gate = FALSE,
      transformed_gate_position = NA
    )
    
    # -------- NEW: prune invalid/missing modifier gates for this batch --------
    ok <- vapply(modifierGates, function(g) !is.null(g) && inherits(g, "filter"), logical(1))
    if (any(!ok)) {
      dropped <- modifierNodes$name[!ok]
      message("Dropping ", sum(!ok), " invalid modifier gate(s) in batch ", batch_, ": ",
              paste(dropped, collapse = ", "))
      modifierNodes <- modifierNodes[ok, , drop = FALSE]
      modifierGates <- modifierGates[ok]
      gatemetadata  <- gatemetadata[ok, , drop = FALSE]
    }
    
    if (length(modifierGates) == 0L) {
      message("No valid modifier gates for batch ", batch_, " — skipping level2Gating for this batch.")
    } else {
      level2_paths = level2Gating(
        gs[[batch_]],
        modifierNodes$node, gatemetadata, modifierGates,
        paste0("/", rootPaths$path[path_]),
        max_depth = 3,
        level1_nodes = modifierNodes$node,
        cellcountCutoff = 0
      )
    }
    # -------------------------------------------------------------------------
  }
  
  # counts and save as before
  counts = getRawCounts(gs, design)
  
  # underscore file name (keep '/' inside matrices for subsets)
  fname = str_replace_all(rootPaths$name[path_], "/", "_")
  dataFileName = paste0(fname, ".Rdata")
  
  # raw counts
  save(counts, file = file.path(rawCountsDir, dataFileName))
  
  # combine with design and build matrices
  counts = merge(counts, design, by.x = "name", by.y = "fcsName")
  colnames(counts)[colnames(counts) == "name"] = "fcsName"
  counts = dplyr::filter(counts, control == FALSE)
  
  nodes = data.frame(path = unique(counts$Population))
  nodes = processGatingPaths(nodes, max_depth = 3)
  nodes = assignRoot(nodes, rootPaths$path[path_])
  nodes$pathInCompartment = str_remove(nodes$pathInCompartment, "^/")
  
  nodes$subset = paste(rootPaths$name[path_], nodes$pathInCompartment, sep = "/")
  nodes$subset = str_remove(nodes$subset, "/$")
  nodes$subset[is.na(nodes$rootPath)] = NA
  
  counts = merge(counts, nodes, by.x = "Population", by.y = "path")
  
  matrices = convertToMatrix(dplyr::filter(counts, !is.na(subset)), depths = c(0, 1, 2, 3))
  matrices[["nodes"]]  = nodes
  matrices[["design"]] = design
  
  # processed counts
  save(matrices, file = file.path(processedCountsDir, dataFileName))
  
  return("Done")
}


gateSupervisedSignatures = function(gs, design, signatures, gateListAll, gatesAll) {
  uniqueSignatures = unique(signatures$path)
  print(paste0("Gating ", length(uniqueSignatures), " unique signatures across all batches."))
  for (signature_ in uniqueSignatures) {
    print(signature_)
    for (batch_ in unique(design$batch)) {
      x_ = getRootPathTcells(signature_, gateListAll, gatesAll[[batch_]], gs[[batch_]])
    }
  }
  for (batch_ in unique(design$batch)) recompute(gs[[batch_]])
  counts = getRawCounts(gs, design)
  return(counts)
}
