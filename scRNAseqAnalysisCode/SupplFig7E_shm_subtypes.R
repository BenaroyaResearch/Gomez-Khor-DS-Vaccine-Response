## Supplemental Figure 7E --- subtypes of somatic hypermutation
##
## "Frequency of subtypes of somatic hypermutations in spike-specific B cells
## from people with versus without DS." Transitions, transversions, silent and
## nonsilent mutations at the IgH locus in IgG+ cells -- the same sequences as
## the top-left panel of main text Figure 5D, broken down by mutation type.
##
## Input:  inputs/imgt_vquest/IgH_IgG.txt  (IMGT/V-QUEST export; see inputs/README.md)
##         pseudobulk_clean_output/step1_library_qc_metrics.csv (group-label check)
## Output: figures/SupplFig7E_shm_subtypes.pdf/.png
##         figure_source_data/SupplFig7E_shm_subtypes.csv
##         figure_source_data/SupplFig7E_donor_medians.csv
##         figure_source_data/SupplFig7E_statistics.csv
##
## Run:  Rscript figure_panels/SupplFig7E_shm_subtypes.R

local({
  for (p in c("figure_panels/_shared.R", "_shared.R")) if (file.exists(p)) return(source(p, chdir = FALSE))
  stop("_shared.R not found; run from the repo root or from figure_panels/")
})
suppressPackageStartupMessages(library(stringr))

## Locus and productive isotype, matching the submitted panel.
LOCUS   <- "IgH"
ISOTYPE <- "IgG"

## Mutation subtypes, in the published left-to-right order.
METRICS <- c(Transitions   = "Transitions",
             Transversions = "Transversions",
             Silent        = "Silent mutations",
             Nonsilent     = "Nonsilent mutations")

## The legend says "(E) Welch two sample t test", on each donor's MEDIAN count
## per sequence -- the same unit of analysis as Figure 5D.
extractCount <- function(x) tidyr::replace_na(as.numeric(str_extract(as.character(x), "^\\d+")), 0)

path <- file.path(REPO_ROOT, "inputs", "imgt_vquest", sprintf("%s_%s.txt", LOCUS, ISOTYPE))
if (!file.exists(path)) stop("missing V-QUEST export: ", path)

vq <- read_delim(path, delim = "\t", show_col_types = FALSE, name_repair = "minimal")
ti <- c("a>g", "g>a", "c>t", "t>c")
tv <- c("a>c", "c>a", "a>t", "t>a", "g>c", "c>g", "g>t", "t>g")
sumCols <- function(d, subs) rowSums(sapply(subs, function(s) extractCount(d[[paste0("V-REGION ", s)]])))

s7e <- tibble::tibble(
  sequenceId  = vq$`Sequence ID`,
  studyGroup  = factor(dplyr::if_else(str_detect(vq$`Sequence ID`, "TC$"), "TC", "DS"),
                       levels = c("TC", "DS")),
  project     = str_extract(vq$`Sequence ID`, "P635-\\d+"),
  Transitions   = sumCols(vq, ti),
  Transversions = sumCols(vq, tv),
  Silent        = extractCount(vq$`V-REGION Nb of silent mutations`),
  Nonsilent     = extractCount(vq$`V-REGION Nb of nonsilent mutations`))
stopifnot(!any(is.na(s7e$project)))

## The group label is parsed from a sequence-ID suffix, so check it against the
## donor table step 1 wrote rather than trusting the string.
step1 <- read_csv(file.path(OUT_DIR, "step1_library_qc_metrics.csv"), show_col_types = FALSE) %>%
  dplyr::distinct(project, studyGroupStep1 = studyGroup)
chk <- s7e %>% dplyr::distinct(project, studyGroup) %>% dplyr::left_join(step1, by = "project")
if (any(is.na(chk$studyGroupStep1)) || any(as.character(chk$studyGroup) != chk$studyGroupStep1))
  stop("V-QUEST sequence-ID study groups disagree with step 1")
cat(sprintf("%d %s sequences from %s+ cells, %d donors; group labels agree with step 1\n",
            nrow(s7e), LOCUS, ISOTYPE, dplyr::n_distinct(s7e$project)))

long <- s7e %>%
  tidyr::pivot_longer(dplyr::all_of(names(METRICS)), names_to = "metric", values_to = "count") %>%
  dplyr::mutate(metric = factor(metric, levels = names(METRICS)),
                metricLabel = factor(METRICS[as.character(metric)], levels = METRICS))

donorMedians <- long %>%
  dplyr::group_by(metric, metricLabel, studyGroup, project) %>%
  dplyr::summarise(nSequences = dplyr::n(), donorMedian = median(count),
                   donorMean = mean(count), .groups = "drop")

s7eStats <- donorMedians %>%
  dplyr::group_split(metric) %>%
  lapply(function(d) {
    x <- d$donorMedian[d$studyGroup == "TC"]; y <- d$donorMedian[d$studyGroup == "DS"]
    p <- suppressWarnings(tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_))
    tibble::tibble(metric = as.character(d$metric[1]),
                   metricLabel = as.character(d$metricLabel[1]),
                   nDonorCtrl = length(x), nDonorDS = length(y),
                   meanOfDonorMediansCtrl = mean(x), meanOfDonorMediansDS = mean(y),
                   pWelchTOnDonorMedians = p)
  }) %>% dplyr::bind_rows() %>%
  dplyr::mutate(metric = factor(metric, levels = names(METRICS))) %>%
  dplyr::arrange(metric)

print(s7eStats %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))) %>% as.data.frame())

fmtP <- function(p) sprintf("p = %.2f", p)
annot <- long %>%
  dplyr::group_by(metric, metricLabel) %>%
  dplyr::summarise(yTop = max(count), .groups = "drop") %>%
  dplyr::left_join(s7eStats %>% dplyr::select(metric, pWelchTOnDonorMedians), by = "metric") %>%
  dplyr::mutate(label = fmtP(pWelchTOnDonorMedians))

pS7E <- ggplot(long, aes(studyGroup, count, colour = studyGroup)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey30", fill = NA, linewidth = 0.3) +
  geom_jitter(width = 0.2, height = 0.15, size = 0.7, alpha = 0.6) +
  geom_text(data = annot, inherit.aes = FALSE, vjust = 0, size = 3,
            aes(x = 1.5, y = yTop * 1.06, label = label)) +
  facet_wrap(~ metricLabel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = GROUP_LABELS) +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Mutations per sequence",
       title = "Supplemental Figure 7E: subtypes of somatic hypermutation",
       subtitle = sprintf("%s locus, %s+ cells | Welch t-test on per-donor medians", LOCUS, ISOTYPE)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank())

print(pS7E)
savePanel(pS7E, "SupplFig7E_shm_subtypes", width = 9, height = 3.6)

saveSource(
  long %>% dplyr::transmute(studyGroup = as.character(studyGroup),
                            groupLabel = GROUP_LABELS[as.character(studyGroup)],
                            project, sequenceId, locus = LOCUS, cellIsotype = ISOTYPE,
                            metric = as.character(metric),
                            metricLabel = as.character(metricLabel), count),
  "SupplFig7E_shm_subtypes",
  "Mutations of each subtype per sequence (the plotted points)",
  panel = "Suppl 7E", role = "plotted")

saveSource(
  donorMedians %>% dplyr::transmute(studyGroup = as.character(studyGroup),
                                    groupLabel = GROUP_LABELS[as.character(studyGroup)],
                                    project, metric = as.character(metric),
                                    metricLabel = as.character(metricLabel),
                                    nSequences, donorMedian, donorMean) %>%
    dplyr::arrange(metric, studyGroup, as.numeric(sub("^P635-", "", project))),
  "SupplFig7E_donor_medians",
  "Per-donor median and mean mutations per sequence (the values the test is computed from)",
  panel = "Suppl 7E", role = "supporting")

saveSource(s7eStats %>% dplyr::mutate(metric = as.character(metric)),
           "SupplFig7E_statistics",
           "Welch t-test on per-donor medians, per mutation subtype",
           panel = "Suppl 7E", role = "statistics")
