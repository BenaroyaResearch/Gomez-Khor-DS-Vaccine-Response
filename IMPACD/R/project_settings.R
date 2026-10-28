# Setup analysis folders
#baseDir = getwd()
baseDir = "/Users/rgomez/Khor_covidvax_response_IMPACD"
dataDir = file.path(baseDir,'data')
resultsDir = file.path(baseDir, "results")
tablesDir = file.path(baseDir, 'tables')
fcsDir = file.path(baseDir, "fcs_files")

rawCountsDir = file.path(dataDir, "raw_counts")
processedCountsDir = file.path(dataDir, "counts")

dir.create(dataDir, recursive = T, showWarnings = F)
dir.create(resultsDir, recursive = T, showWarnings = F)
dir.create(tablesDir, recursive = T, showWarnings = F)
dir.create(fcsDir, recursive = T, showWarnings = F)
dir.create(rawCountsDir, recursive = T, showWarnings = F)
dir.create(processedCountsDir, recursive = T, showWarnings = F)

# Colors ---------------------------------------------------------------
colors = c()

# colors$heatmap = colorRamp2(c(-4, 0, 4), c("#0033a0", "white", "#FE5000")) # colors for the heatmap

col_fun = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")) # colors for the heatmap
colors$heatmap = col_fun

colors$patientType = c("black", "#E6847D")
names(colors$patientType) = c("Control", "Down Syndrome")

# colors$studyGroupShort = c("#40A184", "#E0795C")
colors$studyGroupShort = c("black", "#E6847D")
names(colors$studyGroupShort) = c("TC", "DS")

colors$studyGroup = c("black", "#E6847D")
names(colors$studyGroup) = c("TC", "DS")

# colors$condition = c("#60B1F5", "#F5A01D", "#89F572", "#A0488B")
# names(colors$condition) = unique(design$condition)
# 
# colors$timePoint = c("#CCEDE7", "#8AA19D", "#6D7F7B")
# names(colors$timePoint) = sort(unique(design$timePoint))

colors$sex = c("#F0A8A8", "#60B1F5")
names(colors$sex) = c("female", "male")

colors$compartment = c("#D5EFD3", "#EFD2B9")
names(colors$compartment) = c("global", "antigen_specific")

colors$volcano = c("#2166AC", "#B2182B", "#E7E7E7")
names(colors$volcano) = c("down", "up", "ns")

# make a color map ageAtDraw
# colors$ageAtDraw = colorRamp2(c(min(design$ageAtDraw), max(design$ageAtDraw)), c("#E7E7E7", "#818181")) # colors for the heatmap
colors$ageAtDraw = colorRamp2(c(5, 52), c("#E7E7E7", "#818181")) # colors for the heatmap

# shared paths --------------------------------------------------------
# baseDir = getwd()
# dataDir = file.path(baseDir,'data')
# plotsDir = file.path(baseDir,'plots')
# resultsDir = file.path(baseDir, 'results')
# tablesDir = file.path(baseDir, 'tables')
# fcsDir = file.path(baseDir, "fcs_files")
# 
# dir.create(dataDir, recursive = T, showWarnings = F)
# dir.create(plotsDir, recursive = T, showWarnings = F)
# dir.create(resultsDir, recursive = T, showWarnings = F)
# dir.create(tablesDir, recursive = T, showWarnings = F)
# dir.create(fcsDir, recursive = T, showWarnings = F)