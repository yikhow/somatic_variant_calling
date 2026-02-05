# ==============================================================================
# Script: somatic_visualization.R
# Description: Visualization of somatic variants including Ti/Tv distribution.
# ==============================================================================

library(maftools)

# Get args from environment or set default
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
    PROJECT_DIR <- args[1]
} else {
    PROJECT_DIR <- getwd() 
}

message("[INFO] Project Directory: ", PROJECT_DIR)

# Paths
FILE_MAF   <- file.path(PROJECT_DIR, "results", "somatic_functional.maf")
FILE_CGC   <- file.path(PROJECT_DIR, "resources", "CancerGeneCensus", "Census_allSun_Feb_1_09_04_28_2026.csv")
OUTPUT_DIR <- file.path(PROJECT_DIR, "plots")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# Verify Files
if(!file.exists(FILE_MAF)) stop("MAF file not found: ", FILE_MAF)

# ----------------------- [Analysis] ------------------------------------------ #
message("[INFO] Loading MAF...")
laml_raw <- read.maf(maf = FILE_MAF, verbose = FALSE)

# Optional: Filter by CGC if file exists
if(file.exists(FILE_CGC)){
    cgc_ref  <- read.csv(FILE_CGC, header = TRUE, check.names = FALSE)
    official_genes <- as.character(cgc_ref$`Gene Symbol`)
    laml <- subsetMaf(maf = laml_raw, genes = official_genes, mafObj = TRUE)
    message("[INFO] Filtered by COSMIC CGC.")
} else {
    laml <- laml_raw
    message("[INFO] COSMIC CGC not found, using full MAF.")
}

# ----------------------- [Step 1] MAF Summary Plot --------------------------- #
message("[INFO] Step 1: Generating MAF summary...")
png(filename = file.path(OUTPUT_DIR, "maf_summary.png"), width = 1000, height = 1200, res = 150)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

# ----------------------- [Step 2] Ti/Tv Distribution ------------------------- #
message("[INFO] Step 2: Generating Ti/Tv distribution plot...")
# Calculate Ti/Tv summary
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)

png(filename = file.path(OUTPUT_DIR, "titv_distribution.png"), width = 1000, height = 800, res = 150)
plotTiTv(res = laml.titv)
dev.off()

# ----------------------- [Step 3] Lollipop Plot (KMT2C) ---------------------- #
message("[INFO] Step 3: Generating lollipop plot for KMT2C...")
png(filename = file.path(OUTPUT_DIR, "lollipop_KMT2C.png"), width = 1200, height = 800, res = 150)
lollipopPlot(
    maf = laml, 
    gene = "KMT2C", 
    AACol = "Protein_Change", 
    showDomainLabel = FALSE, 
    labelPos = "all", 
    repel = TRUE
)
dev.off()

message("[SUCCESS] All plots (Summary, Ti/Tv, Lollipop) generated in: ", OUTPUT_DIR)
