# ============================================================================
# 1. LIBRARIES
# ============================================================================

library(minfi)
library(limma)
library(sva)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(missMethyl)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(corrplot)
library(ggrepel)
library(viridis)
library(scales)
library(RColorBrewer)
library(enrichplot)
library(patchwork)
library(matrixStats)
library(readr)
library(tibble)
library(pROC)
library(networkD3)
library(htmlwidgets)

# ============================================================================
# 2. DIRECTORIES AND PATHS
# ============================================================================

IDAT_DIR         <- ""
SAMPLE_CSV       <- ""
MANIFEST_CSV     <- ""
OUT_DIR          <- ""

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STAGE_COLORS <- c(
  Control = "#2166ac", Stage1 = "#74add1", Stage2 = "#fee090",
  Stage3  = "#fdae61", Stage4 = "#d73027", Stage5 = "#a50026"
)
GROUP_ORDER <- c("Control", "Stage1", "Stage2", "Stage3", "Stage4", "Stage5")
STAGE_ORDER <- c("Stage1", "Stage2", "Stage3", "Stage4", "Stage5")

# ============================================================================
# 3. LOAD DATA
# ============================================================================

mapping <- read.csv(SAMPLE_CSV, stringsAsFactors = FALSE)
mapping$Basename <- file.path(IDAT_DIR, mapping$SampleID)
mapping$Age      <- as.numeric(mapping$Age)

if (!"Age" %in% colnames(mapping))
  stop("CSV must have 'Age' column")

RGSet <- read.metharray.exp(targets = mapping, extended = TRUE, force = TRUE)
RGSet@annotation <- c(array      = "IlluminaHumanMethylationEPICv2",
                       annotation = "20a1.hg38")

cat("Loaded samples:", ncol(RGSet), "\n")
print(table(mapping$Group))

# ============================================================================
# 4. QUALITY CONTROL
# ============================================================================

detP <- detectionP(RGSet)

# --- Sample-level QC ---
failed_samples <- colMeans(detP) >  0.05
if (any(failed_samples)) {
  message("Removing failed samples: ",
          paste(colnames(RGSet)[failed_samples], collapse = ", "))
  RGSet   <- RGSet[, !failed_samples]
  mapping <- mapping[!failed_samples, ]
  detP    <- detP[, !failed_samples]
}

# --- Raw QC plot ---
mSet_qc <- preprocessRaw(RGSet)
qc      <- getQC(mSet_qc)
png(file.path(OUT_DIR, "QC_Plot.png"), width = 800, height = 600, res = 600)
plotQC(qc)
dev.off()

cat("Samples after QC:", ncol(RGSet), "\n")

# ============================================================================
# 5. NORMALIZATION
# ============================================================================

# preprocessFunnorm - recommended for multi-group studies (>2 groups)

mSet_norm <- preprocessFunnorm(RGSet)
cat("Normalization complete. Probes:", nrow(mSet_norm), "\n")

# ============================================================================
# 6. PROBE FILTERING
# ============================================================================

# --- Detection p-value filter (strict: ALL samples must pass) ---
keep_detP <- rowSums(detP[rownames(mSet_norm), ] < 0.01) == ncol(mSet_norm)
mSet_detP <- mSet_norm[keep_detP, ]
cat("After detP filter:", nrow(mSet_detP), "probes\n")

# --- SNP filter (both CpG position AND SBE position, MAF > 5%) ---
mSet_snp <- dropLociWithSnps(mSet_detP, snps = c("SBE", "CpG"), maf = 0.05)
cat("After SNP filter:", nrow(mSet_snp), "probes\n")

# --- Sex chromosome removal ---
ann <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
sex_probes <- ann$Name[ann$chr %in% c("chrX", "chrY")]
keep_auto  <- !(featureNames(mSet_snp) %in% sex_probes)
mSet_final <- mSet_snp[keep_auto, ]
cat("After sex chr removal:", nrow(mSet_final), "probes\n")

# --- QC summary ---
qc_summary <- data.frame(
  Step = c("Raw probes", "After sample QC", "After normalization",
           "After detP filter", "After SNP filter", "After sex chr removal"),
  Count = c(nrow(mSet_qc), ncol(RGSet), nrow(mSet_norm),
            nrow(mSet_detP), nrow(mSet_snp), nrow(mSet_final)),
  Note  = c("Start", paste0("Removed ", sum(failed_samples), " samples"),
            "preprocessFunnorm", "Strict: all samples pass",
            "SBE + CpG, MAF>5%", "chrX + chrY removed")
)
write.csv(qc_summary, file.path(OUT_DIR, "QC_Summary.csv"), row.names = FALSE)
print(qc_summary)

# Cross-reactive probes removal from pidsley2024.csv 

PIDSLEY_CSV <- "/mnt/home1/Catherine_Optimization/pidsley2024.csv"
pid_manifest <- read.csv(PIDSLEY_CSV, stringsAsFactors = FALSE,
                           row.names = 1)  

cat("\n=== Cross-Reactive Probe Diagnostic ===\n")
cat("Total probes in manifest:           ", nrow(pid_manifest), "\n")


offtarget_table <- table(pid_manifest$Num_offtargets == 0)
cat("Probes with Num_offtargets == 0 (clean):    ",
    sum(pid_manifest$Num_offtargets == 0, na.rm = TRUE), "\n")
cat("Probes with Num_offtargets  > 0 (cross-reactive): ",
    sum(pid_manifest$Num_offtargets  > 0, na.rm = TRUE), "\n")


cat("Probes with missing position (MissingPos == Y): ",
    sum(pid_manifest$MissingPos == "Y", na.rm = TRUE), "\n")

overlap <- sum(featureNames(mSet_final) %in% rownames(pid_manifest))
cat("Your probes matching manifest:      ", overlap,
    "of", nrow(mSet_final), "\n\n")

# BUILD PROBE REMOVAL LISTS

cross_reactive <- rownames(pid_manifest)[
  !is.na(pid_manifest$Num_offtargets) &
  pid_manifest$Num_offtargets > 0
]
cat("Cross-reactive probes identified:   ", length(cross_reactive), "\n")

# --- Non-mapping probes (MissingPos == "Y") ---

non_mapping <- rownames(pid_manifest)[
  !is.na(pid_manifest$MissingPos) &
  pid_manifest$MissingPos == "Y"
]
cat("Non-mapping probes (chr0):           ", length(non_mapping), "\n")

# --- Combine: all probes to remove ---
probes_to_remove <- unique(c(cross_reactive, non_mapping))
cat("Total probes flagged for removal:   ", length(probes_to_remove), "\n\n")

# APPLY FILTER TO mSet

probes_before_cr <- nrow(mSet_final)

all_probe_names <- featureNames(mSet_final)

keep_cr <- !all_probe_names %in% probes_to_remove

removed_probe_names <- all_probe_names[!keep_cr]

mSet_final <- mSet_final[keep_cr, ]

cat("=== Cross-Reactive Probe Removal Summary ===\n")
cat("Probes before filter:  ", probes_before_cr, "\n")
cat("Probes removed:        ", length(removed_probe_names), "\n")
cat("  of which cross-reactive (Num_offtargets > 0):",
    sum(cross_reactive %in% removed_probe_names), "\n")
cat("  of which non-mapping  (chr0/MissingPos == Y):",
    sum(non_mapping %in% removed_probe_names), "\n")
cat("Probes remaining:      ", nrow(mSet_final), "\n\n")

# UPDATE QC SUMMARY TABLE

qc_summary <- rbind(
  qc_summary,
  data.frame(
    Step = c(
      "After cross-reactive/non-mapping probe removal"
    ),

    Count = c(
      nrow(mSet_final)
    ),

    Note = c(
      paste0(
        "Removed ",
        length(removed_probe_names),
        " probes using Peters et al. 2024 EPICv2 manifest"
      )
    ),

    stringsAsFactors = FALSE
  )
)

write.csv(qc_summary, file.path(OUT_DIR, "QC_Summary_complete.csv"), row.names = FALSE)
