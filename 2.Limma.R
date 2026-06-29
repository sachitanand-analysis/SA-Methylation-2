# ============================================================================
# 1. EXTRACT M-VALUES AND BETA VALUES
# ============================================================================

beta_vals <- getBeta(mSet_final)
m_vals    <- getM(mSet_final)

# Align sample order (critical after QC filtering)
mapping <- mapping[match(colnames(m_vals), mapping$SampleID), ]
cat("Sample alignment check:", all(colnames(m_vals) == mapping$SampleID), "\n")

# Beta density plot
png(file.path(OUT_DIR, "Beta_Density.png"), width = 4800, height = 3600, res = 600)
densityPlot(beta_vals, sampGroups = mapping$Group,
            main = "Beta Value Distribution - Post-Normalization")
dev.off()

# ============================================================================
# 2. DEFINE GROUPS AND COVARIATES
# ============================================================================

group <- factor(mapping$Group, levels = GROUP_ORDER)
age   <- mapping$Age

# ============================================================================
# 3. SVA - SURROGATE VARIABLE ANALYSIS
# ============================================================================

mod_sva  <- model.matrix(~ group + age,
                          data = data.frame(group = group, age = age))
mod0_sva <- model.matrix(~ age,
                          data = data.frame(age = age))

n_sv_leek <- num.sv(m_vals, mod = mod_sva, method = "leek")
cat("SVs estimated (leek, single-thread):", n_sv_leek, "\n")

# Degrees of freedom safety check
n_samples   <- ncol(m_vals)
max_safe_sv <- n_samples - (nlevels(group) + 1) - 2
n_sv_use    <- min(n_sv_leek, max_safe_sv)
cat("SVs used:", n_sv_use, "\n")

if (n_sv_use > 0) {
  sva_obj   <- sva(m_vals, mod = mod_sva, mod0 = mod0_sva, n.sv = n_sv_use)
  sv_matrix <- sva_obj$sv
  colnames(sv_matrix) <- paste0("SV", seq_len(n_sv_use))
  sv_df <- as.data.frame(sv_matrix)
  sv_df$SampleID <- mapping$SampleID
  write.csv(sv_df, file.path(OUT_DIR, "Surrogate_Variables.csv"), row.names = FALSE)

  # SV correlation plot
  known_vars <- data.frame(Group_numeric = as.numeric(group), Age = age)
  sv_cor <- cor(sv_matrix, known_vars, method = "spearman")
  png(file.path(OUT_DIR, "SV_Correlation_Heatmap.png"), width = 800, height = 500)
  corrplot(t(sv_cor), method = "color", tl.cex = 0.9,
           title = "SV Correlations with Known Variables", mar = c(0,0,2,0))
  dev.off()
} else {
  sv_matrix <- NULL
  cat("No SVs estimated - proceeding without SVA correction\n")
}

# ============================================================================
# 4. MANIFEST ANNOTATION (mirrors final_github.R Sections 8-10)
# Uses the raw EPIC v2 manifest CSV for precise N/S Shore/Shelf assignment
# ============================================================================

cat("Loading manifest for refined CpG annotation...\n")

manifest <- read.csv(MANIFEST_CSV, skip = 7, stringsAsFactors = FALSE)

assign_primary_context <- function(x) {
  if (is.na(x) || x == "") return("OpenSea")
  parts <- unique(trimws(unlist(strsplit(x, ";"))))
  if ("Island"  %in% parts) return("Island")
  if ("N_Shore" %in% parts) return("N_Shore")
  if ("S_Shore" %in% parts) return("S_Shore")
  if ("N_Shelf" %in% parts) return("N_Shelf")
  if ("S_Shelf" %in% parts) return("S_Shelf")
  return("OpenSea")
}

assign_simple_context <- function(x) {
  if (is.na(x) || x == "") return("OpenSea")
  parts <- unique(trimws(unlist(strsplit(x, ";"))))
  if ("Island" %in% parts)                     return("Island")
  if (any(c("N_Shore","S_Shore") %in% parts))  return("Shore")
  if (any(c("N_Shelf","S_Shelf") %in% parts))  return("Shelf")
  return("OpenSea")
}

is_multi_context <- function(x) {
  if (is.na(x) || x == "") return(FALSE)
  parts <- unique(trimws(unlist(strsplit(x, ";"))))
  sum(c("Island" %in% parts,
        any(c("N_Shore","S_Shore") %in% parts),
        any(c("N_Shelf","S_Shelf") %in% parts))) > 1
}

raw_vals <- manifest$Relation_to_UCSC_CpG_Island
manifest$Relation_Primary  <- sapply(raw_vals, assign_primary_context)
manifest$Relation_Simple   <- sapply(raw_vals, assign_simple_context)
manifest$Is_Multi_Context  <- sapply(raw_vals, is_multi_context)

idx <- match(rownames(ann), manifest$IlmnID)
ann$Relation_Primary  <- manifest$Relation_Primary[idx]
ann$Relation_Simple   <- manifest$Relation_Simple[idx]
ann$Is_Multi_Context  <- manifest$Is_Multi_Context[idx]

# Fill unmatched probes from existing annotation
na_idx <- is.na(ann$Relation_Primary)
ann$Relation_Primary[na_idx] <- as.character(ann$Relation_to_Island[na_idx])
ann$Relation_Simple[na_idx]  <- as.character(ann$Relation_to_Island[na_idx])
ann$Is_Multi_Context[na_idx] <- FALSE

cat("Manifest annotation complete\n")

# Full annotation columns to use for DMPs
ann_cols <- c("chr","pos","strand","Name","UCSC_RefGene_Name",
              "UCSC_RefGene_Group","Relation_to_Island","Relation_Primary",
              "Relation_Simple","Is_Multi_Context","Islands_Name",
              "Regulatory_Feature_Group")
ann_cols <- ann_cols[ann_cols %in% colnames(ann)]
ann_sub  <- as.data.frame(ann)[, ann_cols]

# ============================================================================
# 5. LIMMA DESIGN MATRIX + CONTRASTS
# ============================================================================

# Build design (intercept-free, with age and SVs as covariates)
if (!is.null(sv_matrix)) {
  design <- cbind(model.matrix(~ 0 + group),
                  age,
                  sv_matrix)
} else {
  design <- cbind(model.matrix(~ 0 + group), age)
}
colnames(design)[1:nlevels(group)] <- levels(group)

cat("Design matrix columns:", colnames(design), "\n")
cat("Residual df:", nrow(design) - ncol(design), "\n")

# All stage vs control + all pairwise contrasts
cont.matrix <- makeContrasts(
  Stage1_vs_Control = Stage1 - Control,
  Stage2_vs_Control = Stage2 - Control,
  Stage3_vs_Control = Stage3 - Control,
  Stage4_vs_Control = Stage4 - Control,
  Stage5_vs_Control = Stage5 - Control,
  Stage1_vs_Stage2  = Stage1 - Stage2,
  Stage1_vs_Stage3  = Stage1 - Stage3,
  Stage1_vs_Stage4  = Stage1 - Stage4,
  Stage1_vs_Stage5  = Stage1 - Stage5,
  Stage2_vs_Stage1  = Stage2 - Stage1,
  Stage2_vs_Stage3  = Stage2 - Stage3,
  Stage2_vs_Stage4  = Stage2 - Stage4,
  Stage2_vs_Stage5  = Stage2 - Stage5,
  Stage3_vs_Stage1  = Stage3 - Stage1,
  Stage3_vs_Stage2  = Stage3 - Stage2,
  Stage3_vs_Stage4  = Stage3 - Stage4,
  Stage3_vs_Stage5  = Stage3 - Stage5,
  Stage4_vs_Stage1  = Stage4 - Stage1,
  Stage4_vs_Stage2  = Stage4 - Stage2,
  Stage4_vs_Stage3  = Stage4 - Stage3,
  Stage4_vs_Stage5  = Stage4 - Stage5,
  Stage5_vs_Stage1  = Stage5 - Stage1,
  Stage5_vs_Stage2  = Stage5 - Stage2,
  Stage5_vs_Stage3  = Stage5 - Stage3,
  Stage5_vs_Stage4  = Stage5 - Stage4,
  levels = design
)

# ============================================================================
# 6. FIT MODEL
# ============================================================================

fit  <- lmFit(m_vals, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
cat("Model fitting complete\n")
