# ============================================================
# 1. PROGRESSIVE CKD SIGNATURE
# ============================================================

# READ OVERALL CKD DMPs

overall_dmp <- read.csv(
  file.path(
    DMP_DIR,
    "Sig_DMPs_annotated_Overall_CKD.csv"
  ),
  stringsAsFactors = FALSE
)

overall_dmp <- overall_dmp %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 0.2
  )

candidate_cpgs <- unique(
  overall_dmp$CpG
)

cat(
  "Overall CKD DMPs:",
  length(candidate_cpgs),
  "\n"
)

# CREATE CKD SEVERITY SCORE

mapping$Severity <- dplyr::case_when(
  mapping$Group == "Control" ~ 0,
  mapping$Group == "Stage1"  ~ 1,
  mapping$Group == "Stage2"  ~ 2,
  mapping$Group == "Stage3"  ~ 3,
  mapping$Group == "Stage4"  ~ 4,
  mapping$Group == "Stage5"  ~ 5
)

severity <- mapping$Severity

stopifnot(
  length(severity) == ncol(beta_vals)
)

# PEARSON CORRELATION

progression_list <- lapply(
  candidate_cpgs,
  function(cpg){

    if(!cpg %in% rownames(beta_vals))
      return(NULL)

    beta <- as.numeric(
      beta_vals[cpg, ]
    )

    if(sd(beta, na.rm = TRUE) == 0)
      return(NULL)

    zbeta <- as.numeric(
      scale(beta)
    )

    r <- cor(
      zbeta,
      severity,
      method = "pearson"
    )

    n <- length(severity)

    tstat <- r * sqrt(
      (n - 2) / (1 - r^2)
    )

    pval <- 2 * pt(
      -abs(tstat),
      df = n - 2
    )

    data.frame(
      CpG = cpg,
      r   = r,
      p   = pval
    )
  }
)

progression_df <- bind_rows(
  progression_list
)

# BH CORRECTION

progression_df$FDR <- p.adjust(
  progression_df$p,
  method = "BH"
)

progression_df$Direction <- ifelse(
  progression_df$r > 0,
  "Progressive_Gain",
  "Progressive_Loss"
)

# ADD GENE

progression_df <- merge(
  progression_df,
  overall_dmp[
    ,
    c(
      "CpG",
      "UCSC_RefGene_Name",
      "logFC",
      "adj.P.Val"
    )
  ],
  by = "CpG",
  all.x = TRUE
)

colnames(progression_df)[
  colnames(progression_df) == "UCSC_RefGene_Name"
] <- "Gene"

# ADD STAGE-SPECIFIC MEAN BETA VALUES

stage_levels <- c(
  "Control",
  "Stage1",
  "Stage2",
  "Stage3",
  "Stage4",
  "Stage5"
)

stage_means_list <- lapply(
  progression_df$CpG,
  function(cpg){

    vals <- as.numeric(
      beta_vals[cpg, ]
    )

    means <- sapply(
      stage_levels,
      function(g){

        mean(
          vals[mapping$Group == g],
          na.rm = TRUE
        )
      }
    )

    data.frame(
      CpG     = cpg,
      Control = means[1],
      Stage1  = means[2],
      Stage2  = means[3],
      Stage3  = means[4],
      Stage4  = means[5],
      Stage5  = means[6]
    )
  }
)

stage_means_df <- bind_rows(
  stage_means_list
)

progression_df <- merge(
  progression_df,
  stage_means_df,
  by = "CpG",
  all.x = TRUE
)

progression_df$Delta_Stage5_Control <-
  progression_df$Stage5 -
  progression_df$Control

progression_df <- progression_df %>%
  dplyr::select(
    CpG,
    Gene,
    logFC,
    adj.P.Val,
    Direction,
    r,
    p,
    FDR,
    Control,
    Stage1,
    Stage2,
    Stage3,
    Stage4,
    Stage5,
    Delta_Stage5_Control
  )

progression_df <- progression_df %>%
  arrange(FDR)

# SAVE ALL RESULTS

write.csv(
  progression_df,
  file.path(
    OUT_DIR,
    "OverallCKD_Progressive_WithStageMeans.csv"
  ),
  row.names = FALSE
)

# SIGNIFICANT PROGRESSIVE CpGs

sig_progressive <- progression_df %>%
  filter(
    FDR < 0.05
  )

write.csv(
  sig_progressive,
  file.path(
    OUT_DIR,
    "OverallCKD_Progressive_Pearson_FDR05.csv"
  ),
  row.names = FALSE
)

# STRONG CORRELATION CpGs

strong_progressive <- progression_df %>%
  filter(
    abs(r) >= 0.7
  )

write.csv(
  strong_progressive,
  file.path(
    OUT_DIR,
    "OverallCKD_Progressive_Pearson_r0.7.csv"
  ),
  row.names = FALSE
)

# SUMMARY

cat("\n====================================\n")
cat("Progressive CKD CpGs\n")
cat("====================================\n")

cat(
  "Total CpGs:",
  nrow(progression_df),
  "\n"
)

cat(
  "FDR < 0.05:",
  nrow(sig_progressive),
  "\n"
)

cat(
  "|r| >= 0.7:",
  nrow(strong_progressive),
  "\n"
)

cat(
  "Gain:",
  sum(strong_progressive$r > 0),
  "\n"
)

cat(
  "Loss:",
  sum(strong_progressive$r < 0),
  "\n"
)

# TOP 20

top20 <- progression_df %>%
  slice(1:20)

write.csv(
  top20,
  file.path(
    OUT_DIR,
    "Top20_overall_Progressive_CpGs.csv"
  ),
  row.names = FALSE
)

print(top20)

# ============================================================
# 2. BIOMARKER
# ============================================================

prog <- read.csv(
  file.path(
    OUT_DIR,
    "OverallCKD_Progressive_WithStageMeans.csv"
  ),
  stringsAsFactors = FALSE
)

top_cpgs <- prog %>%
  filter(FDR < 0.05) %>%
  arrange(desc(abs(r))) %>%
  slice(1:20)

top_cpgs[, c("CpG","Gene","r","FDR")]

roc_df <- data.frame(
  Group = mapping$Group
)

roc_df$CKD_Status <- ifelse(
  roc_df$Group == "Control",
  "Control",
  "CKD"
)

roc_df$CKD_Status <- factor(
  roc_df$CKD_Status,
  levels = c("Control","CKD")
)

auc_results <- lapply(
  top_cpgs$CpG,
  function(cpg){

    vals <- as.numeric(
      beta_vals[cpg, ]
    )

    roc_obj <- roc(
      response = roc_df$CKD_Status,
      predictor = vals,
      levels = c("Control","CKD"),
      quiet = TRUE
    )

    data.frame(
      CpG = cpg,
      AUC = as.numeric(
        auc(roc_obj)
      )
    )
  }
)

auc_results <- bind_rows(
  auc_results
)

auc_results <- merge(
  auc_results,
  top_cpgs[
    ,
    c(
      "CpG",
      "Gene",
      "r",
      "FDR"
    )
  ],
  by = "CpG"
)

auc_results <- auc_results %>%
  arrange(desc(AUC))
  
write.csv(
  auc_results,
  file.path(
    OUT_DIR,
    "Top20_ProgressiveCpGs_AUC_Ranking.csv"
  ),
  row.names = FALSE
)

print(auc_results)

top10 <- auc_results %>%
  arrange(desc(AUC)) %>%
  slice(1:10)
  
png(
  file.path(
    OUT_DIR,
    "Top10_ProgressiveCpGs_ROC.png"
  ),
  width = 1200,
  height = 1000
)

first <- TRUE

for(cpg in top10$CpG){

  vals <- as.numeric(
    beta_vals[cpg, ]
  )

  roc_obj <- roc(
    roc_df$CKD_Status,
    vals,
    quiet = TRUE
  )

  if(first){

    plot(
      roc_obj,
      lwd = 2
    )

    first <- FALSE

  } else {

    plot(
      roc_obj,
      add = TRUE,
      lwd = 2
    )
  }
}

legend(
  "bottomright",
  legend = paste0(
    top10$Gene,
    " (",
    round(top10$AUC,3),
    ")"
  ),
  cex = 0.8
)

dev.off()

# ============================================================
# 3. READ CANDIDATE CpGs
# ============================================================

cpg_matrix <- read.csv(
  file.path(
    OUT_DIR,
    "CandidateGene_CpG_ComparisonMatrix.csv"
  ),
  stringsAsFactors = FALSE
)

candidate_cpgs <- unique(cpg_matrix$CpG)

cat(
  "Candidate CpGs:",
  length(candidate_cpgs),
  "\n"
)

# CREATE CKD SEVERITY SCORE

mapping$Severity <- dplyr::case_when(
  mapping$Group == "Control" ~ 0,
  mapping$Group == "Stage1"  ~ 1,
  mapping$Group == "Stage2"  ~ 2,
  mapping$Group == "Stage3"  ~ 3,
  mapping$Group == "Stage4"  ~ 4,
  mapping$Group == "Stage5"  ~ 5
)

severity <- mapping$Severity

stopifnot(
  length(severity) == ncol(beta_vals)
)

# PEARSON CORRELATION

progression_list <- lapply(
  candidate_cpgs,
  function(cpg){

    if(!cpg %in% rownames(beta_vals))
      return(NULL)

    beta <- as.numeric(
      beta_vals[cpg, ]
    )

    # skip invariant CpGs
    if(sd(beta, na.rm = TRUE) == 0)
      return(NULL)

    zbeta <- as.numeric(
      scale(beta)
    )

    ct <- suppressWarnings(
      cor.test(
        zbeta,
        severity,
        method = "pearson"
      )
    )

    data.frame(
      CpG = cpg,
      r   = unname(ct$estimate),
      p   = ct$p.value
    )
  }
)

progression_df <- bind_rows(
  progression_list
)

# BH CORRECTION

progression_df$FDR <- p.adjust(
  progression_df$p,
  method = "BH"
)

progression_df$Direction <- ifelse(
  progression_df$r > 0,
  "Progressive_Gain",
  "Progressive_Loss"
)

# ADD GENE

progression_df <- merge(
  progression_df,
  unique(
    cpg_matrix[, c("CpG","Gene")]
  ),
  by = "CpG",
  all.x = TRUE
)

# ADD STAGE-SPECIFIC MEAN BETA VALUES

stage_levels <- c(
  "Control",
  "Stage1",
  "Stage2",
  "Stage3",
  "Stage4",
  "Stage5"
)

stage_means_list <- lapply(
  progression_df$CpG,
  function(cpg){

    vals <- as.numeric(
      beta_vals[cpg, ]
    )

    means <- sapply(
      stage_levels,
      function(g){

        mean(
          vals[mapping$Group == g],
          na.rm = TRUE
        )
      }
    )

    data.frame(
      CpG     = cpg,
      Control = means[1],
      Stage1  = means[2],
      Stage2  = means[3],
      Stage3  = means[4],
      Stage4  = means[5],
      Stage5  = means[6]
    )
  }
)

stage_means_df <- bind_rows(
  stage_means_list
)

progression_df <- merge(
  progression_df,
  stage_means_df,
  by = "CpG",
  all.x = TRUE
)

progression_df$Delta_Stage5_Control <-
  progression_df$Stage5 -
  progression_df$Control

progression_df <- progression_df %>%
  dplyr::select(
    CpG,
    Gene,
    Direction,
    r,
    p,
    FDR,
    Control,
    Stage1,
    Stage2,
    Stage3,
    Stage4,
    Stage5,
    Delta_Stage5_Control
  )

progression_df <- progression_df %>%
  arrange(FDR)

# SAVE ALL RESULTS

write.csv(
  progression_df,
  file.path(
    OUT_DIR,
    "CandidateGene_Progressive_WithStageMeans.csv"
  ),
  row.names = FALSE
)

# SIGNIFICANT PROGRESSIVE CpGs

sig_progressive <- progression_df %>%
  filter(
    FDR < 0.05
  )

write.csv(
  sig_progressive,
  file.path(
    OUT_DIR,
    "CandidateGene_Progressive_Pearson_FDR05.csv"
  ),
  row.names = FALSE
)

# STRONG CORRELATION CpGs

strong_progressive <- progression_df %>%
  filter(
    abs(r) >= 0.7
  )

write.csv(
  strong_progressive,
  file.path(
    OUT_DIR,
    "CandidateGene_Progressive_Pearson_r0.7.csv"
  ),
  row.names = FALSE
)

# SUMMARY

cat("\n====================================\n")
cat("Progressive CKD CpGs\n")
cat("====================================\n")

cat(
  "Total CpGs:",
  nrow(progression_df),
  "\n"
)

cat(
  "FDR < 0.05:",
  nrow(sig_progressive),
  "\n"
)

cat(
  "|r| >= 0.7:",
  nrow(strong_progressive),
  "\n"
)

cat(
  "Gain:",
  sum(strong_progressive$r > 0),
  "\n"
)

cat(
  "Loss:",
  sum(strong_progressive$r < 0),
  "\n"
)

# TOP 20

top20 <- progression_df %>%
  slice(1:20)

write.csv(
  top20,
  file.path(
    OUT_DIR,
    "Top20_Progressive_CpGs.csv"
  ),
  row.names = FALSE
)

print(top20)

# ============================================================
# 4. CYP3A5 SUMMARY TABLE
# ============================================================

cpg <- "cg16879596_BC21"

stage_summary <- data.frame(
  Group = c(
    "Control",
    "Stage1",
    "Stage2",
    "Stage3",
    "Stage4",
    "Stage5"
  ),
  Mean_Beta = c(
    mean(beta_vals[cpg, mapping$Group=="Control"]),
    mean(beta_vals[cpg, mapping$Group=="Stage1"]),
    mean(beta_vals[cpg, mapping$Group=="Stage2"]),
    mean(beta_vals[cpg, mapping$Group=="Stage3"]),
    mean(beta_vals[cpg, mapping$Group=="Stage4"]),
    mean(beta_vals[cpg, mapping$Group=="Stage5"])
  ),
  SD_Beta = c(
    sd(beta_vals[cpg, mapping$Group=="Control"]),
    sd(beta_vals[cpg, mapping$Group=="Stage1"]),
    sd(beta_vals[cpg, mapping$Group=="Stage2"]),
    sd(beta_vals[cpg, mapping$Group=="Stage3"]),
    sd(beta_vals[cpg, mapping$Group=="Stage4"]),
    sd(beta_vals[cpg, mapping$Group=="Stage5"])
  )
)

write.csv(
  stage_summary,
  file.path(
    OUT_DIR,
    "CYP3A5_StageMeans.csv"
  ),
  row.names = FALSE
)

# PLOT

cpg <- "cg16879596_BC21"

plot_df <- data.frame(
  Beta = as.numeric(
    beta_vals[cpg, ]
  ),
  Group = factor(
    mapping$Group,
    levels = c(
      "Control",
      "Stage1",
      "Stage2",
      "Stage3",
      "Stage4",
      "Stage5"
    )
  )
)

p1 <- ggplot(
  plot_df,
  aes(
    Group,
    Beta
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = 1),
    linewidth = 1.2,
    colour = "red"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    colour = "red"
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "",
    y = expression(beta~value)
  )

ggsave(
  filename = file.path(
    OUT_DIR,
    "Figure_CYP3A5_Trajectory.png"
  ),
  plot = p1,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(
    OUT_DIR,
    "Figure_CYP3A5_Trajectory.pdf"
  ),
  plot = p1,
  width = 7,
  height = 5
)

# AUC plot

roc_ckd <- roc(
  roc_df$CKD_Status,
  as.numeric(
    beta_vals["cg16879596_BC21", ]
  )
)

png(
  file.path(
    OUT_DIR,
    "Figure_CYP3A5_ROC.png"
  ),
  width = 2200,
  height = 2200,
  res = 300
)

plot(
  roc_ckd,
  col = "#1F77B4",
  lwd = 4,
  main = ""
)

abline(
  a = 0,
  b = 1,
  lty = 2,
  col = "grey50"
)

legend(
  "bottomright",
  legend = paste0(
    "AUC = ",
    round(
      auc(roc_ckd),
      3
    )
  ),
  bty = "n",
  cex = 1.3
)

dev.off()

# COMBINED PLOT

library(patchwork)

roc_plot <- ggroc(
  roc_ckd,
  linewidth = 1.2
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  theme_bw(base_size = 14) +

combined_plot <- p1 + roc_plot +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  file.path(
    OUT_DIR,
    "Figure_CYP3A5_Biomarker.pdf"
  ),
  combined_plot,
  width = 12,
  height = 5
)

ggsave(
  file.path(
    OUT_DIR,
    "Figure_CYP3A5_Biomarker.png"
  ),
  combined_plot,
  width = 12,
  height = 5,
  dpi = 600
)

coords(
  roc_ckd,
  "best",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  )
)

best_cut <- as.numeric(
  coords(
    roc_ckd,
    "best",
    ret = "threshold"
  )
)

pred <- ifelse(
  roc_df$Beta >= best_cut,
  "CKD",
  "Control"
)

length(pred)
length(roc_df$CKD_Status)

table(
  Predicted = pred,
  Actual = roc_df$CKD_Status
)
