# ============================================================================
# 1. GLOBAL F-TEST
# ============================================================================

global_results <- topTable(fit2, coef = NULL, number = Inf,
                            adjust.method = "BH", sort.by = "F")
globally_sig <- global_results[global_results$adj.P.Val < 0.05, ]
cat("Globally significant CpGs (F-test, FDR<0.05):", nrow(globally_sig), "\n")

# Annotate global results
gs_cpg  <- data.frame(CpG = rownames(globally_sig), globally_sig)
gs_ann  <- merge(gs_cpg, ann_sub, by.x = "CpG", by.y = "Name", all.x = TRUE)
write.csv(gs_ann, file.path(OUT_DIR, "Global_ANOVA_significant_annotated.csv"),
          row.names = FALSE)

# ============================================================================
# 2. OVERALL CKD vs CONTROL MODEL (all stages as one CKD group)
# ============================================================================

mapping$CKD_Group <- factor(
  ifelse(mapping$Group == "Control", "Control", "CKD"),
  levels = c("Control", "CKD")
)
keep_sv_ckd <- intersect(colnames(sv_matrix), paste0("SV", seq_len(n_sv_use)))

if (!is.null(sv_matrix)) {

  design_ckd <- cbind(
    model.matrix(~ 0 + CKD_Group,
                 data = mapping),
    age,
    sv_matrix[, keep_sv_ckd, drop=FALSE]
  )

} else {

  design_ckd <- cbind(
    model.matrix(~ 0 + CKD_Group,
                 data = mapping),
    age
  )
}

colnames(design_ckd)[1:2] <- c(
  "Control",
  "CKD"
)

cont_ckd  <- makeContrasts(CKD_vs_Control = CKD - Control, levels = design_ckd)
fit_ckd   <- lmFit(m_vals, design_ckd)
fit_ckd2  <- contrasts.fit(fit_ckd, cont_ckd)
fit_ckd2  <- eBayes(fit_ckd2)
ckd_res   <- topTable(fit_ckd2, coef = 1, number = Inf, adjust.method = "BH")


# ============================================================
# 3. QQ PLOT + GENOMIC INFLATION (LAMBDA)
# Overall CKD vs Control
# ============================================================

pvals_ckd <- ckd_res$P.Value
pvals_ckd <- pvals_ckd[!is.na(pvals_ckd) & pvals_ckd > 0]

# --- QQ Plot ---
png(file.path(OUT_DIR, "QQplot_Overall_CKD.png"),
    width = 4800, height = 3600, res = 600)

expected <- -log10(ppoints(length(pvals_ckd)))
observed <- -log10(sort(pvals_ckd))

plot(expected, observed,
     pch = 16,
     cex = 0.5,
     xlab = "Expected -log10(p)",
     ylab = "Observed -log10(p)")

abline(0, 1, col = "red", lwd = 2)

dev.off()

# --- Genomic Inflation Factor (lambda) ---
lambda_ckd <- median(qchisq(1 - pvals_ckd, df = 1)) /
              qchisq(0.5, df = 1)

cat("\n====================================\n")
cat("Overall CKD vs Control Lambda:\n")
cat("Lambda =", round(lambda_ckd, 3), "\n")
cat("====================================\n")

ckd_res$CpG <- rownames(ckd_res)

ckd_res$Direction  <- ifelse(ckd_res$logFC > 0, "Hyper",
                      ifelse(ckd_res$logFC < 0, "Hypo", "NoChange"))

# Annotate overall CKD results
ckd_all_ann <- merge(ckd_res, ann_sub, by.x = "CpG", by.y = "Name", all.x = TRUE)
write.csv(ckd_all_ann, file.path(OUT_DIR, "Overall_CKD_vs_Control_annotated.csv"),
          row.names = FALSE)
cat("Overall CKD vs Control DMPs:", nrow(subset(ckd_all_ann, adj.P.Val < 0.05 & abs(logFC) > 0.2 )), "\n")

# ============================================================================
# 4. EXTRACT ALL PAIRWISE RESULTS AND BUILD SUMMARY TABLE
# ============================================================================

dmps_list <- list()

for (i in seq_len(ncol(cont.matrix))) {

  cname <- colnames(cont.matrix)[i]

  res <- topTable(
    fit2,
    coef = i,
    number = Inf,
    adjust.method = "BH"
  )

  res$CpG <- rownames(res)

  res$Comparison <- cname
  parts <- strsplit(cname, "_vs_")[[1]]

cols_a <- colnames(beta_vals)[
  mapping$Group == parts[1]
]

cols_b <- colnames(beta_vals)[
  mapping$Group == parts[2]
]

if (length(cols_a) > 0 &&
    length(cols_b) > 0) {

  db <- rowMeans(
    beta_vals[, cols_a, drop = FALSE]
  ) -
  rowMeans(
    beta_vals[, cols_b, drop = FALSE]
  )

  res$delta_beta <- db[
    match(
      res$CpG,
      names(db)
    )
  ]
}
  # -----------------------------------
  # HYPER / HYPO LABEL
  # -----------------------------------

  res$Direction <- ifelse(
    res$logFC > 0,
    "Hyper",
    ifelse(
      res$logFC < 0,
      "Hypo",
      "NoChange"
    )
  )

  # -----------------------------------
  # STORE RESULTS
  # -----------------------------------

  dmps_list[[cname]] <- res

  cat(
    "Stored:",
    cname,
    "| Rows:",
    nrow(res),
    "\n"
  )
}

summary_df <- data.frame(
  Comparison         = names(dmps_list),
  Sig_FDR05          = sapply(dmps_list, function(x) sum(x$adj.P.Val < 0.05)),
  Sig_FDR05_logFC02  = sapply(dmps_list, function(x)
    sum(x$adj.P.Val < 0.05 & abs(x$logFC) > 0.2)),
  Hyper              = sapply(dmps_list, function(x)
    sum(x$adj.P.Val < 0.05 & x$logFC > 0.2)),
  Hypo               = sapply(dmps_list, function(x)
    sum(x$adj.P.Val < 0.05 & x$logFC < -0.2))
)

# Parse comparison names into Group A vs Group B
parse_comparison <- function(cname) {
  parts <- strsplit(cname, "_vs_")[[1]]
  data.frame(GroupA = parts[1], GroupB = parts[2])
}
parsed <- do.call(rbind, lapply(summary_df$Comparison, parse_comparison))
summary_df <- cbind(summary_df, parsed)

# All unique groups
all_groups <- GROUP_ORDER

print(summary_df)
write.csv(summary_df, file.path(OUT_DIR, "Summary_DMPs_per_comparison.csv"),
          row.names = FALSE)

ckd_sig <- subset(
  ckd_res,
  adj.P.Val < 0.05 &
  abs(logFC) > 0.2
)

# Save to Excel
wb <- createWorkbook()
addWorksheet(wb, "Summary");     writeData(wb, "Summary", summary_df)
addWorksheet(wb, "Global_Ftest"); writeData(wb, "Global_Ftest", globally_sig)
addWorksheet(wb, "Overall_CKD")
writeData(wb, "Overall_CKD", ckd_sig)
for (cname in names(dmps_list)) {
  sig <- subset(
  dmps_list[[cname]],
  adj.P.Val < 0.05 &
  abs(logFC) > 0.2
)
  if (nrow(sig) > 0) {
    sheet_name <- gsub("[\\[\\]\\*\\/\\\\\\?:]", "_", substr(cname, 1, 31))

addWorksheet(
  wb,
  sheet_name
)

writeData(
  wb,
  sheet_name,
  sig
   )
  }
}
saveWorkbook(wb, file.path(OUT_DIR, "All_DMPs_SVAcorrected.xlsx"), overwrite=TRUE)
