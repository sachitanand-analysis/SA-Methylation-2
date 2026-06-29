# ============================================================================
# 1. MAIN PER-STAGE ANALYSIS LOOP
# ============================================================================

cat("\n========= Starting per-stage analysis =========\n")

lambda_summary <- data.frame()
all_comparisons <- names(dmps_list)

for (cname in all_comparisons) {
  stage_name <- cname
  cat("\n--- Processing:", stage_name, "---\n")

  stage_dir <- file.path(OUT_DIR, stage_name)
  dir.create(stage_dir, showWarnings=FALSE, recursive=TRUE)

  all_res <- dmps_list[[cname]]
  sig_res <- subset(
  all_res,
  adj.P.Val < 0.05 &
  abs(logFC) > 0.2
)
  cat("  Significant DMPs:", nrow(sig_res), "\n")

  all_ann  <- merge(all_res, ann_sub, by.x="CpG", by.y="Name", all.x=TRUE)
  all_ann  <- add_regulation(all_ann)
  all_ann  <- add_annotation_confidence(all_ann)
  sig_ann <- subset(
  all_ann,
  adj.P.Val < 0.05 &
  abs(logFC) > 0.2
)

  write.csv(all_ann, file.path(stage_dir, paste0("All_DMPs_annotated_", stage_name, ".csv")),
            row.names=FALSE)
  write.csv(sig_ann, file.path(stage_dir, paste0("Sig_DMPs_annotated_", stage_name, ".csv")),
            row.names=FALSE)

  if (nrow(sig_ann) > 0) {
    reg_sum <- sig_ann %>% count(Regulation, Direction) %>% filter(Regulation != "Unclear")
    write.csv(reg_sum,
              file.path(stage_dir, paste0("Regulation_Summary_", stage_name, ".csv")),
              row.names=FALSE)
    cat("  Predicted Activation:", sum(sig_ann$Regulation=="Predicted Activation"),
        "| Predicted Silencing:", sum(sig_ann$Regulation=="Predicted Silencing"), "\n")
  }

  lam <- plot_qq_stage(all_res, stage_name, stage_dir)
  lambda_summary <- rbind(lambda_summary,
                           data.frame(Comparison=cname, Lambda=round(lam,3)))
  cat("  Lambda:", round(lam, 3), "\n")

  if (nrow(sig_ann) < 2) {
    cat("  Fewer than 2 significant DMPs - skipping plots\n")
    next
  }
   
   stage_label <- strsplit(cname, "_vs_")[[1]][1]
   
  plot_volcano_stage(all_ann, sig_ann, stage_name, stage_dir)
  cat("  Volcano saved\n")

  plot_heatmap_stage(
  sig_ann,
  stage_name,
  stage_dir,
  group_col = strsplit(cname, "_vs_")[[1]][2],
  stage_col = stage_label
)

  cat("  Heatmap saved\n")

  plot_violin_stage(sig_ann, stage_name, stage_dir)
  cat("  Violin saved\n")

  plot_lollipop_stage(sig_ann, stage_name, stage_dir)
  cat("  Lollipop saved\n")

  plot_chromosomal_stage(sig_ann, stage_name, stage_dir)
  cat("  Chromosomal distribution saved\n")

  plot_cpg_proportions_stage(sig_ann, stage_name, stage_dir)
  cat("  CpG proportions saved\n")

  # Pathway analysis
  run_pathway_stage(sig_ann, all_ann, stage_name, stage_dir)
  cat("  Pathway analysis done\n")

}

write.csv(lambda_summary, file.path(OUT_DIR, "Lambda_Summary_All_Stages.csv"),
          row.names=FALSE)
          
# ============================================================================
# 2. OVERALL CKD vs CONTROL 
# ============================================================================

cat("\n--- Overall CKD vs Control analysis ---\n")

ckd_dir <- file.path(OUT_DIR, "Overall_CKD_vs_Control")
dir.create(ckd_dir, showWarnings=FALSE, recursive=TRUE)

ckd_all_ann <- add_regulation(ckd_all_ann)
ckd_all_ann <- add_annotation_confidence(ckd_all_ann)
ckd_sig_ann <- subset(
  ckd_all_ann,
  adj.P.Val < 0.05 &
  abs(logFC) > 0.2
)

write.csv(ckd_sig_ann,
          file.path(ckd_dir, "Sig_DMPs_Overall_CKD_annotated.csv"),
          row.names=FALSE)

if (nrow(ckd_sig_ann) >= 2) {
  plot_volcano_stage(ckd_all_ann, ckd_sig_ann, "Overall_CKD", ckd_dir)
  plot_heatmap_stage(ckd_sig_ann, "Overall_CKD", ckd_dir, stage_col="CKD")
  plot_violin_stage(ckd_sig_ann, "Overall_CKD", ckd_dir)
  plot_lollipop_stage(ckd_sig_ann, "Overall_CKD", ckd_dir)
  plot_chromosomal_stage(ckd_sig_ann, "Overall_CKD", ckd_dir)
  plot_cpg_proportions_stage(ckd_sig_ann, "Overall_CKD", ckd_dir)
  run_pathway_stage(ckd_sig_ann, ckd_all_ann, "Overall_CKD", ckd_dir)
}
