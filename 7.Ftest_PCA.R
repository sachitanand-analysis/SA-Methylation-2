# ============================================================================
# 1. GLOBAL HEATMAP (top globally significant CpGs - F-test)
# ============================================================================

if (nrow(globally_sig) >= 2) {
  n_global <- min(500, nrow(globally_sig))
  top_g    <- rownames(globally_sig)[1:n_global]
  heat_g   <- beta_vals[intersect(top_g, rownames(beta_vals)),]
  heat_g_sc <- t(scale(t(heat_g)))
  heat_g_sc[heat_g_sc>3] <- 3; heat_g_sc[heat_g_sc< -3] <- -3

  ann_all_col <- data.frame(Group=mapping$Group, row.names=colnames(heat_g_sc))
  grp_colors  <- list(Group=STAGE_COLORS)

  png(file.path(OUT_DIR, "Global_Heatmap_Top500.png"), width=4800, height=3600, res=600)
  tryCatch(pheatmap(heat_g_sc, annotation_col=ann_all_col,
                    annotation_colors=grp_colors,
                    show_rownames=FALSE, fontsize_col=8,
                    clustering_method="ward.D2",
                    color=colorRampPalette(c("#313695","white","#a50026"))(100),
                    ),
           finally=dev.off())
}

# ============================================================================
# 2. PCA PLOTS 
# ============================================================================

pca_raw <- prcomp(t(beta_vals), scale.=FALSE)
pct_raw <- round(100 * pca_raw$sdev^2 / sum(pca_raw$sdev^2), 1)

pca_df <- data.frame(PC1=pca_raw$x[,1], PC2=pca_raw$x[,2], Group=mapping$Group)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color=Group)) +
  geom_point(size=4, alpha=0.88, stroke=0) +
  scale_color_manual(values=STAGE_COLORS, name=NULL) +
  labs(x=paste0("PC1 (",pct_raw[1],"%)"), y=paste0("PC2 (",pct_raw[2],"%)")) +
  theme_bw(base_size=13)
ggsave(file.path(OUT_DIR, "PCA_Beta.png"), p_pca, width=8, height=6, dpi=600)
