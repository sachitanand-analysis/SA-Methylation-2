# ============================================================================
# 1HELPER FUNCTIONS 
# ============================================================================

# --- Regulation interpretation ---
add_regulation <- function(df) {
  df$Regulation <- "Unclear"
  promoter_pattern <- "TSS200|TSS1500|5'UTR|1stExon|5UTR|exon_1"
  df$Regulation[df$Direction == "Hyper" &
    grepl(promoter_pattern, df$UCSC_RefGene_Group)] <- "Predicted Silencing"
  df$Regulation[df$Direction == "Hypo" &
    grepl(promoter_pattern, df$UCSC_RefGene_Group)] <- "Predicted Activation"
  df
}

# --- Annotation confidence ---
add_annotation_confidence <- function(df) {
  df$Annotation_Confidence <- "Low"
  df$Annotation_Confidence[
    grepl("TSS1500|5'UTR|1stExon|5UTR|exon_1", df$UCSC_RefGene_Group)
  ] <- "Medium"
  df$Annotation_Confidence[
    grepl("TSS200|1stExon|5'UTR|5UTR|exon_1", df$UCSC_RefGene_Group) &
    df$Relation_Primary == "Island"
  ] <- "High"
  df
}

# --- Genomic feature assignment ---
add_genomic_feature <- function(df) {
  df$Genomic_Feature <- "Other"
  df$Genomic_Feature[grepl("TSS200",  df$UCSC_RefGene_Group)] <- "TSS200"
  df$Genomic_Feature[grepl("TSS1500", df$UCSC_RefGene_Group) &
    df$Genomic_Feature == "Other"] <- "TSS1500"
  df$Genomic_Feature[grepl("5'UTR|5UTR", df$UCSC_RefGene_Group) &
    df$Genomic_Feature == "Other"] <- "5'UTR"
  df$Genomic_Feature[grepl("1stExon|exon_1", df$UCSC_RefGene_Group) &
    df$Genomic_Feature == "Other"] <- "1st Exon"
  df$Genomic_Feature[grepl("exon_[2-9]|exon_[1-9][0-9]|Body",
    df$UCSC_RefGene_Group) & df$Genomic_Feature == "Other"] <- "Body"
  df$Genomic_Feature[grepl("3'UTR|3UTR", df$UCSC_RefGene_Group) &
    df$Genomic_Feature == "Other"] <- "3'UTR"
  df$Genomic_Feature <- factor(df$Genomic_Feature,
    levels = c("TSS200","TSS1500","5'UTR","1st Exon","Body","3'UTR","Other"))
  df
}

# --- Volcano plot  ---
plot_volcano_stage <- function(all_df, sig_df, stage_name, out_dir,
                                fdr_thresh = 0.05, lfc_thresh = 0.2) {
  all_df$Category <- "Not significant"
  all_df$Category[all_df$adj.P.Val < fdr_thresh & all_df$logFC >  lfc_thresh] <- "Hypermethylated"
  all_df$Category[all_df$adj.P.Val < fdr_thresh & all_df$logFC < -lfc_thresh] <- "Hypomethylated"

  # Labels from sig genes
  ann_gene <- ann_sub[, c("Name","UCSC_RefGene_Name")]
  ann_gene$Gene <- sapply(strsplit(ann_gene$UCSC_RefGene_Name, ";"),
                           function(x) if (length(x)==0||x[1]=="") NA else x[1])
  all_df$Gene <- ann_gene$Gene[match(all_df$CpG, ann_gene$Name)]

  lh <- subset(all_df, Category=="Hypermethylated" & !is.na(Gene) & Gene!="")
  lh <- lh %>%
  arrange(adj.P.Val) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = 8)
  lp <- subset(all_df, Category=="Hypomethylated"  & !is.na(Gene) & Gene!="")
  lp <- lp %>%
  arrange(adj.P.Val) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = 8)
  label_df <- rbind(lh, lp)

  x_max <- max(abs(all_df$logFC), na.rm=TRUE)
  y_max <- max(-log10(all_df$adj.P.Val[is.finite(-log10(all_df$adj.P.Val))]),
               na.rm=TRUE)

  p <- ggplot(all_df[order(all_df$Category=="Not significant",
                             decreasing=TRUE),],
              aes(logFC, -log10(adj.P.Val), color=Category,
                  size=Category, alpha=Category)) +
    geom_point(stroke=0) +
    geom_hline(yintercept=-log10(fdr_thresh), linetype="dashed",
               linewidth=0.5, color="grey40") +
    geom_vline(xintercept=c(-lfc_thresh, lfc_thresh), linetype="dashed",
               linewidth=0.5, color="grey40") +
    {if(nrow(label_df)>0)
      geom_label_repel(data=label_df, aes(label=Gene), size=2.8,
                       fontface="italic", color="black",
                       fill=alpha("white",0.85), label.size=0.25,
                       box.padding=0.4, max.overlaps=20,
                       segment.size=0.3, show.legend=FALSE)} +
    scale_color_manual(values=c(Hypermethylated="#B2182B",
                                 Hypomethylated="#2166AC",
                                 "Not significant"="grey75"), name=NULL,
                       guide=guide_legend(override.aes=list(size=3,alpha=1))) +
    scale_size_manual(values=c(Hypermethylated=1.0,Hypomethylated=1.0,
                                "Not significant"=0.4), guide="none") +
    scale_alpha_manual(values=c(Hypermethylated=0.8,Hypomethylated=0.8,
                                 "Not significant"=0.3), guide="none") +
    scale_x_continuous(name="M-value difference (logFC)",
                       limits=c(-ceiling(x_max*10)/10-0.2,
                                 ceiling(x_max*10)/10+0.2)) +
    scale_y_continuous(name=expression(-log[10]~"(FDR-adjusted"~italic(p)~")"),
                       limits=c(0, ceiling(y_max)+0.3),
                       expand=expansion(mult=c(0,0.1))) +
    annotate("text", x=ceiling(x_max*10)/10, y=(ceiling(y_max)+0.3)*0.97,
             label=paste0("Hyper: ",sum(all_df$Category=="Hypermethylated")),
             color="#B2182B", size=3.5, fontface="bold", hjust=1) +
    annotate("text", x=-ceiling(x_max*10)/10, y=(ceiling(y_max)+0.3)*0.97,
             label=paste0("Hypo: ",sum(all_df$Category=="Hypomethylated")),
             color="#2166AC", size=3.5, fontface="bold", hjust=0) +
    theme_classic(base_size=13) +
    theme(axis.title=element_text(face="bold"),
          axis.text=element_text(color="black"),
          legend.position="top",
          panel.grid.major=element_line(color="grey93", linewidth=0.3))

  ggsave(file.path(out_dir, paste0("Volcano_", stage_name, ".png")),
         p, width=9, height=10, dpi=600)
  ggsave(file.path(out_dir, paste0("Volcano_", stage_name, ".pdf")),
         p, width=9, height=10)
}

# --- Heatmap ---
plot_heatmap_stage <- function(sig_df, stage_name, out_dir,
                                group_col="Control", stage_col,
                                n_top=100, fdr_thresh=0.05) {
  sig_sorted <- sig_df[sig_df$adj.P.Val < fdr_thresh,]
  sig_sorted <- sig_sorted[order(sig_sorted$adj.P.Val),]
  if (nrow(sig_sorted) < 2) {
    cat("Skipping heatmap for", stage_name, "- fewer than 2 sig DMPs\n"); return(NULL)
  }
  top_cpgs <- head(sig_sorted$CpG, n_top)
  top_cpgs <- intersect(top_cpgs, rownames(beta_vals))

  ctrl_cols <- colnames(beta_vals)[
  mapping$Group == group_col
]
  if (stage_col == "CKD") {

  stage_cols <- colnames(beta_vals)[mapping$CKD_Group == "CKD"]

} else {

  stage_cols <- colnames(beta_vals)[mapping$Group == stage_col]
}
  ordered    <- c(ctrl_cols, stage_cols)
  mat        <- beta_vals[top_cpgs, ordered]

  mat_sc <- t(scale(t(mat)))
  mat_sc[mat_sc > 4]  <-  4
  mat_sc[mat_sc < -4] <- -4

  annotation_col <- data.frame(
    Group = factor(
  c(
    rep(group_col, length(ctrl_cols)),
    rep(stage_col, length(stage_cols))
  )
),
    row.names = ordered
  )
  ann_colors <- list(
  Group = c(
    setNames("#2166AC", group_col),
    setNames("#D6604D", stage_col)
  )
)

  hc_rows <- hclust(dist(mat_sc, method="euclidean"), method="complete")

  png(file.path(out_dir, paste0("Heatmap_", stage_name, ".png")),
      width=4800, height=3600, res=600)
  tryCatch(pheatmap(mat_sc, cluster_rows=hc_rows, cluster_cols=FALSE,
                    treeheight_row=30, treeheight_col=0,
                    color=colorRampPalette(c("#2166AC","#F7F7F7","#D6604D"))(100),
                    breaks=seq(-4,4,length.out=101),
                    show_rownames=TRUE, show_colnames=FALSE,
                    annotation_col=annotation_col, annotation_colors=ann_colors,
                    border_color=NA, fontsize=5, annotation_legend=TRUE,
                    main=paste0("Top ", length(top_cpgs),
                                " DMPs - ", stage_name)),
           finally=dev.off())
  pdf(file.path(out_dir, paste0("Heatmap_", stage_name, ".pdf")),
      width=10, height=14)
  tryCatch(pheatmap(mat_sc, cluster_rows=hc_rows, cluster_cols=FALSE,
                    color=colorRampPalette(c("#2166AC","#F7F7F7","#D6604D"))(100),
                    breaks=seq(-4,4,length.out=101),
                    show_rownames=TRUE, show_colnames=FALSE,
                    annotation_col=annotation_col, annotation_colors=ann_colors,
                    border_color=NA, fontsize=5,
                    main=paste0("Top ", length(top_cpgs),
                                " DMPs - ", stage_name)),
           finally=dev.off())
}

# --- Violin plot ---
plot_violin_stage <- function(sig_ann, stage_name, out_dir) {
  df <- add_genomic_feature(sig_ann)
  df <- df[df$Genomic_Feature != "Other",]
  if (nrow(df) < 5) return(NULL)
  p <- ggplot(df, aes(x=Genomic_Feature, y=logFC, fill=Genomic_Feature)) +
    geom_violin(trim=TRUE, adjust=0.8, color="grey25", linewidth=0.3) +
    coord_cartesian(ylim=c(-1,1)) +
    scale_fill_manual(values=c("TSS20"="#8B3FBF","TSS1500"="#1E5AA8",
                                "5'UTR"="#006400","1st Exon"="#B22222",
                                "Body"="#008B8B","3'UTR"="#808000")) +
    theme_minimal(base_size=11) +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) +
    labs(x=NULL, y="logFC (M-value)")
  ggsave(file.path(out_dir, paste0("Violin_", stage_name, ".png")),
         p, width=8, height=10, dpi=600)
  ggsave(file.path(out_dir, paste0("Violin_", stage_name, ".pdf")),
         p, width=8, height=10)
}

# --- Lollipop plot ---
plot_lollipop_stage <- function(sig_ann, stage_name, out_dir, n_top=30) {
  df <- sig_ann %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 0.2) %>%
    mutate(
  Gene = sapply(
    strsplit(as.character(UCSC_RefGene_Name), ";"),
    `[`,
    1
  )
) %>%
    filter(!is.na(Gene) & Gene != "") %>%
    arrange(adj.P.Val) %>%
    distinct(Gene, .keep_all=TRUE) %>%
    slice_head(n=n_top) %>%
    mutate(Direction   = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated"),
           neg_log10p  = -log10(adj.P.Val),
           Gene        = factor(Gene, levels=rev(Gene)))
  if (nrow(df) < 2) return(NULL)
  p <- ggplot(df, aes(x=logFC, y=Gene)) +
    geom_segment(aes(x=0, xend=logFC, y=Gene, yend=Gene, color=Direction),
                 linewidth=0.9, alpha=0.8) +
    geom_point(aes(color=Direction, size=neg_log10p), alpha=0.9) +
    geom_vline(xintercept=0, color="#444444", linewidth=0.6, linetype="dashed") +
    scale_color_manual(values=c(Hypermethylated="#d73027",Hypomethylated="#4575b4")) +
    scale_size_continuous(name="-log₁₀(FDR)", range=c(2.5,7)) +
    labs(x="logFC vs Control", y="Gene") +
    theme_bw(base_size=11) +
    theme(axis.text.y=element_text(size=8, face="italic"), legend.position="right")
  ggsave(file.path(out_dir, paste0("Lollipop_", stage_name, ".png")),
         p, width=10, height=8, dpi=600)
  ggsave(file.path(out_dir, paste0("Lollipop_", stage_name, ".pdf")),
         p, width=10, height=8)
}

# --- Chromosomal distribution ---
plot_chromosomal_stage <- function(sig_ann, stage_name, out_dir) {
  df <- sig_ann[!is.na(sig_ann$chr) & !is.na(sig_ann$pos),]
  df$chr <- gsub("chr","",df$chr)
  df <- df[df$chr %in% as.character(1:22),]
  df$chr <- factor(df$chr, levels=rev(as.character(1:22)))
  df$Methylation_Status <- ifelse(df$logFC>0,"Hypermethylated","Hypomethylated")

  chr_lengths <- data.frame(
    chr    = factor(as.character(1:22), levels=rev(as.character(1:22))),
    length = c(248956422,242193529,198295559,190214555,181538259,
               170805979,159345973,145138636,138394717,133797422,
               135086622,133275309,114364328,107043718,101991189,
               90338345,83257441,80373285,58617616,64444167,
               46709983,50818468)
  )

  p <- ggplot() +
    geom_segment(data=chr_lengths, aes(x=0, xend=length, y=chr, yend=chr),
                 linewidth=3, color="grey65") +
    geom_point(data=df, aes(x=pos, y=chr, color=Methylation_Status),
               size=0.7, alpha=0.8) +
    scale_color_manual(values=c(Hypomethylated="blue",Hypermethylated="red")) +
    scale_x_continuous(labels=label_number(scale_cut=cut_short_scale())) +
    theme_minimal(base_size=11) +
    theme(panel.grid.major.y=element_blank(), legend.title=element_blank(),
          axis.text.y=element_text(face="bold")) +
    labs(x="Genomic Position", y="Chromosome")
  ggsave(file.path(out_dir, paste0("Chromosomal_", stage_name, ".png")),
         p, width=12, height=10, dpi=600)
  ggsave(file.path(out_dir, paste0("Chromosomal_", stage_name, ".pdf")),
         p, width=12, height=10)
}

# --- CpG region proportions ---
plot_cpg_proportions_stage <- function(sig_ann, stage_name, out_dir) {
  df <- sig_ann[!is.na(sig_ann$Relation_Primary),]
  df$Methylation_Status <- ifelse(df$logFC>0,"Hypermethylated","Hypomethylated")
  df$CpG_Region <- recode(df$Relation_Primary,
    "Island"="Island","N_Shore"="N Shore","S_Shore"="S Shore",
    "N_Shelf"="N Shelf","S_Shelf"="S Shelf","OpenSea"="Open Sea")
  df$CpG_Region <- factor(df$CpG_Region,
    levels=c("Island","N Shore","S Shore","N Shelf","S Shelf","Open Sea"))
  prop_df <- df %>%
    group_by(CpG_Region, Methylation_Status) %>%
    summarise(n=n(), .groups="drop") %>%
    group_by(CpG_Region) %>%
    mutate(Proportion=n/sum(n))
  if (nrow(prop_df) < 2) return(NULL)
  p <- ggplot(prop_df, aes(x=CpG_Region, y=Proportion, fill=Methylation_Status)) +
    geom_bar(stat="identity", width=0.7, color="black", linewidth=0.2) +
    scale_fill_manual(values=c(Hypomethylated="#4C72B0",Hypermethylated="#D62728")) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    theme_minimal(base_size=11) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
          legend.title=element_blank()) +
    labs(x=NULL, y="Proportion of DMPs")
  ggsave(file.path(out_dir, paste0("CpG_Proportions_", stage_name, ".png")),
         p, width=7, height=10, dpi=600)
  ggsave(file.path(out_dir, paste0("CpG_Proportions_", stage_name, ".pdf")),
         p, width=7, height=10)
}

# --- QQ / Lambda plot ---
plot_qq_stage <- function(res, stage_name, out_dir) {
  pvals <- res$P.Value

pvals <- pvals[
  !is.na(pvals) &
  is.finite(pvals) &
  pvals > 0 &
  pvals <= 1
]
  chisq  <- qchisq(1-pvals, df=1)
  lambda <- median(chisq)/qchisq(0.5,1)
  exp_p  <- -log10(ppoints(length(pvals)))
  obs_p  <- -log10(sort(pvals))
  qq_df  <- data.frame(Expected=exp_p, Observed=obs_p)
  p <- ggplot(qq_df, aes(Expected,Observed)) +
    geom_point(size=0.8, alpha=0.5, color="steelblue") +
    geom_abline(slope=1, intercept=0, color="red", linetype=2, linewidth=1) +
    theme_bw(base_size=13) +
    labs(subtitle=paste0("λ = ", round(lambda,3)),
         x="Expected -log10(P)", y="Observed -log10(P)")
  ggsave(file.path(out_dir, paste0("QQ_", stage_name, ".png")),
         p, width=6, height=6, dpi=600)
  return(lambda)
}

# --- GO/KEGG/Reactome ---
run_pathway_stage <- function(sig_ann, all_ann, stage_name, out_dir) {
  extract_genes <- function(gc) {
    g <- unique(trimws(unlist(strsplit(gc, ";"))))
    g[g != "" & !is.na(g)]
  }
  sig_genes <- extract_genes(
    paste(sig_ann$UCSC_RefGene_Name[!is.na(sig_ann$UCSC_RefGene_Name)],
          collapse=";"))
  bg_genes  <- extract_genes(
    paste(all_ann$UCSC_RefGene_Name[!is.na(all_ann$UCSC_RefGene_Name)],
          collapse=";"))

  bitr_safe <- function(genes) {
    tryCatch(bitr(genes, fromType="SYMBOL", toType="ENTREZID",
                  OrgDb=org.Hs.eg.db)$ENTREZID,
             error=function(e) character(0))
  }
  sig_e <- unique(bitr_safe(sig_genes))
  bg_e  <- unique(bitr_safe(bg_genes))
  if (length(sig_e) < 5) {
    cat("  Fewer than 5 genes for", stage_name, "- skipping pathway analysis\n")
    return(invisible(NULL))
  }

  run_enrich <- function(enrich_fn, ...) {
    tryCatch(enrich_fn(...), error=function(e) { cat("  Error:", e$message, "\n"); NULL })
  }

  # GO
  for (ont in c("BP","CC","MF")) {
    go_res <- run_enrich(enrichGO, gene=sig_e, universe=bg_e, OrgDb=org.Hs.eg.db,
                          keyType="ENTREZID", ont=ont, pAdjustMethod="BH",
                          qvalueCutoff=0.1, readable=TRUE)
    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      write.csv(as.data.frame(go_res),
                file.path(out_dir, paste0("GO_",ont,"_",stage_name,".csv")),
                row.names=FALSE)
      p <- dotplot(go_res, showCategory=20)
      ggsave(file.path(out_dir, paste0("GO_",ont,"_",stage_name,".png")),
             p, width=10, height=10, dpi=600)
      ggsave(file.path(out_dir, paste0("GO_",ont,"_",stage_name,".pdf")),
             p, width=10, height=10)
    }
  }

  # KEGG
  kegg_res <- run_enrich(enrichKEGG, gene=sig_e, universe=bg_e, organism="hsa",
                          keyType="ncbi-geneid", pAdjustMethod="BH",
                          pvalueCutoff=0.05, qvalueCutoff=0.1)
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    write.csv(as.data.frame(kegg_res),
              file.path(out_dir, paste0("KEGG_",stage_name,".csv")),
              row.names=FALSE)
    kegg_df <- as.data.frame(kegg_res) %>%
      arrange(p.adjust) %>% slice_head(n=20) %>%
      mutate(neglog10FDR=-log10(p.adjust),
             Description=factor(Description, levels=rev(Description)))
    p_k <- ggplot(kegg_df, aes(neglog10FDR, Description, size=Count, color=neglog10FDR)) +
      geom_point(alpha=0.9) +
      scale_color_gradient(low="#56B1F7",high="#D73027") +
      theme_bw(base_size=13) +
      labs(x=expression(-log[10]("FDR")), y=NULL)
    ggsave(file.path(out_dir, paste0("KEGG_",stage_name,".png")), p_k,
           width=11, height=10, dpi=600)
  }

  # Reactome
  react_res <- run_enrich(enrichPathway, gene=sig_e, universe=bg_e, organism="human",
                           pAdjustMethod="BH", qvalueCutoff=0.1, readable=TRUE)
  if (!is.null(react_res) && nrow(as.data.frame(react_res)) > 0) {
    write.csv(as.data.frame(react_res),
              file.path(out_dir, paste0("Reactome_",stage_name,".csv")),
              row.names=FALSE)
    p_r <- dotplot(react_res, showCategory=20)
    ggsave(file.path(out_dir, paste0("Reactome_",stage_name,".png")), p_r,
           width=10, height=10, dpi=600)
  }
}
