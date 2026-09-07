# ============================================================================
# PLOTS
# ============================================================================

EPLOT_DIR <- file.path(OUT_DIR, "Plots")
dir.create(EPLOT_DIR, showWarnings = FALSE)

# ---- SHARED THEME ---------------------------------------
theme_pub <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", size = base_size + 2,
                                        hjust = 0, colour = "#1a1a2e"),
      plot.subtitle     = element_text(size = base_size - 1, colour = "#555577",
                                        hjust = 0, margin = margin(b = 8)),
      plot.caption      = element_text(size = 9, colour = "#888888", hjust = 1),
      axis.title        = element_text(face = "bold", colour = "#1a1a2e"),
      axis.text         = element_text(colour = "#333333"),
      panel.grid.major  = element_line(colour = "#e8e8f0", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = "#ccccdd", fill = NA),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      strip.background  = element_rect(fill = "#f0f0f8", colour = "#ccccdd"),
      strip.text        = element_text(face = "bold", size = base_size - 1)
    )
}

# ============================================================================
# PLOT 1: TREND / ESCALATION LINE PLOT
# ============================================================================
cat("Generating Plot 1: Trend / Escalation line plot...\n")

ctrl_df <- summary_df %>%
  filter(GroupB == "Control") %>%
  mutate(
    Stage       = factor(GroupA, levels = STAGE_ORDER),
    Total       = Sig_FDR05_logFC02,
    Hyper       = Hyper,
    Hypo        = Hypo,
    Pct_Hyper   = ifelse(Total > 0, Hyper / Total * 100, 50)
  )
s3_count <- ctrl_df$Total[
  ctrl_df$Stage == "Stage3"
]

s4_count <- ctrl_df$Total[
  ctrl_df$Stage == "Stage4"
]

p1 <- ggplot(ctrl_df, aes(x = Stage)) +

  geom_area(aes(y = Total, group = 1),
            fill = "#d73027", alpha = 0.12) +

  geom_line(aes(y = Total, group = 1),
            colour = "#a50026", linewidth = 1.2, linetype = "solid") +

  geom_point(aes(y = Total, fill = Stage, size = Total),
             shape = 21, colour = "white", stroke = 1.5) +

  geom_text(aes(y = Total,
                label = ifelse(Total > 0,
                               formatC(Total, format = "d", big.mark = ","),
                               "0")),
            vjust = -1.3, fontface = "bold", size = 3.8,
            colour = "#1a1a2e") +

  annotate(
  "segment",
  x = 3.05,
  xend = 3.85,
  y = max(s3_count * 5, 1.5),
  yend = max(s4_count * 0.7, 2),
  arrow = arrow(
    length = unit(0.25, "cm"),
    type = "closed"
  ),
  colour = "#a50026",
  linewidth = 0.8
) +

  scale_y_log10(
    breaks  = c(1, 10, 100, 1000, 10000),
    limits  = c(0.5, 50000)
  ) +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  scale_size_continuous(range = c(4, 14), guide = "none") +

  labs(
    x        = "CKD Stage",
    y        = "Number of significant DMPs (log₁₀ scale)"
  ) +
  theme_pub()

ggsave(file.path(EPLOT_DIR, "Plot1_Trend_Escalation.png"),
       p1, width = 9, height = 6, dpi = 600)
cat("  Saved Plot1_Trend_Escalation.png\n")

# ============================================================================
# PLOT 2: LOLLIPOP PLOT — Top DMPs per comparison
# ============================================================================
cat("Generating Plot 2: Lollipop plots for top DMPs...\n")

make_lollipop <- function(df, stage_name, n_top = 10, color_up, color_dn) {

  sig <- df %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 0.2) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = n_top) %>%
    mutate(
      Direction = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated"),
      CpG_label = gsub("_TC21|_BC21|_TC11|_BC11", "", CpG),
      neg_log10p = -log10(adj.P.Val)
    ) %>%
    mutate(
      Gene = sapply(
        strsplit(
          as.character(UCSC_RefGene_Name),
          ";"
        ),
        `[`,
        1
      )
    ) %>%
    arrange(logFC) %>%
    mutate(CpG_label = factor(CpG_label, levels = CpG_label))

  if (nrow(sig) == 0) return(NULL)

  ggplot(sig, aes(x = logFC, y = CpG_label)) +
    geom_segment(aes(x = 0, xend = logFC,
                     y = CpG_label, yend = CpG_label,
                     colour = Direction),
                 linewidth = 0.9, alpha = 0.8) +
    geom_point(aes(colour = Direction, size = neg_log10p),
               alpha = 0.9) +
    geom_vline(xintercept = 0, colour = "#444444",
               linewidth = 0.6, linetype = "dashed") +
    scale_colour_manual(
      values  = c("Hypermethylated" = color_up,
                  "Hypomethylated"  = color_dn),
      name    = "Direction"
    ) +
    scale_size_continuous(
      name   = "-log₁₀(FDR)",
      range  = c(2.5, 7),
      breaks = c(2, 5, 10, 20)
    ) +
    labs(
      x        = "logFC vs Control (positive = hypermethylated in stage)",
      y        = "CpG site"
    ) +
    theme_pub(base_size = 11) +
    theme(
      axis.text.y      = element_text(size = 8, family = "mono"),
      legend.position  = "right",
      legend.direction = "vertical"
    )
}

stage4_ann <- merge(
  dmps_list[["Stage4_vs_Control"]],
  ann_sub,
  by.x = "CpG",
  by.y = "Name",
  all.x = TRUE
)

stage5_ann <- merge(
  dmps_list[["Stage5_vs_Control"]],
  ann_sub,
  by.x = "CpG",
  by.y = "Name",
  all.x = TRUE
)

p2a <- make_lollipop(stage4_ann,
                     "Stage4", n_top = 20,
                     color_up = "#d73027", color_dn = "#4575b4")
p2b <- make_lollipop(stage5_ann,
                     "Stage5", n_top = 20,
                     color_up = "#a50026", color_dn = "#313695")

if (!is.null(p2a)) {
  ggsave(file.path(EPLOT_DIR, "Plot2a_Lollipop_Stage4.png"),
         p2a, width = 10, height = 8, dpi = 600)
  cat("  Saved Plot2a_Lollipop_Stage4.png\n")
}
if (!is.null(p2b)) {
  ggsave(file.path(EPLOT_DIR, "Plot2b_Lollipop_Stage5.png"),
         p2b, width = 10, height = 8, dpi = 600)
  cat("  Saved Plot2b_Lollipop_Stage5.png\n")
}

# ============================================================================
# PLOT 3: STACKED PROPORTION BAR — Hyper/Hypo balance
# ============================================================================
cat("Generating Plot 3: Stacked proportion bar...\n")

prop_df <- summary_df %>%
  filter(Sig_FDR05_logFC02 > 0) %>%
  mutate(
    Pct_Hyper = Hyper / Sig_FDR05_logFC02 * 100,
    Pct_Hypo  = Hypo  / Sig_FDR05_logFC02 * 100,
    Comparison = factor(Comparison, levels = Comparison[order(GroupA, GroupB)])
  ) %>%
  pivot_longer(cols      = c(Pct_Hyper, Pct_Hypo),
               names_to  = "Direction",
               values_to = "Proportion") %>%
  mutate(
    Direction = recode(Direction,
                       "Pct_Hyper" = "Hypermethylated in numerator",
                       "Pct_Hypo"  = "Hypomethylated in numerator"),
    IsCtrl = grepl("Control", Comparison)
  )

p3 <- ggplot(prop_df, aes(x = Comparison, y = Proportion, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed",
             colour = "white", linewidth = 0.8) +
  geom_text(
    data = summary_df %>%
      filter(Sig_FDR05_logFC02 > 0) %>%
      mutate(
        Comparison = factor(Comparison, levels = levels(prop_df$Comparison)),
        IsCtrl     = grepl("Control", Comparison)
      ),
    aes(x = Comparison, y = 103,
        label = formatC(Sig_FDR05, format = "d", big.mark = ",")),
    inherit.aes = FALSE,
    size = 2.8, colour = "#333333", fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Hypermethylated in numerator" = "#d73027",
               "Hypomethylated in numerator"  = "#4575b4"),
    name   = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 107),
    breaks = c(0, 25, 50, 75, 100),
    labels = paste(c(0, 25, 50, 75, 100), "%")
  ) +
  facet_wrap(~ IsCtrl,
             scales   = "free_x",
             labeller = labeller(IsCtrl = c(
               "TRUE"  = "Stage vs Control",
               "FALSE" = "Stage vs Stage"
             ))) +
  labs(
    x        = NULL,
    y        = "Proportion of significant DMPs (%)"
  ) +
  theme_pub() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(EPLOT_DIR, "Plot3_Proportion_HyperHypo.png"),
       p3, width = 14, height = 7, dpi = 600)
cat("  Saved Plot3_Proportion_HyperHypo.png\n")

# ============================================================================
# PLOT 4: PER-STAGE logFC STRIP / HEATMAP
# ============================================================================
cat("Generating Plot 4: Per-stage logFC strip heatmap...\n")

if (nrow(globally_sig) > 0) {

  lfc_cols <- c("Stage1_vs_Control", "Stage2_vs_Control",
                 "Stage3_vs_Control", "Stage4_vs_Control",
                 "Stage5_vs_Control")

  n_top    <- min(60, nrow(globally_sig))
  top_cpgs <- rownames(
  head(
    globally_sig,
    n_top
  )
)

lfc_mat <- data.frame(
  CpG = top_cpgs
)

for (cc in lfc_cols) {

  tmp <- dmps_list[[cc]][,
    c("CpG", "logFC")
  ]

  colnames(tmp)[2] <- cc

  lfc_mat <- merge(
    lfc_mat,
    tmp,
    by = "CpG",
    all.x = TRUE
  )
}

  strip_df <- lfc_mat %>%
    pivot_longer(cols      = all_of(lfc_cols),
                 names_to  = "Stage_Comparison",
                 values_to = "logFC") %>%
    mutate(
      Stage = gsub("_vs_Control", "", Stage_Comparison),
      Stage = factor(Stage, levels = STAGE_ORDER),
      CpG_clean = gsub("_TC21|_BC21|_TC11|_BC11", "", CpG)
    )

  cpg_order <- lfc_mat %>%
    arrange(desc(!!sym("Stage5_vs_Control"))) %>%
    mutate(CpG_clean = gsub("_TC21|_BC21|_TC11|_BC11", "", CpG)) %>%
    pull(CpG_clean)

  strip_df$CpG_clean <- factor(strip_df$CpG_clean, levels = rev(cpg_order))

  lfc_cap <- quantile(abs(strip_df$logFC), 0.99, na.rm = TRUE)
  strip_df$logFC_capped <- pmin(pmax(strip_df$logFC, -lfc_cap), lfc_cap)

  p5 <- ggplot(strip_df,
               aes(x = Stage, y = CpG_clean, fill = logFC_capped)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low      = "#313695",
      mid      = "#fffff0",
      high     = "#a50026",
      midpoint = 0,
      name     = "logFC vs\nControl",
      limits   = c(-lfc_cap, lfc_cap)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x        = NULL,
      y        = "CpG site"
    ) +
    theme_pub(base_size = 10) +
    theme(
      axis.text.y     = element_text(size = 6.5, family = "mono"),
      axis.text.x     = element_text(face = "bold", size = 11),
      legend.position = "right",
      panel.border    = element_rect(colour = "#aaaaaa")
    )

  ggsave(file.path(EPLOT_DIR, "Plot4_Stage_logFC_Strip.png"),
         p5, width = 8, height = 12, dpi = 600)
  cat("  Saved Plot4_Stage_logFC_Strip.png\n")

} else {
  cat("  Skipped Plot 4 — no globally significant CpGs\n")
}

# ============================================================================
# PLOT 5: MANHATTAN PLOT (ggplot2)
# ============================================================================
cat("Generating Plot 5: Manhattan plot (Stage5 vs Control)...\n")

make_manhattan <- function(result_df, stage_name,
                            fdr_threshold   = 0.05,
                            anno_object     = ann,
                            top_n_label     = 10) {

  anno_pos <- as.data.frame(anno_object[, c("Name", "chr", "pos")])
  df_ann   <- merge(result_df, anno_pos,
                    by.x = "CpG", by.y = "Name", all.x = TRUE)
  df_ann   <- df_ann[!is.na(df_ann$chr) & !is.na(df_ann$pos), ]
  df_ann   <- df_ann[!df_ann$chr %in% c("chrX", "chrY", "chrM"), ]

  chr_nums <- as.numeric(gsub("chr", "", df_ann$chr))
  df_ann   <- df_ann[!is.na(chr_nums),  ]
  chr_nums <- as.numeric(gsub("chr", "", df_ann$chr))
  df_ann   <- df_ann[order(chr_nums, df_ann$pos), ]
  df_ann$chr_num <- as.numeric(gsub("chr", "", df_ann$chr))

  chr_offsets <- df_ann %>%
    group_by(chr_num) %>%
    summarise(max_pos = max(pos), .groups = "drop") %>%
    arrange(chr_num) %>%
    mutate(offset = cumsum(lag(max_pos, default = 0)) +
                    (chr_num - 1) * 5e7)

  df_ann <- df_ann %>%
    left_join(chr_offsets, by = "chr_num") %>%
    mutate(
      cum_pos    = pos + offset,
      neg_log10p = pmin(-log10(adj.P.Val), 20),   # cap at 20 for plotting
      Sig        = adj.P.Val < fdr_threshold & abs(logFC) > 0.2,
      chr_colour = ifelse(chr_num %% 2 == 0, "even", "odd")
    )

  axis_df <- df_ann %>%
    group_by(chr_num) %>%
    summarise(centre = mean(cum_pos), .groups = "drop") %>%
    mutate(chr_label = as.character(chr_num))

  top_label <- df_ann %>%
    filter(Sig) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = top_n_label) %>%
    mutate(CpG_label = gsub("_TC21|_BC21|_TC11|_BC11", "", CpG))

  sig_line <- -log10(fdr_threshold)

  ggplot(df_ann, aes(x = cum_pos, y = neg_log10p)) +

    geom_point(data = df_ann[!df_ann$Sig, ],
               aes(colour = chr_colour),
               size = 0.5, alpha = 0.5) +

    geom_point(data = df_ann[df_ann$Sig, ],
               colour = "#a50026", size = 1.2, alpha = 0.8) +

    geom_hline(yintercept = sig_line,
               linetype = "dashed", colour = "#333333",
               linewidth = 0.7) +

    ggrepel::geom_label_repel(
      data          = top_label,
      aes(label     = CpG_label),
      size          = 2.5,
      colour        = "#a50026",
      fill          = "white",
      box.padding   = 0.4,
      label.padding = 0.2,
      max.overlaps  = 15,
      segment.colour = "#d73027",
      segment.size   = 0.4
    ) +

    annotate("text",
             x = max(df_ann$cum_pos) * 0.02,
             y = sig_line + 0.3,
             label = paste0("FDR = ", fdr_threshold),
             hjust = 0, size = 3, colour = "#333333") +

    scale_colour_manual(
      values = c(even = "#74add1", odd = "#2166ac"),
      guide  = "none"
    ) +
    scale_x_continuous(
      breaks = axis_df$centre,
      labels = axis_df$chr_label,
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 2, 5, 10, 15, 20),
      labels = c("0", "2", "5", "10", "15", "≥20"),
      expand = c(0.02, 0)
    ) +
    labs(
      x        = "Chromosome",
      y        = expression(-log[10](FDR))
    ) +
    theme_pub() +
    theme(
      axis.text.x     = element_text(size = 9),
      panel.grid.major.x = element_blank()
    )
}

p6a <- make_manhattan(dmps_list[["Stage5_vs_Control"]],
                       "Stage5", fdr_threshold = 0.05)
ggsave(file.path(EPLOT_DIR, "Plot5a_Manhattan_Stage5.png"),
       p6a, width = 14, height = 6, dpi = 600)
cat("  Saved Plot5a_Manhattan_Stage5.png\n")

p6b <- make_manhattan(dmps_list[["Stage4_vs_Control"]],
                       "Stage4", fdr_threshold = 0.05)
ggsave(file.path(EPLOT_DIR, "Plot5b_Manhattan_Stage4.png"),
       p6b, width = 14, height = 6, dpi = 600)
cat("  Saved Plot5b_Manhattan_Stage4.png\n")

# ============================================================================
# PLOT 6: PAIRWISE DMP COUNT MATRIX
# ============================================================================

cat("Plot 6: Pairwise DMP count matrix...\n")

mat_df <- expand.grid(GroupA = all_groups, GroupB = all_groups,
                      stringsAsFactors = FALSE)
mat_df <- mat_df[mat_df$GroupA != mat_df$GroupB, ]

mat_df <- merge(mat_df, summary_df[, c("GroupA","GroupB","Sig_FDR05_logFC02")],
                by = c("GroupA","GroupB"), all.x = TRUE)
mat_df$Sig_FDR05_logFC02[is.na(mat_df$Sig_FDR05_logFC02)] <- 0

mat_df$GroupA <- factor(mat_df$GroupA, levels = all_groups)
mat_df$GroupB <- factor(mat_df$GroupB, levels = rev(all_groups))

p1 <- ggplot(mat_df, aes(x = GroupA, y = GroupB, fill = Sig_FDR05_logFC02)) +
  geom_tile(color = "black", linewidth = 0.8) +
  geom_text(aes(label = ifelse(Sig_FDR05_logFC02 > 0, Sig_FDR05_logFC02, "")),
            size = 3.5, fontface = "bold", color = "black") +
  scale_fill_gradientn(
    colors  = c("#f7f7f7", "#d1e5f0", "#4393c3", "#2166ac", "#053061"),
    values  = rescale(c(0, 1, 50, 200, max(mat_df$Sig_FDR05_logFC02, 1))),
    name    = "Significant\nDMPs\n(FDR<0.05)",
    na.value = "grey90"
  ) +
  labs(
    x        = "Numerator Group (A in A vs B)",
    y        = "Denominator Group (B in A vs B)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x      = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 15),
    legend.position  = "right"
  )

ggsave(file.path(EPLOT_DIR, "Plot6_DMP_Pairwise_Heatmap.png"),
       p1, width = 8, height = 7, dpi = 600)
cat("  -> Saved Plot 6\n")
