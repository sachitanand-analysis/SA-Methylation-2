# ============================================================
# CHROD PLOT
# ============================================================
library(circlize)

S4_FILE <- file.path(
  DMP_DIR,
  "Sig_DMPs_annotated_Stage4_vs_Control.csv"
)

S5_FILE <- file.path(
  DMP_DIR,
  "Sig_DMPs_annotated_Stage5_vs_Control.csv"
)

PLOT_OUT <- file.path(OUT_DIR, "ChordPlot_Stage4_Stage5_SharedCpGs.png")

# READ DMP FILES

s4 <- read.csv(S4_FILE, stringsAsFactors = FALSE)
s5 <- read.csv(S5_FILE, stringsAsFactors = FALSE)

cat("Stage4 DMPs loaded:", nrow(s4), "\n")
cat("Stage5 DMPs loaded:", nrow(s5), "\n")

# FIX CpG COLUMN

fix_cpg_column <- function(df) {
  if (!"CpG" %in% colnames(df)) {
    cpg_col <- intersect(
      c("Name", "ProbeID", "Probe", "probe", "cpg", "ID"),
      colnames(df)
    )
    if (length(cpg_col) > 0) {
      df$CpG <- df[[cpg_col[1]]]
    } else {
      df$CpG <- rownames(df)
    }
  }
  return(df)
}

s4 <- fix_cpg_column(s4)
s5 <- fix_cpg_column(s5)


s4 <- s4 %>% mutate(Direction_S4 = ifelse(logFC > 0, "Hyper", "Hypo"))
s5 <- s5 %>% mutate(Direction_S5 = ifelse(logFC > 0, "Hyper", "Hypo"))

# FIND SHARED AND UNIQUE

shared_cpgs <- intersect(s4$CpG, s5$CpG)
only_s4     <- setdiff(s4$CpG, s5$CpG)
only_s5     <- setdiff(s5$CpG, s4$CpG)

cat("\nShared CpGs:", length(shared_cpgs), "\n")
cat("Stage4 only:", length(only_s4), "\n")
cat("Stage5 only:", length(only_s5), "\n")

# CLASSIFY SHARED CpGs

shared_df <- s4 %>%
  filter(CpG %in% shared_cpgs) %>%
  select(CpG, Direction_S4) %>%
  left_join(s5 %>% select(CpG, Direction_S5), by = "CpG") %>%
  mutate(
    SharedCategory = case_when(
      Direction_S4 == "Hyper" & Direction_S5 == "Hyper" ~ "Shared_BothHyper",
      Direction_S4 == "Hypo"  & Direction_S5 == "Hypo"  ~ "Shared_BothHypo",
      Direction_S4 == "Hypo"  & Direction_S5 == "Hyper" ~ "Shared_HypoS4_HyperS5",
      Direction_S4 == "Hyper" & Direction_S5 == "Hypo"  ~ "Shared_HyperS4_HypoS5",
      TRUE ~ "Unknown"
    )
  )

cat("\nShared direction breakdown:\n")
print(table(shared_df$SharedCategory))

# COUNTS PER CATEGORY

s4_only_hyper    <- sum(s4$CpG %in% only_s4 & s4$Direction_S4 == "Hyper")
s4_only_hypo     <- sum(s4$CpG %in% only_s4 & s4$Direction_S4 == "Hypo")
s5_only_hyper    <- sum(s5$CpG %in% only_s5 & s5$Direction_S5 == "Hyper")
s5_only_hypo     <- sum(s5$CpG %in% only_s5 & s5$Direction_S5 == "Hypo")
shared_both_hyper    <- sum(shared_df$SharedCategory == "Shared_BothHyper")
shared_both_hypo     <- sum(shared_df$SharedCategory == "Shared_BothHypo")
shared_hypo_to_hyper <- sum(shared_df$SharedCategory == "Shared_HypoS4_HyperS5")
shared_hyper_to_hypo <- sum(shared_df$SharedCategory == "Shared_HyperS4_HypoS5")

cat("\n--- Category counts ---\n")
cat("S4 only Hyper:        ", s4_only_hyper, "\n")
cat("S4 only Hypo:         ", s4_only_hypo,  "\n")
cat("S5 only Hyper:        ", s5_only_hyper, "\n")
cat("S5 only Hypo:         ", s5_only_hypo,  "\n")
cat("Shared Both Hyper:    ", shared_both_hyper,    "\n")
cat("Shared Both Hypo:     ", shared_both_hypo,     "\n")
cat("Hypo(S4)->Hyper(S5):  ", shared_hypo_to_hyper, "\n")
cat("Hyper(S4)->Hypo(S5):  ", shared_hyper_to_hypo, "\n")

# NODE NAMES WITH DMP COUNTS EMBEDDED

nodes <- c(
  paste0("S4 Hyper only\n(n=", s4_only_hyper, ")"),
  paste0("S4 Hypo only\n(n=", s4_only_hypo, ")"),
  paste0("S5 Hyper only\n(n=", s5_only_hyper, ")"),
  paste0("S5 Hypo only\n(n=", s5_only_hypo, ")"),
  paste0("Shared Both Hyper\n(n=", shared_both_hyper, ")"),
  paste0("Shared Both Hypo\n(n=", shared_both_hypo, ")"),
  paste0("Hypo\u2192Hyper (S4\u2192S5)\n(n=", shared_hypo_to_hyper, ")"),
  paste0("Hyper\u2192Hypo (S4\u2192S5)\n(n=", shared_hyper_to_hypo, ")")
)

n   <- length(nodes)
mat <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))

mat[nodes[1], nodes[1]] <- s4_only_hyper
mat[nodes[2], nodes[2]] <- s4_only_hypo
mat[nodes[3], nodes[3]] <- s5_only_hyper
mat[nodes[4], nodes[4]] <- s5_only_hypo

mat[nodes[5], nodes[1]] <- shared_both_hyper    # Shared Both Hyper  <-> S4 Hyper
mat[nodes[5], nodes[3]] <- shared_both_hyper    # Shared Both Hyper  <-> S5 Hyper

mat[nodes[6], nodes[2]] <- shared_both_hypo     # Shared Both Hypo   <-> S4 Hypo
mat[nodes[6], nodes[4]] <- shared_both_hypo     # Shared Both Hypo   <-> S5 Hypo

mat[nodes[7], nodes[2]] <- shared_hypo_to_hyper # Hypo->Hyper switch <-> S4 Hypo
mat[nodes[7], nodes[3]] <- shared_hypo_to_hyper # Hypo->Hyper switch <-> S5 Hyper

mat[nodes[8], nodes[1]] <- shared_hyper_to_hypo # Hyper->Hypo switch <-> S4 Hyper
mat[nodes[8], nodes[4]] <- shared_hyper_to_hypo # Hyper->Hypo switch <-> S5 Hypo

# ============================================================
# COLOURS
#   S4 Hyper  : vivid orange
#   S4 Hypo   : pistachio green
#   S5 Hyper  : deep purple
#   S5 Hypo   : turquoise
#   Shared BH : gold
#   Shared ByH: teal
#   Switch H>H: bright red
#   Switch H>y: coral red
# ============================================================

sector_colors <- c(
  "#E07B00",   # S4 Hyper only      — vivid orange
  "#6DB33F",   # S4 Hypo only       — pistachio green
  "#6A0DAD",   # S5 Hyper only      — deep purple
  "#00B4D8",   # S5 Hypo only       — turquoise blue
  "#FFD166",   # Shared Both Hyper  — gold/amber
  "#06D6A0",   # Shared Both Hypo   — teal
  "#FF0000",   # Hypo->Hyper switch — bright red
  "#FF6B35"    # Hyper->Hypo switch — coral orange
)
names(sector_colors) <- nodes

# PLOT

png(PLOT_OUT, width = 7200, height = 7200, res = 600)

par(
  mar = c(1, 1, 1, 1),   
  oma = c(2, 2, 12, 2)    
)

circos.clear()

circos.par(
  "track.height"           = 0.10,
  "gap.after"              = c(rep(2, n - 1), 6),
  "start.degree"           = 100,      
  "clock.wise"             = TRUE,
  canvas.ylim = c(-1.2, 1.4)
)

chordDiagram(
  mat,
  grid.col          = sector_colors,
  transparency      = 0.30,
  annotationTrack   = "grid",
  preAllocateTracks = list(track.height = mm_h(10)),
  self.link         = 1,
  scale             = FALSE,
  directional       = 0,
  reduce            = 0,
  link.lwd          = 2,
  link.border       = "black"
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun   = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    xlim        <- get.cell.meta.data("xlim")
    ylim        <- get.cell.meta.data("ylim")
    circos.text(
      mean(xlim),
      ylim[1] + 0.45,       # push labels outward
      sector_name,
      facing     = "clockwise",
      niceFacing = TRUE,
      adj        = c(0, 0.5),
      cex        = 0.68,
      font       = 2
    )
  },
  bg.border = NA
)

# Legend
legend(
  x      = "bottomleft",
  legend = c(
    paste0("Stage 4 only — Hypermethylated (n=", s4_only_hyper, ")"),
    paste0("Stage 4 only — Hypomethylated  (n=", s4_only_hypo,  ")"),
    paste0("Stage 5 only — Hypermethylated (n=", s5_only_hyper, ")"),
    paste0("Stage 5 only — Hypomethylated  (n=", s5_only_hypo,  ")"),
    paste0("Shared — Both Hypermethylated  (n=", shared_both_hyper,    ")"),
    paste0("Shared — Both Hypomethylated   (n=", shared_both_hypo,     ")"),
    paste0("Switch: Hypo(S4)\u2192Hyper(S5)     (n=", shared_hypo_to_hyper, ")"),
    paste0("Switch: Hyper(S4)\u2192Hypo(S5)     (n=", shared_hyper_to_hypo, ")")
  ),
  fill   = unname(sector_colors),
  border = NA,
  bty    = "n",
  cex    = 0.82,
  title  = "Category",
  title.font = 2
)

circos.clear()
dev.off()

cat("\nPlot saved to:", PLOT_OUT, "\n")

summary_out <- data.frame(
  Category = c(
    "Stage4 only — Hyper",
    "Stage4 only — Hypo",
    "Stage5 only — Hyper",
    "Stage5 only — Hypo",
    "Shared — Both Hyper",
    "Shared — Both Hypo",
    "Switch: Hypo(S4) to Hyper(S5)",
    "Switch: Hyper(S4) to Hypo(S5)"
  ),
  Count = c(
    s4_only_hyper, s4_only_hypo,
    s5_only_hyper, s5_only_hypo,
    shared_both_hyper, shared_both_hypo,
    shared_hypo_to_hyper, shared_hyper_to_hypo
  )
)

write.csv(
  summary_out,
  file.path(OUT_DIR, "ChordPlot_Category_Summary_v2.csv"),
  row.names = FALSE
)

write.csv(
  shared_df,
  file.path(OUT_DIR, "Shared_CpGs_Stage4_Stage5_Classified.csv"),
  row.names = FALSE
)

cat("Summary and shared CpG tables saved.\n")
