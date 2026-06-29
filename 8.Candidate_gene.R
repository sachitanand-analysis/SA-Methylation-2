# ============================================================
# 1. INPUTS
# ============================================================

# Folder containing ALL significant DMP CSV files
DMP_DIR <- ""

# Candidate gene list
GENE_FILE <- ""

# Output folder
OUT_DIR <- file.path(DMP_DIR, "Candidate_Gene_Overlap")
dir.create(OUT_DIR, showWarnings = FALSE)

# READ GENE LIST

gene_df <- read.csv(GENE_FILE, stringsAsFactors = FALSE)

candidate_genes <- unique(
  toupper(trimws(gene_df$Gene))
)

cat("Total candidate genes:",
    length(candidate_genes), "\n")

# GET ALL DMP FILES

dmp_files <- list.files(
  DMP_DIR,
  pattern = "Sig_DMPs_annotated_.*\\.csv$",
  full.names = TRUE
)

cat("Total DMP files found:",
    length(dmp_files), "\n")

# EMPTY SUMMARY OBJECT

summary_list <- list()

# LOOP THROUGH EACH DMP FILE

for (f in dmp_files) {

  cat("\n====================================\n")
  cat("Processing:", basename(f), "\n")
  cat("====================================\n")

  dmp <- read.csv(f, stringsAsFactors = FALSE)
  
  required_cols <- c("UCSC_RefGene_Name", "logFC")

  missing_cols <- setdiff(required_cols, colnames(dmp))

  if (length(missing_cols) > 0) {

    cat("Skipping file - missing columns:",
        paste(missing_cols, collapse = ", "),
        "\n")

    next
  }
  
  dmp_expanded <- dmp %>%
    mutate(
      UCSC_RefGene_Name =
        as.character(UCSC_RefGene_Name)
    ) %>%
    separate_rows(
      UCSC_RefGene_Name,
      sep = ";"
    ) %>%
    mutate(
      Gene = toupper(trimws(UCSC_RefGene_Name))
    ) %>%
    filter(
      Gene != "",
      !is.na(Gene)
    )
  
  overlap <- dmp_expanded %>%
    filter(Gene %in% candidate_genes)
  
  overlap <- overlap %>%
    mutate(
      Methylation =
        ifelse(logFC > 0,
               "Hypermethylated",
               "Hypomethylated")
    )
  
  overlap_unique <- overlap %>%
    distinct(
      Gene,
      CpG,
      logFC,
      adj.P.Val,
      Methylation,
      .keep_all = TRUE
    )
 
  comparison_name <- gsub(
    "Sig_DMPs_annotated_|\\.csv",
    "",
    basename(f)
  )

  out_file <- file.path(
    OUT_DIR,
    paste0(
      "CandidateGeneOverlap_",
      comparison_name,
      ".csv"
    )
  )

  write.csv(
    overlap_unique,
    out_file,
    row.names = FALSE
  )

  hyper_n <- sum(
    overlap_unique$Methylation ==
      "Hypermethylated"
  )

  hypo_n <- sum(
    overlap_unique$Methylation ==
      "Hypomethylated"
  )

  total_n <- nrow(overlap_unique)

  unique_genes <- length(
    unique(overlap_unique$Gene)
  )

  summary_list[[comparison_name]] <-
    data.frame(
      Comparison = comparison_name,
      Total_DMPs = total_n,
      Unique_Genes = unique_genes,
      Hyper = hyper_n,
      Hypo = hypo_n
    )

  cat("Candidate overlaps found:",
      total_n, "\n")
}

summary_df <- bind_rows(summary_list)

summary_df <- summary_df %>%
  arrange(desc(Total_DMPs))

write.csv(
  summary_df,
  file.path(
    OUT_DIR,
    "CandidateGeneOverlap_Summary.csv"
  ),
  row.names = FALSE
)

cat("\n====================================\n")
cat("FINAL SUMMARY\n")
cat("====================================\n")

print(summary_df)

# ============================================================
# 2. CREATE GENE x COMPARISON MATRIX
# ============================================================

all_overlap_files <- list.files(
  OUT_DIR,
  pattern = "CandidateGeneOverlap_.*\\.csv$",
  full.names = TRUE
)

gene_matrix_list <- list()

for (f in all_overlap_files) {

  df <- read.csv(f, stringsAsFactors = FALSE)

  comp <- gsub(
    "CandidateGeneOverlap_|\\.csv",
    "",
    basename(f)
  )

 if (nrow(df) == 0) next

required_cols <- c("Gene", "Methylation")

if (!all(required_cols %in% colnames(df))) {

  cat("Skipping malformed file:",
      basename(f), "\n")

  next
}

tmp <- df %>%
  select(Gene, Methylation) %>%
  distinct()

  tmp[[comp]] <- tmp$Methylation

  tmp <- tmp %>%
    select(Gene, all_of(comp))

  gene_matrix_list[[comp]] <- tmp
}

if (length(gene_matrix_list) > 0) {

  gene_matrix <- Reduce(
    function(x, y)
      full_join(x, y, by = "Gene"),
    gene_matrix_list
  )

  write.csv(
    gene_matrix,
    file.path(
      OUT_DIR,
      "CandidateGene_ComparisonMatrix.csv"
    ),
    row.names = FALSE
  )
}

# ============================================================
# 3. CREATE CpG x COMPARISON MATRIX
# ============================================================

cpg_matrix_list <- list()

for (f in all_overlap_files) {

  df <- read.csv(f, stringsAsFactors = FALSE)

  comp <- gsub(
    "CandidateGeneOverlap_|\\.csv",
    "",
    basename(f)
  )
  
  if (nrow(df) == 0) next

  required_cols <- c(
    "CpG",
    "Gene",
    "Methylation"
  )

  if (!all(required_cols %in% colnames(df))) {

    cat("Skipping malformed file:",
        basename(f), "\n")

    next
  }

  tmp <- df %>%
    select(CpG, Gene, Methylation) %>%
    distinct()

  tmp[[comp]] <- tmp$Methylation

  tmp <- tmp %>%
    select(CpG, Gene, all_of(comp))

  cpg_matrix_list[[comp]] <- tmp
}

if (length(cpg_matrix_list) > 0) {

  cpg_matrix <- Reduce(
    function(x, y)
      full_join(x, y,
                by = c("CpG", "Gene")),
    cpg_matrix_list
  )

  write.csv(
    cpg_matrix,
    file.path(
      OUT_DIR,
      "CandidateGene_CpG_ComparisonMatrix.csv"
    ),
    row.names = FALSE
  )
}

cat("\nCpG comparison matrix saved.\n")
cat("\nDONE.\n")

# ============================================================
# 4. SANKEY FLOW PLOT
# ============================================================

summary_df <- data.frame(

Comparison = c(
"Stage4",
"Stage5",
"Overall CKD"
),

Total_DMPs = c(10,33,176),

Hyper = c(4,13,80),

Hypo = c(6,20,96)
)

comparison_nodes <- summary_df$Comparison

total_nodes <- paste0(
summary_df$Comparison,
"\nTotal DMPs (n=",
summary_df$Total_DMPs,
")"
)

hyper_nodes <- paste0(
summary_df$Comparison,
"\nHyper (n=",
summary_df$Hyper,
")"
)

hypo_nodes <- paste0(
summary_df$Comparison,
"\nHypo (n=",
summary_df$Hypo,
")"
)

nodes <- data.frame(
name = c(
comparison_nodes,
total_nodes,
hyper_nodes,
hypo_nodes
)
)

links <- data.frame()

for(i in 1:nrow(summary_df)) {

comp  <- comparison_nodes[i]
total <- total_nodes[i]
hyper <- hyper_nodes[i]
hypo  <- hypo_nodes[i]

links <- rbind(
links,
data.frame(
source = match(comp, nodes$name) - 1,
target = match(total, nodes$name) - 1,
value  = summary_df$Total_DMPs[i],
group  = "Total"
)
)

links <- rbind(
links,
data.frame(
source = match(total, nodes$name) - 1,
target = match(hyper, nodes$name) - 1,
value  = summary_df$Hyper[i],
group  = "Hyper"
)
)

links <- rbind(
links,
data.frame(
source = match(total, nodes$name) - 1,
target = match(hypo, nodes$name) - 1,
value  = summary_df$Hypo[i],
group  = "Hypo"
)
)
}

my_color <- '
d3.scaleOrdinal()

.domain([

"Stage4",
"Stage5",
"Overall CKD",

"Total",

"Hyper",
"Hypo"

])

.range([

"#1B9E77",   /* Stage4 green */
"#D95F02",   /* Stage5 orange */
"#7570B3",   /* Overall CKD purple */

"#BDBDBD",   /* Total grey */

"#B2182B",   /* Hyper red */
"#2166AC"    /* Hypo blue */

])
'

nodes$group <- c(

rep("Comparison", 3),

rep("Total", 3),

rep("Hyper", 3),

rep("Hypo", 3)

)

links$group <- c(

"Stage4",
"Stage5",
"Overall CKD",

"Hyper",
"Hyper",
"Hyper",

"Hypo",
"Hypo",
"Hypo"

)

p_sankey <- sankeyNetwork(

Links = links,
Nodes = nodes,

Source = "source",
Target = "target",
Value  = "value",

NodeID = "name",

NodeGroup = "group",
LinkGroup = "group",

sinksRight = FALSE,

fontSize = 15,

nodeWidth = 18,

nodePadding = 30,

colourScale = my_color

)

htmlwidgets::saveWidget(

p_sankey,

file.path(
OUT_DIR,
"CKD_Sankey.html"
),

selfcontained = TRUE
)

cat("Saved improved Sankey plot.\n")

# ============================================================
# 5. CANDIDATE GENE HEATMAP
# ============================================================

beta_vals <- getBeta(mSet_final)

mapping <- read.csv("")

STAGE_COLORS <- c(
"Control" = "#1b9e77",
"Stage1"  = "#d95f02",
"Stage2"  = "#7570b3",
"Stage3"  = "#e7298a",
"Stage4"  = "#66a61e",
"Stage5"  = "#e6ab02"
)

sig_cpgs <- unique(cpg_matrix$CpG)

heat_beta <- beta_vals[
rownames(beta_vals) %in% sig_cpgs,
]

heat_beta <- heat_beta[
complete.cases(heat_beta),
]

heat_scaled <- t(scale(t(heat_beta)))

heat_scaled[heat_scaled > 3]  <- 3
heat_scaled[heat_scaled < -3] <- -3

ann_col <- data.frame(
Group = mapping$Group
)

rownames(ann_col) <- colnames(heat_scaled)

grp_colors <- list(
Group = STAGE_COLORS
)

row_labels <- paste(
  cpg_matrix$Gene[
    match(
      rownames(heat_scaled),
      cpg_matrix$CpG
    )
  ],
  rownames(heat_scaled),
  sep = " | "
)

png(
file.path(
OUT_DIR,
"CandidateGene_Heatmap.png"
),
width = 4800,
height = 3600,
res = 600
)

pheatmap(
heat_scaled,

annotation_col = ann_col,
annotation_colors = grp_colors,

clustering_method = "ward.D2",

color = colorRampPalette(
c("#313695", "white", "#a50026")
)(100),

show_colnames = FALSE,
show_rownames = TRUE,

labels_row = row_labels,

fontsize_row = 5
)

dev.off()

cat("Candidate gene heatmap saved.\n")
