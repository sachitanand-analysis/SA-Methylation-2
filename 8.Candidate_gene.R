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
# 4. SANKEY FLOW PLOT - overall CKD
# ============================================================

overall <- read.csv(
file.path(
OUT_DIR,
"CandidateGeneOverlap_Overall_CKD.csv"
)
)

summary_df <- data.frame(

Comparison="Overall CKD",

Total_DMPs=nrow(overall),

Hyper=sum(
overall$Methylation=="Hypermethylated"
),

Hypo=sum(
overall$Methylation=="Hypomethylated"
)

)

comparison_nodes <- paste0(
  summary_df$Comparison,
  "\nTotal CpGs (n=",
  summary_df$Total_DMPs,
  ")"
)

total_nodes <- "Distribution"

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

"Overall CKD",

"Total",

"Hyper",
"Hypo"

])

.range([


"#7570B3",   /* Overall CKD purple */

"#BDBDBD",   /* Total grey */

"#B2182B",   /* Hyper red */
"#2166AC"    /* Hypo blue */

])
'

nodes$group <- c(

"Comparison",

"Total",

"Hyper",

"Hypo"

)

links$group <- c(

"Overall CKD",

"Hyper",

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
"overall_CKD_Sankey.html"
),

selfcontained = TRUE
)

cat("Saved Sankey plot.\n")

# ============================================================
# STAGE 4 + STAGE 5 CANDIDATE GENE SANKEY
# ============================================================

library(networkD3)
library(htmlwidgets)

# ------------------------------------------------------------
# Read candidate gene overlap files
# ------------------------------------------------------------

stage4 <- read.csv(
  file.path(
    OUT_DIR,
    "CandidateGeneOverlap_Stage4_vs_Control.csv"
  ),
  stringsAsFactors = FALSE
)

stage5 <- read.csv(
  file.path(
    OUT_DIR,
    "CandidateGeneOverlap_Stage5_vs_Control.csv"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary_df <- data.frame(

  Comparison = c(
    "Stage4",
    "Stage5"
  ),

  Total_DMPs = c(
    nrow(stage4),
    nrow(stage5)
  ),

  Hyper = c(
    sum(stage4$Methylation == "Hypermethylated"),
    sum(stage5$Methylation == "Hypermethylated")
  ),

  Hypo = c(
    sum(stage4$Methylation == "Hypomethylated"),
    sum(stage5$Methylation == "Hypomethylated")
  )

)

print(summary_df)

# ------------------------------------------------------------
# Node labels
# ------------------------------------------------------------

comparison_nodes <- paste0(
  summary_df$Comparison,
  "\nTotal CpGs (n=",
  summary_df$Total_DMPs,
  ")"
)

total_nodes <- paste0(
  summary_df$Comparison,
  "\nDistribution"
)

hyper_nodes <- paste0(
  "Hyper\n(n=",
  summary_df$Hyper,
  ")"
)

hypo_nodes <- paste0(
  "Hypo\n(n=",
  summary_df$Hypo,
  ")"
)

# ------------------------------------------------------------
# Nodes
# ------------------------------------------------------------

nodes <- data.frame(

  name = c(
    comparison_nodes,
    total_nodes,
    hyper_nodes,
    hypo_nodes
  ),

  stringsAsFactors = FALSE

)

# ------------------------------------------------------------
# Links
# ------------------------------------------------------------

links <- data.frame()

for(i in 1:nrow(summary_df)){

  comp  <- comparison_nodes[i]
  total <- total_nodes[i]
  hyper <- hyper_nodes[i]
  hypo  <- hypo_nodes[i]

  links <- rbind(
    links,
    data.frame(
      source = match(comp, nodes$name)-1,
      target = match(total, nodes$name)-1,
      value = summary_df$Total_DMPs[i]
    )
  )

  links <- rbind(
    links,
    data.frame(
      source = match(total, nodes$name)-1,
      target = match(hyper, nodes$name)-1,
      value = summary_df$Hyper[i]
    )
  )

  links <- rbind(
    links,
    data.frame(
      source = match(total, nodes$name)-1,
      target = match(hypo, nodes$name)-1,
      value = summary_df$Hypo[i]
    )
  )

}

# ------------------------------------------------------------
# Node groups
# ------------------------------------------------------------

nodes$group <- c(

  rep("Comparison",2),

  rep("Total",2),

  rep("Hyper",2),

  rep("Hypo",2)

)

# ------------------------------------------------------------
# Link groups
# ------------------------------------------------------------

links$group <- c(

  "Stage4",
  "Hyper",
  "Hypo",

  "Stage5",
  "Hyper",
  "Hypo"

)

# ------------------------------------------------------------
# Colours
# ------------------------------------------------------------

my_color <- '

d3.scaleOrdinal()

.domain([

"Stage4",

"Stage5",

"Total",

"Hyper",

"Hypo"

])

.range([

"#1B9E77",

"#D95F02",

"#BDBDBD",

"#B2182B",

"#2166AC"

])

'

# ------------------------------------------------------------
# Draw Sankey
# ------------------------------------------------------------

p_sankey <- sankeyNetwork(

  Links = links,

  Nodes = nodes,

  Source = "source",

  Target = "target",

  Value = "value",

  NodeID = "name",

  NodeGroup = "group",

  LinkGroup = "group",

  sinksRight = FALSE,

  fontSize = 14,

  nodeWidth = 25,

  nodePadding = 35,

  colourScale = my_color

)

# ------------------------------------------------------------
# Save
# ------------------------------------------------library(circlize)------------

saveWidget(
  p_sankey,
  file.path(
    OUT_DIR,
    "Stage4_Stage5_CandidateGene_Sankey.html"
  ),
  selfcontained = TRUE
)

cat("Stage-specific Sankey saved.\n")
