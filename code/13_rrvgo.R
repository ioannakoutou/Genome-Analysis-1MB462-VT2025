if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("rrvgo", quietly = TRUE)) {
  BiocManager::install("rrvgo")
}

if (!requireNamespace("org.At.tair.db", quietly = TRUE)) {
  BiocManager::install("org.At.tair.db")
}

if (!requireNamespace("GOstats", quietly = TRUE)) {
  BiocManager::install("GOstats")
}

if (!requireNamespace("GSEABase", quietly = TRUE)) {
  BiocManager::install("GSEABase")
}

library(GOstats)
library(GSEABase)
library(rrvgo)
library(dplyr)
library(stringr)
library(tidyr)

setwd("/Users/ioannakoutou/Documents/2nd Semester/Genome Analysis - 1MB462/Project_Files/DEA/")

# Load gff
gff <- read.delim("eggnog_anot.emapper.decorated.gff", sep = "\t", header  = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(gff) <- c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Load DESeq2 results
res_df <- read.delim("deseq2_all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sig_genes <- res_df$GeneName[ res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 2 ]
gene_universe <- res_df$GeneName

# GO terms from gff
go_df <- gff %>%
  filter(type == "mRNA", source == "AUGUSTUS") %>%
  mutate(
    gene_id = str_extract(attributes, "Parent=[^;]+") %>% str_remove("Parent="),
    go_terms = str_extract(attributes, "GO:[^;]+") %>% str_extract_all("GO:\\d+")
  ) %>%
  unnest(go_terms) %>%
  filter(!is.na(go_terms)) %>%
  transmute(
    go_id = as.character(go_terms),
    Evidence = "IEA",
    gene_id = as.character(gene_id)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

goframe <- GOFrame(go_df, organism = "N. japonicum")
goallframe <- GOAllFrame(goframe)
gsc <- GeneSetCollection(goallframe, setType = GOCollection())

# hyperGTest parameters
params <- GSEAGOHyperGParams(
  name = "GOstats enrichment",
  geneSetCollection = gsc,
  geneIds           = sig_genes,
  universeGeneIds   = gene_universe,
  ontology          = "BP",
  pvalueCutoff      = 0.05,
  conditional = FALSE,
  testDirection     = "over"
)

# Gene enrichment
go_result <- hyperGTest(params)
go_table_stats <- summary(go_result)

# Run rrvgo
simMatrix <- calculateSimMatrix(
  go_table_stats$GOBPID,
  orgdb = "org.At.tair.db",
  ont = "BP",
  method = "Rel"
)

scores <- setNames(-log10(go_table_stats$Pvalue), go_table_stats$GOBPID)

reducedTerms <- reduceSimMatrix(
  simMatrix,
  scores,
  orgdb = "org.At.tair.db",
  threshold = 0.7
)

# Plots
p_scatter <- scatterPlot(simMatrix, reducedTerms, addLabel = TRUE)
ggsave("GO_scatter.png", p_scatter, width = 8, height = 6,bg = "white")

png("BP_GO_treemap.png", width = 8, height = 6, units  = "in", res = 300)
treemapPlot(reducedTerms)
dev.off()

png(filename = "GO_wordcloud.png", width = 8, height = 6, units    = "in", res = 300)
wordcloudPlot(reducedTerms, min.freq = 1, colors = "black")
dev.off()