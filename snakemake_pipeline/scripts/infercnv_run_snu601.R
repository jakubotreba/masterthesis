library(infercnv)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

load_matrix <- read.table(snakemake@input$matrix, sep = "\t", header = TRUE, row.names = 1)
load_matrix <- load_matrix[!duplicated(load_matrix[, 1]), ]

rownames(load_matrix) <- load_matrix[, 1]

load_matrix <- load_matrix[, -1]
annotation <- snakemake@input$annot
gene_ordering <- snakemake@input$gene_pos

output_file <- snakemake@output$cnv_file

load_matrix <- read.table("../snu601_data/SNU601_matrix_.txt", sep = "\t", header = TRUE, row.names = 1)

annotations <- "../snu601_data/SNU601_annotation.txt"

formatted_annotation <- read.table(annotations, header = FALSE, stringsAsFactors = FALSE)

gene_ordering <- "../snu601_data/hg38_gencode_v27.txt"

reference_cells <- NULL
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=load_matrix,
                                annotations_file = annotations,
                                delim="\t",
                                gene_order_file=gene_ordering,
                                ref_group_names = reference_cells)
infercnv_obj = infercnv::run(infercnv_obj,
                         cutoff=0.1,
                        out_dir = "./results/output_snu601/infercnv",
                         cluster_by_groups=TRUE,
                         denoise=TRUE,
                         HMM=T, 
                             HMM_type="i6",
                             analysis_mode="subclusters"
                         )
# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

