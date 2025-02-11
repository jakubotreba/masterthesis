 library(infercnv)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

load_matrix <- read.table(snakemake@input$matrix, sep = "\t", header = TRUE, row.names = 1)
annotation <- snakemake@input$annot
gene_ordering <- snakemake@input$gene_pos

output_file <- snakemake@output$cnv_file

gene_expr_matrix <- load_matrix[-c(1, 2), ]
gene_names <- rownames(gene_expr_matrix)
gene_expr_matrix <- as.data.frame(lapply(gene_expr_matrix, as.numeric))
rownames(gene_expr_matrix) <- gene_names
reference_cells <- c("Normal")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=gene_expr_matrix,
                                annotations_file = annotation,
                                delim="\t",
                                gene_order_file=gene_ordering,
                                ref_group_names = reference_cells)
infercnv_obj = infercnv::run(infercnv_obj,
                         cutoff=0.1,
                        out_dir = "./results/output_tnbc1/infercnv",
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

