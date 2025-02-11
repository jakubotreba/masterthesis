library(devtools)
library(SCEVAN)

 Logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Read input data
load_matrix <- snakemake@input$load_matrix
annotations <- snakemake@input$annot

output_file <- snakemake@output$cnv_file
output_pred <- snakemake@output$pred_file


load_matrix <- read.table("../tnbc1_data/GSM4476486_combined_UMIcount_CellTypes_TNBC1.txt", sep = "\t", header = TRUE, row.names = 1)

gene_expr_matrix <- load_matrix[-c(1, 2), ]
gene_names <- rownames(gene_expr_matrix)
gene_expr_matrix <- as.data.frame(lapply(gene_expr_matrix, as.numeric))
rownames(gene_expr_matrix) <- gene_names

annotations <- "../tnbc1_data/TNBC1_scevan_annotation.txt"

formatted_annotation <- read.table(annotations, header = FALSE, stringsAsFactors = FALSE)

# Ensure the annotation file has the correct number of columns
if (ncol(formatted_annotation) < 2) {
    stop("Annotation file does not have at least 2 columns. Check the format.")
}

# Extract normal cells
normal_umis <- formatted_annotation[formatted_annotation[, 2] == "Normal", 1]
normal_cells <- as.character(normal_umis)

# Run SCEVAN pipeline
results <- SCEVAN::pipelineCNA(gene_expr_matrix,
                               sample = "TNBC1",
                               norm_cell = normal_cells,
                               SUBCLONES = TRUE,
                               par_cores = 1,
                               plotTree = FALSE)

# Save results correctly (avoid absolute path issue)
#write.table(results, file=output_pred, sep="\t", quote=FALSE)

# ------------------------------------------------------------------------------
print("SessionInfo:")
sessionInfo()
