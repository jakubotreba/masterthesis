library(data.table)
library(copykat)

# Logging
log <- file(snakemake@log[[1]], open="wt")
ink(log)
sink(log, type="message")

# Read input data
load_matrix <- read.table(snakemake@input$matrix, sep = "\t", header = TRUE, row.names = 1)
annotations <- snakemake@input$annot
output_file <- snakemake@output$cnv_file

# Convert gene expression matrix to numeric
gene_expr_matrix <- as.data.frame(lapply(load_matrix, as.numeric))
rownames(gene_expr_matrix) <- rownames(load_matrix)  # Keep rownames intact
#####
load_matrix <- read.table("../tnbc1_data/GSM4476486_combined_UMIcount_CellTypes_TNBC1.txt", sep = "\t", header = TRUE, row.names = 1)

annotations <- "../tnbc1_data/TNBC1_copykat_annotation.txt"

formatted_annotation <- read.table(annotations, header = FALSE, stringsAsFactors = FALSE)

gene_expr_matrix <- load_matrix
#####
if (ncol(formatted_annotation) < 2) {
    stop("Error: Annotation file does not have at least 2 columns. Check the format.")
}

# Extract normal cell UMIs correctly
normal_umis <- formatted_annotation[formatted_annotation[, 2] == "Normal", 1]
normal_cells <- as.character(normal_umis)

# Run CopyKat analysis
copykat.run <- copykat(rawmat=gene_expr_matrix, 
                    id.type="S",
                    cell.line="no",
                    ngene.chr=5, 
                    LOW.DR = 0.05,
                    UP.DR = 0.1,
                    win.size=25, 
                    KS.cut=0.1, 
                    sam.name="TNBC1", 
                    distance="euclidean", 
                    norm.cell.names=normal_cells,
                    output.seg="FALSE", 
                    plot.genes="TRUE", 
                    genome="hg20",
                    n.cores=1)