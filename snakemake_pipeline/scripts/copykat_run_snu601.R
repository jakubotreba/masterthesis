library(data.table)
library(copykat)


load_matrix <- read.table(snakemake@input$matrix, sep = "\t", header = TRUE, row.names = 1)

gene_expr_matrix <- load_matrix

annotations <- snakemake@input$annot

formatted_annotation <- read.table(annotations, header = FALSE, stringsAsFactors = FALSE)

gene_expr_matrix <- load_matrix

output_file <- snakemake@output$cnv_file

#####

if (ncol(formatted_annotation) < 2) {
    stop("Error: Annotation file does not have at least 2 columns. Check the format.")
}

# Extract normal cell UMIs correctly
normal_umis <- formatted_annotation[formatted_annotation[, 2] == "Normal", 1]
normal_cells <- as.character(normal_umis)

normal_cells <- NULL
## analysis
copykat.run <- copykat(rawmat=gene_expr_matrix, 
            id.type="S",
                    ngene.chr=5, 
                    LOW.DR = 0.05,
                    UP.DR = 0.1,
            win.size=25, 
                    KS.cut=0.1, 
            sam.name="SNU601", 
                    distance="euclidean", 
            norm.cell.names=normal_cells,
                    output.seg="FALSE", 
            plot.genes="TRUE", 
                    genome="hg20")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()


