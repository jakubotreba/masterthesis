library(devtools)
library(SCEVAN)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

load_matrix <- snakemake@input$load_matrix

input_clones <- as.logical(snakemake@params$find_clones)

output_file <- snakemake@output$cnv_file
output_pred <- snakemake@output$pred_file

load_matrix <- read.table("../snu601_data/SNU601_matrix_.txt", sep = "\t", header = TRUE, row.names = 1)

normal_cells <- NULL

# Run SCEVAN pipeline
results <- SCEVAN::pipelineCNA(load_matrix,
                               sample = "SNU601",
                               norm_cell = normal_cells,
                               SUBCLONES = TRUE,
                               par_cores = 1,
                               plotTree = FALSE)

# Save expected outputs
write.table(results, file=output_file, sep="\t", quote=FALSE)
write.table(results, file=output_pred, sep="\t", quote=FALSE)

print("SessionInfo:")
sessionInfo()
