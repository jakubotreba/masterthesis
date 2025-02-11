import xclone
import numpy as np
import pandas as pd
import sys
import logging

log_file = snakemake.log[0]

logging.basicConfig(
    filename=log_file,                  
    level=logging.INFO,                 
    format='%(asctime)s - %(levelname)s - %(message)s',  
    filemode='w'                        
)

sys.stdout = open(log_file, 'a')        
sys.stderr = open(log_file, 'a')        

RDR_file = snakemake.input[0]
mtx_barcodes_file = snakemake.input[1]
regions_anno_file = snakemake.input[2]
gene_expr_matrix = pd.read_table(snakemake.input[3], sep="\t")
AD_file = snakemake.input[4]
DP_file = snakemake.input[5]


RDR_adata = xclone.pp.xclonedata(RDR_file, 'RDR', mtx_barcodes_file, genome_mode = "hg38_genes")

BAF_adata = xclone.pp.xclonedata([AD_file, DP_file], 'BAF',
                                 mtx_barcodes_file,
                                 genome_mode = "hg38_genes")

hg38_genes = xclone.pp.load_anno(genome_mode = "hg38_genes")
hg38_blocks = xclone.pp.load_anno(genome_mode = "hg38_blocks")
dataset = "TNBC1"
out_dir = "results/output_tnbc1/xclone/"
copykat_pred_row = gene_expr_matrix.loc['copykat.pred']
RDR_adata_test = xclone.data.tnbc1_rdr()
clusters = RDR_adata_test.obs.cluster
normalized_indices = RDR_adata.obs.index.str.replace("-1$", "", regex=True)
RDR_adata.obs['Tumor'] = normalized_indices.map(copykat_pred_row)
normalized_indices = BAF_adata.obs.index.str.replace("-1$", "", regex=True)
BAF_adata.obs['Tumor'] = normalized_indices.map(copykat_pred_row)
RDR_adata.obs['cluster'] = clusters.reindex(RDR_adata.obs.index)
BAF_adata.obs['cluster'] = clusters.reindex(BAF_adata.obs.index)
RDR_adata.obs['cluster'] = RDR_adata.obs.index.map(clusters)
BAF_adata.obs['cluster'] = BAF_adata.obs.index.map(clusters)
## analysis RDR module
xconfig = xclone.XCloneConfig(dataset_name = dataset, module = "RDR")
xconfig.set_figure_params(xclone= True, fontsize = 18)
xconfig.outdir = out_dir
xconfig.cell_anno_key = "Tumor"
xconfig.ref_celltype = "N"
xconfig.marker_group_anno_key = "Tumor"
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "cluster"
xconfig.exclude_XY = True
xconfig.display()
RDR_Xdata = xclone.model.run_RDR(RDR_adata,
            config_file = xconfig)
## analysis BAF module
xconfig = xclone.XCloneConfig(dataset_name = dataset, module = "BAF")
xconfig.set_figure_params(xclone= True, fontsize = 18)
xconfig.outdir = out_dir
xconfig.cell_anno_key = "Tumor"
xconfig.ref_celltype = "N"
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "cluster"
xconfig.exclude_XY = True
xconfig.display()
BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata,
            config_file = xconfig)
rdr_chr = RDR_Xdata.var["chr"].drop_duplicates().reset_index(drop = True)
baf_chr = BAF_merge_Xdata.var["chr"].drop_duplicates().reset_index(drop = True)
## analysis combine modul
xconfig = xclone.XCloneConfig(dataset_name = dataset, module = "Combine")
xconfig.set_figure_params(xclone= True, fontsize = 18)
xconfig.outdir = out_dir
xconfig.cell_anno_key = "Tumor"
xconfig.ref_celltype = "N"
xconfig.copygain_correct= False
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "cluster"
xconfig.merge_loss = False
xconfig.merge_loh = True
xconfig.BAF_denoise = True
xconfig.display()
combine_Xdata = xclone.model.run_combine(RDR_Xdata,
                BAF_merge_Xdata,
                verbose = True,
                run_verbose = True,
                config_file = xconfig)
















