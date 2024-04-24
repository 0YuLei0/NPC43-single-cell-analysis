# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import anndata

# set a working directory
wdir = "/import/home/lyuah/Py_content/npc_patient_senic"
os.chdir( wdir )
f_loom_path_scenic = "input_scenic.loom"
## cdna and SCT
counts=pd.read_csv("/import/home/lyuah/R/npc_patient/counts.csv")
counts.rename(columns={counts.columns[0]: "gene" }, inplace = True)
counts.index = counts['gene'].values
del counts['gene']
counts = counts.sort_index()

## Load gene annotation
gene = pd.DataFrame(counts.index)
gene.rename(columns={gene.columns[0]: "gene" }, inplace = True)
gene.index = gene['gene'].values

## Creat anndata
adata = anndata.AnnData(X = counts.transpose(),
                        #obs = anno,
                        var = gene,
                        )
# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)

################ in command line 
arboreto_with_multiprocessing.py input_scenic.loom ~/GBM/pysenic/database/hs_hgnc_tfs.txt --method genie3 --output genie3_adj.tsv --seed 66 --num_workers 20
pyscenic ctx genie3_adj.tsv ~/GBM/pysenic/database/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname ~/GBM/pysenic/database/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname input_scenic.loom --mode "dask_multiprocessing" --output genie3_reg.csv --mask_dropouts
pyscenic aucell input_scenic.loom genie3_reg.csv --output scenic_genie3_out.loom --num_workers 20
