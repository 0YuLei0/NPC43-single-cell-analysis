import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
import os
import urllib.request
import anndata
import palantir
import scprep

warnings.simplefilter(action="ignore", category=Warning)
# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
# sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300)

###########################
# Load gene annotation
scale_gene = pd.read_csv("/import/home/lyuah/R/npc43_trajectory/scale_gene.csv")
scale_gene = scale_gene.drop(scale_gene.columns[0], axis=1)
scale_gene.rename(columns={scale_gene.columns[0]: "gene" }, inplace = True)
scale_gene.index = scale_gene['gene'].values
# SCT scaled data (Seurat SCTransform normlzied Data)
scale_data=pd.read_csv("/import/home/lyuah/R/npc43_trajectory/scale_data.csv")
scale_data.index = scale_gene.index
scale_data = scale_data.transpose()
## load sample annotation
anno = pd.read_csv("/import/home/lyuah/R/npc43_trajectory/metadata.csv")
anno.index = anno['cellid'].values
# Load all gene annotation
gene = pd.read_csv("/import/home/lyuah/R/npc43_trajectory/gene.csv")
gene = gene.drop(gene.columns[0], axis=1)
gene.rename(columns={gene.columns[0]: "gene" }, inplace = True)
gene.index = gene['gene'].values
# counts (counts from raw scRNA-seq data)
counts=pd.read_csv("/import/home/lyuah/R/npc43_trajectory/count_data.csv")
counts.index = gene.index
counts = counts.transpose()
counts_filted = counts.loc[scale_data.index, scale_data.columns]
counts_filted = counts_filted.astype('float64')
## Creat anndata
adata = anndata.AnnData(X = counts_filted,
                        obs = anno,
                        var = scale_gene
                        )
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_per_cell(adata)
palantir.preprocess.log_transform(adata)
sc_counts_filted = adata.X
adata.layers["scanpy_scaled_counts"] = sc_counts_filted
adata.X = scale_data
## add NTRK2 and SOX2 expression rank
expression_df = sc.get.obs_df(adata, keys=["NTRK2","SOX2"], layer = "scanpy_scaled_counts")
expression_df["Batch"] = adata.obs["Batch3"]
expression_df = expression_df.sort_values(by=["NTRK2","SOX2"], ascending=False)
expression_df = expression_df[expression_df['Batch'] == 'UT']
adata.obs['highlight'] = adata.obs.index.isin(expression_df.index[0:499]).astype(int)
###### SCANPY with SCT normalized dataset from seurat v5
## PCA with normalized data from SCT seurat v5
sc.tl.pca(adata, n_comps=100,svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca(
    adata,
    color=["Batch", "Batch", "cell_annotation", "cell_annotation"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)
###### SCANPY with SCT normalized dataset from seurat v5
## PCA with normalized data from SCT seurat v5
sc.tl.pca(adata, n_comps=100,svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca(
    adata,
    color=["Batch", "Batch", "cell_annotation", "cell_annotation"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)
## Magic imputation; Optional
dm_res = palantir.utils.run_diffusion_maps(adata, n_components=50)
ms_data = palantir.utils.determine_multiscale_space(adata)
adata.X = sc_counts_filted
imputed_X = palantir.utils.run_magic_imputation(adata,n_jobs=16)
## Dimension reduction with UMAP and FA; neighbors parameters is same as in R
adata.X = scale_data
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=13,metric = "cosine", random_state = 65) # repeat R UMAP in python
## neighbors with denosie using diffmap rather than PCA
sc.tl.diffmap(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50, metric = "cosine", random_state = 65, use_rep="X_diffmap") #sc.pp.neighbors(adata, n_neighbors=30, n_pcs=13,metric = "cosine", random_state = 65, use_rep="X_diffmap")
sc.tl.umap(adata, min_dist=0.1, spread=1)
sc.pl.umap(adata, color="Batch", legend_loc="on data",size = 25, save="leiden_umap_plot.png")
sc.pl.umap(adata, color=["cell_batch_annotation","NTRK2","SOX2"], legend_loc="on data",size = 25,legend_fontsize= 'xx-small',layer="MAGIC_imputed_data")
sc.pl.umap(adata, edges=True, color = 'cell_batch_annotation', legend_loc= 'on data', legend_fontsize= 'xx-small',size = 25)
## Better to use leiden cluster rather than cellular states, as cellular states sometimes are too sparsely distributed
sc.tl.leiden(adata, resolution=1)
sc.tl.paga(adata, groups="leiden")
sc.pl.paga(adata, color='cell_batch_annotation', edge_width_scale = 0.3)
sc.tl.umap(adata, init_pos="paga", min_dist=0.1, spread=1)
sc.pl.umap(adata, color = ["leiden","cell_batch_annotation","highlight","NTRK2","SOX2"], legend_loc= 'on data' ,size = 30,legend_fontsize= 'xx-small',layer="MAGIC_imputed_data")
sc.pl.embedding(adata, color=["leiden","cell_batch_annotation","NTRK2","SOX2"], legend_loc="on data",size = 25,legend_fontsize= 'xx-small',layer="MAGIC_imputed_data", basis = "2d_umap")
## Use dpt to estimate the pseudotime
adata.uns["iroot"] = np.flatnonzero(adata.obs["highlight"] == 1)[0]
sc.tl.dpt(adata, n_dcs = 50)
sc.pl.umap(adata, color=["cell_batch_annotation", "dpt_pseudotime"], legend_loc="on data",size = 25,legend_fontsize= 'xx-small')
## plot flow based on DPT
import cellrank as cr
pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix()
sc.settings.set_figure_params(dpi=100)
pk.plot_projection(basis="X_umap", recompute=True,color="cell_batch_annotation",legend_fontsize= 'xx-small')
