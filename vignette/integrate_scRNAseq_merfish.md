## Integrate scRNAseq data with spatial MERFISH data

We show how to integrate high-plex RNA imaging-based spatially resolved MERFISH data with scRNA-seq data using [uniPort](https://www.biorxiv.org/content/10.1101/2024.03.19.585796v1). This integration tutorial is adapted from Uniport. The MERFISH data includes 64,373 cells with 155 genes, and the scRNA-seq data includes 30,370 cells with 21,030 genes from six mice measured with dissociated scRNA-seq (10X).

```
import uniport as up
import scanpy as sc
import pandas as pd
import numpy as np
```

We can use scanpy to read MERFISH and scRNA-seq data. 

```
adata_merfish = sc.read_h5ad('MERFISH/merfish0.h5ad')
adata_rna = sc.read_h5ad('MERFISH/rna0.h5ad')
```

However, for this tutorial purpose, we read MERFISH and scRNA-seq into AnnData objects using load_file fucntion in uniport for convenience. Then, we use 'domain_id' to identify the modality using a number category and ‘source’ to identify the modality using its name.

```
adata_merfish = up.load_file('MERFISH/MERFISH_mouse1.txt')
adata_merfish.obs['domain_id'] = 0
adata_merfish.obs['domain_id'] = adata_merfish.obs['domain_id'].astype('category')
adata_merfish.obs['source'] = 'MERFISH'

adata_rna = up.load_file('MERFISH/RNA_count.txt')
adata_rna.obs['domain_id'] = 1
adata_rna.obs['domain_id'] = adata_rna.obs['domain_id'].astype('category')
adata_rna.obs['source'] = 'RNA'
```

Optionally, we can also read the cell types for better visualization later. 
```
labels_merfish = pd.read_csv('MERFISH/MERFISH_st_filter_cluster.txt', sep='\t')
celltype_merfish = labels_merfish['cluster_main'].values
adata_merfish.obs['cell_type'] = celltype_merfish

labels_rna = pd.read_csv('MERFISH/MERFISH_scRNA_filter_cluster.txt', sep='\t')
celltype_rna = labels_rna['cluster_main'].values
adata_rna.obs['cell_type'] = celltype_rna
```

According to the recommended pipeline in Uniport, we preprocess data using functions normalize_total, log1p and highly_variable_genes in scanpy and batch_scale in uniport 
```
sc.pp.normalize_total(adata_merfish)
sc.pp.log1p(adata_merfish)
sc.pp.highly_variable_genes(adata_merfish, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(adata_merfish)

sc.pp.normalize_total(adata_rna)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(adata_rna, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(adata_rna)

sc.pp.normalize_total(adata_cm)
sc.pp.log1p(adata_cm)
sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000, batch_key='domain_id', inplace=False, subset=True)
up.batch_scale(adata_cm)
```

At first, concatenate MERFISH and scRNA-seq with common genes using AnnData.concatenate. 
```
adata_cm = adata_merfish.concatenate(adata_rna, join='inner', batch_key='domain_id')
```
Then, we use the advanced technique of integrating using dataset-specific genes by Run() function in uniport. The latent representations of data are stored in adata.obs['latent'].

```
adata = up.Run(adatas=[adata_merfish, adata_rna], adata_cm=adata_cm, lambda_kl=5.0)
```

To impute genes using the integrated scRNA-seq data, use the same Run() command but with additional parameters: out='predict' to project data into the latent space, and pred_id=1 to set the modality to scRNA-seq (as it has domain id 1). Then, save the data at the desired location.

````
adata_impute = up.Run(adata_cm=spatial_data_partial, out='predict', pred_id=1)
adata_impute.write('data/MERFISH_imputed/MERFISH_and_RNA.h5ad', compression='gzip')
````

We recommend users to see the original Uniport site for the details on gene imputation by integrating scRNAseq data. 

Then we call NEST's preprocess, run, post-process, and visualization steps as follows:
```
nest preprocess --data_name='MERFISH_imputed' --data_from='data/MERFISH_imputed/' --tissue_position_file='data/MERFISH_imputed/MERFISH_tissue_positions_list.csv'
```

```
nohup nest run --data_name='MERFISH_imputed'  --num_epoch 60000 --model_name='NEST_MERFISH_imputed' --run_id=1 > output_MERFISH_imputed_run1.log &
nohup nest run  --data_name='MERFISH_imputed'  --num_epoch 60000 --model_name='NEST_MERFISH_imputed' --run_id=2 > output_MERFISH_imputed_run2.log &
nohup nest run  --data_name='MERFISH_imputed'  --num_epoch 60000 --model_name='NEST_MERFISH_imputed' --run_id=3 > output_MERFISH_imputed_run3.log &
nohup nest run  --data_name='MERFISH_imputed'  --num_epoch 60000 --model_name='NEST_MERFISH_imputed' --run_id=4 > output_MERFISH_imputed_run4.log &
nohup nest run  --data_name='MERFISH_imputed'  --num_epoch 60000 --model_name='NEST_MERFISH_imputed' --run_id=5 > output_MERFISH_imputed_run5.log &
```

```
nest postprocess --data_name='MERFISH_imputed' --model_name='NEST_MERFISH_imputed' --total_runs=5 
```

```
nest visualize --data_name='MERFISH_imputed' --model_name='NEST_MERFISH_imputed'
```
