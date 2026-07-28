# 细胞、基因、独立小鼠、原始counts、标准化矩阵、PCA的区别
## **细胞**是`adata.obs`，**基因**是`adata.var`，**独立小鼠**是`adata.obs['sample']`，**原始counts**是`adata.X`，**标准化矩阵**是对**原始counts**进行标准化处理得到的覆盖`adata.X`的矩阵，**PCA**是对**标准化矩阵**进行线性降维得到的(n_cells,n_PCs)矩阵储存在`adata.obsm['X_PCA']`和`adata.varm['PCs']`中
# Anndata
## `.X`是存储原始数据的稀疏矩阵  
## `.var` `.obs`对应genes和cells
## `.layers`是用于存储其他标准化或非标准化数据的结构
## `.obsm`是用于存储不同维度的对齐`n_obs`的数据的结构