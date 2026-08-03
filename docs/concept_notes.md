# 细胞、基因、独立小鼠、原始counts、标准化矩阵、PCA的区别
## **细胞**是`adata.obs`，**基因**是`adata.var`，**独立小鼠**是`adata.obs['sample']`，**原始counts**是`adata.X`，**标准化矩阵**是对**原始counts**进行标准化处理得到的覆盖`adata.X`的矩阵，**PCA**是对**标准化矩阵**进行线性降维得到的(n_cells,n_PCs)矩阵储存在`adata.obsm['X_PCA']`和`adata.varm['PCs']`中
# Anndata
## `.X`是存储原始数据的稀疏矩阵  
## `.var` `.obs`对应genes和cells
## `.layers`是用于存储其他标准化或非标准化数据的结构
## `.obsm`是用于存储不同维度的对齐`n_obs`的数据的结构
# UCE(universal cell embedding)通用细胞嵌入模型
## 特点：可以整合多个单细胞RNA测序数据集，通用性好，针对不同来源的数据集都能适用而不用单独进行模型训练和调整，针对训练阶段未曾遇见的数据也能较好的分析，能够将所有细胞全部标注
## 模型建立：将基因根据表达情况进行加权，使用蛋白质语言模型将其通过其表达产物进行表征，之后将其结合其在染色体上的位置，一同送入Transformer神经网络进行训练。
## 输入：scRNA-seq数据
## 输出：样本集各细胞UCE空间向量
