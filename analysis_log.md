# 分析日志：00_anndata_smoke_test.ipynb

> **项目**：AnnData 对象基础操作 Smoke Test
> **Notebook**：`00_anndata_smoke_test.ipynb`
> **作者**：暴璐玮
> **开始日期**：2026-07-28
> **最后更新**：2026-07-28

---

## 1. 环境信息

- **运行平台**：Google Colab（Python 3.12）
- **记录日期**：2026-07-28

### 核心包

| 包名 | 版本 |
|------|------|
| Python | 3.12 |
| scanpy | 1.12.3 |
| anndata | 0.13.2 |
| numpy | 2.4.6 |
| pandas | 2.3.3 |
| scipy | 1.16.3 |
| scikit-learn | 1.6.1 |
| matplotlib | 3.10.0 |
| seaborn | 0.13.2 |
| statsmodels | 0.14.6 |

### 降维 / 聚类相关

| 包名 | 版本 |
|------|------|
| umap-learn | 0.5.12 |
| pynndescent | 0.6.0 |
| numba | 0.66.0 |
| llvmlite | 0.48.0 |
| networkx | 3.6.1 |

### IO / 存储相关

| 包名 | 版本 |
|------|------|
| h5py | 3.16.0 |
| zarr | 3.2.1 |
| numcodecs | 0.16.5 |

### 其他依赖

| 包名 | 版本 |
|------|------|
| joblib | 1.5.3 |
| tqdm | 4.67.3 |
| pillow | 11.3.0 |
| pydantic | 2.13.4 |
| typing-extensions | 4.16.0 |
| packaging | 26.2 |

---

## 2. 数据来源

| 项目 | 内容 |
|------|------|
| 数据集 | 10x Genomics PBMC 3k |
| 获取方式 | `sc.datasets.pbmc3k()` |
| 下载大小 | 5.58 MB |
| 细胞数 | 2700 |
| 基因数 | 32738 |
| 参考 | https://support.10xgenomics.com/single-cell-gene-expression/datasets |

---

## 3. 分析步骤

### Step 1：安装依赖 & 导入包

- **日期**：2026-07-28
- **操作**：
  ```python
  !pip install scanpy
  !pip install AnnData
  !pip install numpy
  !pip install pandas

  import scanpy as sc
  import anndata as ad
  import numpy as np
  import pandas as pd
  ```
- **结果**：所有包已预装，无需额外下载 ✅

---

### Step 2：下载并保存原始数据

- **日期**：2026-07-28
- **操作**：
  ```python
  adata = sc.datasets.pbmc3k()
  adata.write_h5ad("pbmc3k_raw.h5ad")
  ```
- **输出文件**：`pbmc3k_raw.h5ad`
- **结果**：下载成功（5.58 MB），保存为 h5ad 格式 ✅

---

### Step 3：读取 h5ad 文件

- **日期**：2026-07-28
- **操作**：
  ```python
  adata = ad.read_h5ad("pbmc3k_raw.h5ad")
  print(adata)
  ```
- **输出**：
  ```
  AnnData object with n_obs × n_vars = 2700 × 32738
      var: 'gene_ids'
      layers: None (.X)
  ```

---

### Step 4：查看 AnnData 对象结构

- **日期**：2026-07-28
- **操作**：对 AnnData 各属性进行详细检查

| 属性 | 内容 |
|------|------|
| `adata.obs` | Empty DataFrame，shape (2700, 0)，无列 |
| `adata.var` | 1 列：`gene_ids`（ENSG ID） |
| `adata.obsm` | 空（keys: []） |
| `adata.varm` | 空（keys: []） |
| `adata.uns` | 空（keys: []） |
| `adata.X` | `scipy.sparse._csr.csr_matrix`，shape (2700, 32738)，dtype float32 |
| **稀疏度** | **97.41%** |

- **var 前 5 行示例**：

| index | gene_ids |
|-------|----------|
| MIR1302-10 | ENSG00000243485 |
| FAM138A | ENSG00000237613 |
| OR4F5 | ENSG00000186092 |
| RP11-34P13.7 | ENSG00000238009 |
| RP11-34P13.8 | ENSG00000239945 |

---

### Step 5：QC 指标计算 + 标准化 + HVG + PCA

- **日期**：2026-07-28
- **操作**：
  ```python
  # QC 指标
  sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

  # 标准化
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)

  # 高变基因
  sc.pp.highly_variable_genes(adata, n_top_genes=2000)

  # PCA
  sc.tl.pca(adata, n_comps=50)
  ```
- **QC 后 obs 列**：`['n_genes_by_counts', 'total_counts']`
  - 注：未计算线粒体比例（未指定 `qc_vars='mt'`）
- **PCA 后 obsm keys**：`['X_pca']`
- **X_pca shape**：(2700, 50) ✅

---

### Step 6：筛选 total_counts 最高的 100 个细胞

- **日期**：2026-07-28
- **操作**：
  ```python
  top100_idx = adata.obs['total_counts'].nlargest(100).index
  adata_100 = adata[top100_idx].copy()
  ```
- **筛选策略**：按 `total_counts` 降序取前 100
- **决策理由**：Smoke test 目的为验证 AnnData 操作流程，取高表达细胞可确保数据非空、后续操作有意义
- **结果验证**：

| 属性 | 值 |
|------|-----|
| n_obs | 100 ✅ |
| n_vars | 32738（基因数不变）✅ |
| obsm keys | `['X_pca']` ✅ |
| X_pca shape | (100, 50) ✅ |
| obs 列 | `n_genes_by_counts`, `total_counts` |
| var 列 | `gene_ids`, `n_cells_by_counts`, `mean_counts`, `pct_dropout_by_counts`, `total_counts`, `highly_variable`, `means`, `dispersions`, `dispersions_norm` |
| uns keys | `log1p`, `hvg`, `pca` |
| varm keys | `PCs` |

- **关键发现**：`.copy()` 后所有附属信息（obs、var、obsm、varm、uns）均正确跟随 ✅

---

### Step 7：保存 & 验证输出文件

- **日期**：2026-07-28
- **操作**：
  ```python
  output_path = "pbmc3k_100cells.h5ad"
  adata_100.write_h5ad(output_path)

  # 验证
  adata_check = ad.read_h5ad(output_path)
  print(adata_check)  # 100 × 32738 ✅

  import os
  size_mb = os.path.getsize(output_path) / 1024 / 1024
  print(f"文件大小: {size_mb:.2f} MB")
  ```
- **输出文件**：`pbmc3k_100cells.h5ad`
- **文件大小**：18.12 MB
- **验证结果**：重新读取后确认为 100 × 32738，结构完整 ✅

---

## 4. 输出文件清单

| 文件名 | 内容 | 大小 |
|--------|------|------|
| `pbmc3k_raw.h5ad` | 原始 PBMC 3k 数据（2700 cells × 32738 genes） | ~5.58 MB |
| `pbmc3k_100cells.h5ad` | 筛选后 100 细胞子集（含 QC + 标准化 + HVG + PCA） | 18.12 MB |
| `00_anndata_smoke_test.ipynb` | 分析 Notebook | — |
| `analysis_log.md` | 本日志 | — |

---

## 5. 遇到的问题 & 解决方案

| 日期 | 问题 | 原因 | 解决方案 |
|------|------|------|---------|
| 2026-07-28 | `ImportError: cannot import name '_center' from 'numpy._core.umath'` | numpy 内部文件版本不一致（部分 2.x 部分 1.x 残留） | `!pip install --force-reinstall --no-cache-dir numpy` + 重启 Colab 运行时 |

---

## 6. 结论

✅ **Smoke Test 全部通过**：

1. 成功下载 PBMC 3k 数据并保存为 h5ad 格式
2. 成功读取 h5ad 文件，确认 AnnData 结构完整
3. 详细检查了 `X`、`obs`、`var`、`obsm`、`varm`、`uns` 各属性
4. 完成 QC → 标准化 → HVG → PCA 全流程
5. 成功筛选 100 个细胞（按 total_counts 降序），所有附属信息正确跟随
6. 输出文件可正常读写，数据完整性验证通过

---

## 7. 参考文献

1. Luecken MD, Theis FJ. *Current best practices in single-cell RNA-seq analysis: a tutorial.* Molecular Systems Biology, 2019.
2. Scanpy 官方教程：https://scanpy.readthedocs.io/en/stable/tutorials.html
3. 10x Genomics PBMC 3k：https://support.10xgenomics.com/single-cell-gene-expression/datasets
4. Zheng GXY et al. *Massively parallel digital transcriptional profiling of single cells.* Nature Communications, 2017.
```
