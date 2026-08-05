# 数据字典 (Metadata Dictionary)

> **数据来源**: Tabula Muris Senis (TMS) - Droplet
> **生成时间**: 2026-08-05 16:05
> **组织数量**: 6
> **Counts来源**: `.raw.X`（原始 `.X` 为 log-normalized）

## 1. 文件总览

| 组织 | 细胞数 | 基因数 | Counts位置 | 原始Counts | mouse_id |
|------|--------|--------|-----------|-----------|----------|
| spleen/marrow | 7,784 | 20,138 | .X (原始counts) | ✅ | ✅ |
| Fat | 3,368 | 20,138 | .X (原始counts) | ✅ | ✅ |
| Liver | 3,260 | 20,138 | .X (原始counts) | ✅ | ✅ |
| Lung | 18,334 | 20,138 | .X (原始counts) | ✅ | ✅ |
| Kidney | 1,729 | 20,138 | .X (原始counts) | ✅ | ✅ |
| Heart_and_Aorta | 6,138 | 20,138 | .X (原始counts) | ✅ | ✅ |

## 2. 字段定义与映射

| 标准字段 | 原始列名 | 位置 | 含义 | 关键约束 |
|----------|----------|------|------|----------|
| `mouse_id` | `mouse.id` | `obs` | 独立小鼠个体编号 | 统计重复单位；缺少则不可用于样本级统计 |
| `age_months` | `age` | `obs` | 小鼠月龄 | 需建立 age_group；分组边界第6周冻结 |
| `sex` | `sex` | `obs` | 性别 | pseudobulk模型协变量 |
| `tissue` | `tissue` | `obs` | 组织来源 | 需标准化为 tissue_std |
| `cell_type` | `cell_ontology_class` | `obs` | 细胞类型（本体论分类） | 保留三列: cell / cell_ontology_class / free_annotation |
| `assay` | `method` | `obs` | 测序技术/方法 | 可能与年龄混杂，需作为批次协变量 |
| `counts` | `.raw.X` → `.X` | `.X` | 原始UMI表达矩阵 | 必须为原始整数counts，非log-normalized |

## 3. 逐组织字段详情

### spleen/marrow

- 细胞数: 7,784
- 基因数: 20,138
- Counts: ✅ 原始counts (dtype=float32)

| 标准字段 | 实际列名 | 类型 | 唯一值 | 缺失率 | 示例值 |
|----------|----------|------|--------|--------|--------|
| mouse_id | mouse.id | category | 19 | 0.0% | 18-F-50, 18-F-51, 21-F-54 |
| age_months | age | category | 6 | 0.0% | 18m, 21m, 24m |
| sex | sex | category | 2 | 0.0% | female, male |
| tissue | tissue | category | 1 | 0.0% | spleen/marrow |
| cell_type | cell_ontology_class | category | 4 | 0.0% | macrophage, macrophage dendritic cell progenitor, monocyte |
| assay | method | category | 1 | 0.0% | droplet |

### Fat

- 细胞数: 3,368
- 基因数: 20,138
- Counts: ✅ 原始counts (dtype=float32)

| 标准字段 | 实际列名 | 类型 | 唯一值 | 缺失率 | 示例值 |
|----------|----------|------|--------|--------|--------|
| mouse_id | mouse.id | category | 6 | 0.0% | 18-F-50, 18-M-52, 18-M-53 |
| age_months | age | category | 3 | 0.0% | 18m, 21m, 30m |
| sex | sex | category | 2 | 0.0% | female, male |
| tissue | tissue | category | 1 | 0.0% | Fat |
| cell_type | cell_ontology_class | category | 2 | 0.0% | endothelial cell, mesenchymal stem cell of adipose |
| assay | method | category | 1 | 0.0% | droplet |

### Liver

- 细胞数: 3,260
- 基因数: 20,138
- Counts: ✅ 原始counts (dtype=float32)

| 标准字段 | 实际列名 | 类型 | 唯一值 | 缺失率 | 示例值 |
|----------|----------|------|--------|--------|--------|
| mouse_id | mouse.id | category | 12 | 0.0% | 18-F-51, 21-F-54, 24-M-58 |
| age_months | age | category | 6 | 0.0% | 18m, 21m, 24m |
| sex | sex | category | 2 | 0.0% | female, male |
| tissue | tissue | category | 1 | 0.0% | Liver |
| cell_type | cell_ontology_class | category | 3 | 0.0% | endothelial cell of hepatic sinusoid, Kupffer cell, hepatic stellate cell |
| assay | method | category | 1 | 0.0% | droplet |

### Lung

- 细胞数: 18,334
- 基因数: 20,138
- Counts: ✅ 原始counts (dtype=float32)

| 标准字段 | 实际列名 | 类型 | 唯一值 | 缺失率 | 示例值 |
|----------|----------|------|--------|--------|--------|
| mouse_id | mouse.id | category | 16 | 0.0% | 18-F-50, 18-F-51, 18-M-52 |
| age_months | age | category | 5 | 0.0% | 18m, 21m, 30m |
| sex | sex | category | 2 | 0.0% | female, male |
| tissue | tissue | category | 1 | 0.0% | Lung |
| cell_type | cell_ontology_class | category | 13 | 0.0% | non-classical monocyte, alveolar macrophage, classical monocyte |
| assay | method | category | 1 | 0.0% | droplet |

### Kidney

- 细胞数: 1,729
- 基因数: 20,138
- Counts: ✅ 原始counts (dtype=float32)

| 标准字段 | 实际列名 | 类型 | 唯一值 | 缺失率 | 示例值 |
|----------|----------|------|--------|--------|--------|
| mouse_id | mouse.id | category | 16 | 0.0% | 18-F-50, 18-F-51, 18-M-52 |
| age_months | age | category | 6 | 0.0% | 18m, 21m, 24m |
| sex | sex | category | 2 | 0.0% | female, male |
| tissue | tissue | category | 1 | 0.0% | Kidney |
| cell_type | cell_ontology_class | category | 3 | 0.0% | macrophage, kidney capillary endothelial cell, fibroblast |
| assay | method | category | 1 | 0.0% | droplet |

### Heart_and_Aorta

- 细胞数: 6,138
- 基因数: 20,138
- Counts: ✅ 原始counts (dtype=float32)

| 标准字段 | 实际列名 | 类型 | 唯一值 | 缺失率 | 示例值 |
|----------|----------|------|--------|--------|--------|
| mouse_id | mouse.id | category | 11 | 0.0% | 18-F-50, 18-M-52, 18-M-53 |
| age_months | age | category | 6 | 0.0% | 18m, 21m, 24m |
| sex | sex | category | 2 | 0.0% | female, male |
| tissue | tissue | category | 1 | 0.0% | Heart_and_Aorta |
| cell_type | cell_ontology_class | category | 3 | 0.0% | fibroblast of cardiac tissue, smooth muscle cell, endothelial cell of coronary artery |
| assay | method | category | 1 | 0.0% | droplet |

## 4. 辅助列（未标准化但保留）

| 列名 | 说明 |
|------|------|
| `cell` | 原始细胞类型标签（比 cell_ontology_class 更细） |
| `cell_ontology_id` | 细胞本体论 ID |
| `free_annotation` | 自由注释（研究者手动标注） |
| `subtissue` | 子组织区域 |
| `tissue_free_annotation` | 组织自由注释 |
| `n_genes` | 每个细胞检测到的基因数（QC指标） |
| `n_counts` | 每个细胞的总UMI数（QC指标） |
| `louvain` | Louvain聚类结果 |
| `leiden` | Leiden聚类结果 |

## 5. 关键发现

1. **Counts位置**: 原始 `.X` 为 log-normalized 数据（非负但非整数），原始 counts 存放在 `.raw.X` 中。已提取 `.raw.X` 覆盖保存至 `../data_interim/`。
2. **mouse_id**: 所有组织均存在 `mouse.id` 列，可用于样本级统计。
3. **age**: 以月龄（months）为单位，需后续离散化为 age_group。
4. **cell_type**: 主注释使用 `cell_ontology_class`，同时保留 `cell`（原始标签）和 `free_annotation`（自由注释）用于交叉验证。
5. **assay/method**: 所有数据均为 droplet (10x Chromium)，批次效应主要来自不同小鼠和不同组织。
