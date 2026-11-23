# X-Chromosome DNA Methylation and Cognitive Variation in Centenarians
## 1. Project overview

研究对象为 Dutch 100-plus Study 百岁老人。
使用 whole-blood PacBio methylation 数据。
研究问题为 chrX CpG methylation 与 MMSE 的关联。

本 Repository 覆盖三个 Step、七个 Python 脚本。

## 2. Repository structure

| 目录 | 文件 | 用途 |
|---|---|---|
| 根目录 | `README.md` | 项目说明及运行顺序 |
| 根目录 | `environment.yml` | 记录实际 `xci_0` Environment |
| 根目录 | `.gitignore` | 排除本地数据及 Output |
| `Step_1/` | `make_centenarian_sample_list.py` | Build sample list |
| `Step_1/` | `make_chrX_methylation_matrix_centenarian.py` | Build methylation matrix |
| `Step_1/` | `make_chrX_coverage_matrix_centenarian.py` | Build coverage matrix |
| `Step_2/` | `Coverage_and_CpG_completeness_exploration.py` | 比较 QC thresholds |
| `Step_2/` | `make_chrX_usable_matrix.py` | Apply QC，Export usable matrix |
| `Step_3/` | `run_primary_association.py` | Run Primary analysis |
| `Step_3/` | `sex_stratified_chrX.py` | Run Female、Male analyses |

`Input/` 和 `Output/` 为本地运行目录。
实际样本数据及运行产物不 Push 至 GitHub。

`environment.yml`

`.gitignore`

## 3. Software environment

使用 `xci_0` Environment。

主要 Dependency：

| Dependency | 用途 |
|---|---|
| Python | Run 脚本 |
| pandas | Read、Join 和 Export 表格 |
| NumPy | Matrix 操作 |
| SciPy | 回归及统计检验 |
| Matplotlib | QC plots |

脚本未记录 Dependency 版本。
实际版本须从 `xci_0` Export。

运行准备：

1. Activate `xci_0`。
2. Export 实际 Environment。
3. Update 各脚本中的 paths。
4. Create 本地运行目录。
5. 按 Step_1 → Step_2 → Step_3 顺序 Run。

```bash
conda activate xci_0
conda env export --no-builds > environment.yml

mkdir -p Step_1/Input Step_1/Output
mkdir -p Step_2/Output
mkdir -p Step_3/Input Step_3/Output
```

上述命令在 Repository 根目录 Run。

当前脚本使用 hardcoded paths，无 CLI parameters。
默认工作目录：

```text
/project/holstegelab/Share/yixuan/workflow
```

各脚本的 paths 必须指向同一工作目录。

## 4. Input data

### 4.1 原始 Input

| Input | 必需 columns / 内容 |
|---|---|
| `Step_1/Input/pb_samples_phenotypes_20240425.txt` | TSV；`pacbio_id`、`ID_GWAS`、`ID_100plus`、`PHENOTYPE`、`SEX` |
| `100plus_mmse.csv` | `studysubjectid`、`harmonized_contact_nr`、`mmse_score_imputed` |
| 每个样本的 `*.cpg.combined.bed` | 第 1、2、4、6 columns：chromosome、start、methylation、coverage |
| `Step_3/Input/final_model_samples_correct_age_cellprops.tsv` | `pacbio_id`、`SEX`、`mmse_score_imputed`、`PC1`–`PC10`、education、age、六种 cell proportions |

MMSE 默认路径：

```text
/project/holstegelab/Share/yaran/sample/100plus_mmse.csv
```

BED 默认路径：

```text
/project/holstegelab/Share/yaran/epigenetics/phased/denovo/{pacbio_id}/{pacbio_id}.cpg.combined.bed
```

### 4.2 ID 对应与格式要求

1. 使用 `pacbio_id` Align 样本与 matrix。
2. 使用 `ID_100plus ↔ studysubjectid` Join MMSE。
3. `SEX` 统一使用 `F`、`M`。
4. 各样本表中的 `pacbio_id` 必须唯一。
5. 每个样本的 chrX `start` 必须唯一。
6. MMSE contact 0 中，`studysubjectid` 必须唯一。

### 4.3 Step_3 covariates

| 模型变量 | Input column / 规则 |
|---|---|
| MMSE | `mmse_score_imputed` |
| Education | 优先 `education`；缺失该 column 时用 `verhage_education_level_numeric` |
| Age | 优先 `age`；缺失该 column 时用 `age_entrance_interview` |
| Cell proportions | `CD8T`、`CD4T`、`NK`、`Bcell`、`Mono`、`Neu` |
| B-cell alias | 缺失 `Bcell` column 时，可用 `Bcel` |

如果 Input 已含 `age`，则其值必须来自 `age_entrance_interview`。

七个脚本不 Build Step_3 covariate table。
该文件须事先准备。

### 4.4 Matrix 格式

两个原始 matrix 均为 **sample × CpG**。

前三个 columns：

```text
pacbio_id    ID_GWAS    ID_100plus
```

其余 columns 为：

```text
chrX:{start}
```

CpG coordinate 直接沿用 BED start。
脚本不执行 `+1` 转换。

**Methylation 单位必须一致。**

Step_1 原样 Read BED 第 4 column。
Step_3 将其除以 100。
因此，Input 必须使用 0–100 的 methylation 百分数。

## 5. Step_1 — Sample selection and matrix construction

按以下顺序 Run。
本节 Output 均位于 `Step_1/Output/`。

### 5.1 `make_centenarian_sample_list.py`

**用途：** Build 百岁老人 sample list。

**Input：** phenotype TSV、MMSE CSV、BED 文件路径。

处理步骤：

1. 保留 `PHENOTYPE == "Centenarian"`。
2. Exclude `pacbio_id` 含 `brain` 的样本。
3. Exclude BED 文件不存在的样本。
4. Join `harmonized_contact_nr == "0"` 的 MMSE。
5. Export sample list。

**Output：** `centenarian_sample_list.tsv`。

Output columns：

```text
pacbio_id
ID_GWAS
ID_100plus
SEX
mmse_score_imputed
bed_path
```

此步保留 MMSE 缺失的样本。

### 5.2 `make_chrX_methylation_matrix_centenarian.py`

**用途：** Build 原始 chrX methylation matrix。

**Input：** sample list 及其 `bed_path` 指向的 BED。

处理步骤：

1. Collect 所有样本的 chrX CpG 并集。
2. 按 coordinate 排序。
3. 用 BED 第 4 column Build matrix。
4. 未观测值保留 NaN。
5. Export methylation matrix。

**Output：** `chrX_methylation_matrix.tsv`。

### 5.3 `make_chrX_coverage_matrix_centenarian.py`

**用途：** Build 与 methylation 对齐的 coverage matrix。

**Input：** sample list、methylation matrix、BED。

处理步骤：

1. 沿用 methylation matrix 的 sample 和 CpG 顺序。
2. 用 BED 第 6 column Build coverage matrix。
3. Validate coverage 为非负、有限整数。
4. Export coverage matrix。

**Output：** `chrX_coverage_matrix.tsv`。

两个 matrix 的 rows 和 columns 一一对应。

### 5.4 运行顺序

```bash
python Step_1/make_centenarian_sample_list.py
python Step_1/make_chrX_methylation_matrix_centenarian.py
python Step_1/make_chrX_coverage_matrix_centenarian.py
```

## 6. Step_2 — Coverage quality control

两个脚本均先 Exclude 样本 `8783353`。
本节 Output 均位于 `Step_2/Output/`。

### 6.1 `Coverage_and_CpG_completeness_exploration.py`

**用途：** 比较 coverage 与 completeness thresholds。

**Input：** coverage matrix、sample list。

处理步骤：

1. Exclude `8783353`。
2. Verify 剩余样本数为 247。
3. 比较两组 coverage thresholds。
4. 比较 completeness：50%、80%、90%。
5. Export QC plots 和 summary。

Coverage 设置：

| 模式 | Female | Male |
|---|---|---|
| Same absolute threshold | 5、10、15、20× | 5、10、15、20× |
| Sex-adjusted threshold | 5、10、15、20× | 2.5、5、7.5、10× |

Completeness 为达标人数除以全部 247 人。
NaN coverage 不达标。

**Output：**

| 文件 | 内容 |
|---|---|
| `QC_coverage_thresholds_247samples.png` | 各样本达标的 CpG 数 |
| `QC_completeness_thresholds_247samples.png` | 各 threshold 组合的保留 CpG 数 |
| `QC_threshold_summary_247samples.tsv` | 24 个组合的统计结果 |

### 6.2 `make_chrX_usable_matrix.py`

**用途：** Apply 最终 QC，Export 正式分析 matrix。

**Input：** methylation matrix、coverage matrix、sample list。

处理步骤：

1. Exclude `8783353`。
2. Apply Female ≥10×、Male ≥5×。
3. Retain coverage completeness ≥80% 的 CpGs。
4. 将低 coverage 对应的 methylation Mask 为 NaN。
5. 将 NaN coverage 对应的 methylation Mask 为 NaN。
6. 保留原有 methylation NaN。
7. Export 正式 NaN matrix。

**正式 Output：**

```text
chrX_usable_matrix_10X_80pct_NaN.tsv
```

以 247 人计算时，至少需 198 人 coverage 达标。

Completeness 基于 coverage。
它不等于非缺失 methylation 比例。

**当前代码差异：**

脚本随后仍执行 column-mean imputation。
它另行 Export：

```text
chrX_usable_matrix_10X_80pct_imputed.tsv
```

该文件不用于正式 association。
它不会覆盖已 Export 的 NaN matrix。

此脚本使用实际排除后人数计算 completeness。
它未固定 Verify 分母为 247。

### 6.3 运行顺序

```bash
python Step_2/Coverage_and_CpG_completeness_exploration.py
python Step_2/make_chrX_usable_matrix.py
```

## 7. Step_3 — Association analysis

两个脚本均直接 Read：

```text
Step_2/Output/chrX_usable_matrix_10X_80pct_NaN.tsv
Step_3/Input/final_model_samples_correct_age_cellprops.tsv
```

共同处理步骤：

1. 按 `pacbio_id` Align 两个 Input。
2. Exclude 缺失 MMSE 或 covariates 的样本。
3. 将 methylation 百分数除以 100。
4. 每个 CpG 使用该位点的 complete cases。
5. Fit 回归模型。
6. Export 全量结果和候选结果。

不对 methylation 作 imputation。
本节 Output 均位于 `Step_3/Output/`。

### 7.1 `run_primary_association.py`

合并 Female 与 Male，Run Primary model：

```text
MMSE ~ CpG methylation + sex + age + education
       + PC1–PC10
       + CD8T + CD4T + NK + Bcell + Mono + Neu
```

模型含 intercept。
Sex 编码为 Female=0、Male=1。

代码要求 MMSE 与 covariates 完整的样本数为 235。

**Output：**

| 文件 | 内容 |
|---|---|
| `primary_all_CpGs.tsv` | 所有成功估计的 CpGs |
| `Table_1_primary_suggestive_CpGs.tsv` | `P < 1e-5` 的候选，按 P 排序 |

### 7.2 `sex_stratified_chrX.py`

分别在 Female 与 Male 中 Run：

```text
MMSE ~ CpG methylation + age + education
       + PC1–PC10
       + CD8T + CD4T + NK + Bcell + Mono + Neu
```

模型含 intercept，不含 sex。

代码要求 Input 含 1,159,396 个 CpG columns。

**Output：**

| 文件 | 内容 |
|---|---|
| `female_all_CpGs.tsv` | Female 全部 CpGs |
| `male_all_CpGs.tsv` | Male 全部 CpGs |
| `Table_female_suggestive_CpGs.tsv` | Female `P < 1e-5` 候选 |
| `Table_male_suggestive_CpGs.tsv` | Male `P < 1e-5` 候选 |

Stratified 全量表保留无法估计的 CpG。
这些 CpG 的统计值为 NaN。

Primary 全量表不保留无法估计的 CpG。

### 7.3 结果 columns

| Column | 含义 |
|---|---|
| `CpG position` | CpG ID |
| `beta` | 调整 covariates 后的 methylation 回归系数；使用 0–1 fraction 尺度 |
| `SE` | beta 的 standard error |
| `P-value` | 双侧 t-test 的 unadjusted P |
| `N` | 该 CpG 的实际 complete-case 样本数 |

N 随 CpG 改变，不等于固定 cohort size。

235 和 1,159,396 为代码中的检查值。

### 7.4 运行顺序

```bash
python Step_3/run_primary_association.py
python Step_3/sex_stratified_chrX.py
```

## 8. Analysis settings

| 项目 | 正式分析设置 |
|---|---|
| Sample QC | Exclude `8783353` |
| Coverage | Female ≥10×；Male ≥5× |
| Completeness | ≥80%；目标分母 247 |
| Methylation missing values | 使用 NaN matrix；不作 imputation |
| Association | 每个 CpG 使用 complete cases |
| Age | `age_entrance_interview` |
| Primary covariates | sex、age、education、PC1–PC10、六种 cell proportions |
| Stratified covariates | 与 Primary 相同，去掉 sex |
| 候选 threshold | Unadjusted `P < 1 × 10⁻⁵` |
| 执行方式 | 一次 Load 全部 CpGs，逐位点 Fit；不分 chunk |
| 不纳入 | variance filtering、FDR、leave-one-out、效应换算 |


回归会跳过样本不足或无法估计的位点。
这不属于预先的 variance filtering。

## 9. Data and code availability

代码范围为 Step_1–Step_3。
不含 covariate table 的构建或后续 annotation。
不增加 `.sh`、独立 plotting 脚本或 Step_4。

原始 BED、phenotype、covariates 和样本级 matrix 不随代码公开。
