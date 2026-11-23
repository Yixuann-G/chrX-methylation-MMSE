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

