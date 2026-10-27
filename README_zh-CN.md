# Genetical Mapping

当前版本：**1.1**。

Genetical Mapping 是一款面向基因家族整体序列刻画与样本比较的工具。它汇集同一目标家族的成员序列，建立重叠短片段库，统计片段在成员集合中的出现情况，生成样本层面的特征矩阵、聚类树和热图，帮助比较不同样本的家族序列特征。

本命令行版读取 FASTA 格式的核苷酸序列，导出特征矩阵、聚类结果和片段统计表。

[English](README.md)

## 主要功能

- 整体分析输入的基因家族成员集合，同一样本可包含多条成员序列。
- 将各成员的重叠短片段汇入共同的 k-mer 库，统计片段在全部输入成员中的出现次数。
- 按样本和序列位置汇总片段统计，或生成样本层面的 k-mer 频率矩阵。
- 按片段总出现次数选择分析范围，比较不同频次范围下的家族特征。
- 对位置特征进行行归一化，或合并互为反向互补的 k-mer。
- 导出分组结果、矩阵和聚类诊断；位置模式同时生成交互式热图。
- 记录运行参数、输入校验值、序列数量和依赖版本。

## 1.1 版

本次修订提供独立命令行版、明确的特征表示与记录筛选选项、离线交互热图、输入输出检查及运行清单，并随附合成示例和回归测试。

运行 `python main.py --version` 可查看版本。更新概况见 [CHANGELOG.md](CHANGELOG.md)。

## 安装与快速开始

使用 Python 3.12。在 Windows PowerShell 中进入解压后的 `Genetical-Mapping` 目录，运行：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py --fasta examples/synthetic_demo.fasta --out outputs/demo --k 20 --clusters 2
```

在浏览器中打开 `outputs/demo/interactive_heatmap.html`，即可查看热图和树状图。HTML 已内嵌绘图库，可离线查看。

示例包含 4 个演示组的 8 条人工序列。运行 `python generate_demo.py` 可用固定种子 20261027 重新生成示例。

## 输入格式

每次分析输入一个已确定的目标基因家族成员集合，同一样本可提供多个成员。默认的 `nonempty` 模式保留每条非空成员记录，包括序列内容相同的不同记录；其他 `--variant` 选项提供可选的序列筛选规则。

输入为包含 A、C、G、T 的 FASTA 文件。程序统一转换为大写。每条序列使用唯一标识符，第一个下划线之前的部分作为组名：

```text
>DemoA_sequence1 type:complete
ACGTTGCAACGTTGCA
>DemoB_sequence1 type:complete
ACGTCGCAACGTCGCA
```

这段内容用于说明标题格式。快速开始命令使用随软件提供的完整示例。

k-mer 长度可设为 1–64。筛选后每个组都需要保留可用序列窗口。输入至少包含 3 个组，聚类数介于 2 和组数减 1 之间。输出目录使用新目录或空目录。

## 常用参数

| 参数 | 用途 | 默认值 |
| --- | --- | --- |
| `--fasta` | 输入 FASTA 文件 | 必填 |
| `--out` | 输出目录 | 必填 |
| `--k` | k-mer 长度，1–64 | `20` |
| `--clusters` | 聚类数 | `2` |
| `--lower` | k-mer 总出现次数下限，含端点 | `1` |
| `--upper` | k-mer 总出现次数上限，含端点 | 全部次数 |
| `--representation` | 特征表示方式，见下文 | `position_raw` |
| `--canonical` | 合并互为反向互补的 k-mer | 关闭 |
| `--variant` | 序列筛选规则，见下文 | `nonempty` |
| `--legacy-end-drop` | 在位置模式中将每组最后一个可用位置置零 | 关闭 |

三种表示方式：

- `position_raw`：对每组每个起始位置出现的不同 k-mer，累加其在全部输入中的总出现次数。
- `position_l1`：将上述矩阵的每行除以该行总和。
- `kmer_frequency`：将每组 k-mer 次数除以筛选前的有效窗口总数，以频率平方根的欧氏距离比较各组。

两个位置模式使用欧氏距离。三种模式均采用完全连接法进行层次聚类。位置以输入序列中的零基起始坐标表示，位置模式直接比较各条未对齐序列的输入位置。

序列筛选规则：

- `nonempty`：保留全部非空序列。
- `exact_unique`：在每个组内按完整序列去重。
- `complete_only`：保留 FASTA 标题带有 `type:complete` 的序列。
- `complete_unique`：先按完整标记筛选，再在组内去重。
- `gene_representative`：从标识符末尾去掉 `_i<数字>` 或 `_i<数字>.p<数字>` 后分组；依次优先选择带完整标记的序列、较长序列、标识符字典序较前的序列。

例如，使用频率模式并合并反向互补 k-mer：

```powershell
.\.venv\Scripts\python.exe main.py --fasta examples/synthetic_demo.fasta --out outputs/frequency --k 20 --clusters 2 --representation kmer_frequency --canonical
```

运行 `.\.venv\Scripts\python.exe main.py --help` 可查看命令行选项。

## 输出文件

| 文件 | 内容 |
| --- | --- |
| `cluster_assignments.tsv` | 组名与聚类标签 |
| `distance_matrix.tsv` | 组间距离矩阵 |
| `feature_matrix.tsv` | 位置特征矩阵 |
| `feature_matrix.npz` | 频率模式的稀疏特征矩阵 |
| `linkage.npy` | 层次聚类连接矩阵 |
| `retained_kmers.tsv` | 保留的 k-mer 及总出现次数 |
| `position_support.tsv` | 每个起始位置的序列支持数 |
| `run_manifest.json` | 参数、数量、输入校验值、聚类诊断与环境版本 |
| `interactive_heatmap.html` | 位置模式的离线交互式热图和树状图 |

## 运行测试

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s . -p "test_*.py" -v
```

## 许可

采用 MIT 许可证，详见 [LICENSE](LICENSE)。各依赖保留其自身许可证。
