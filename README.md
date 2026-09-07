# PIEZO1 Virtual Screening

面向第一届全球大学生生命科学挑战赛赛道三的可复现计算流程。仓库为每一步保留输入、版本、阈值、SHA256、原始输出和标准化结果，目标是产出**可追溯的候选数据**，而不是把计算预测包装成实验结论。

```text
4 个结构分别处理：8YEZ / 8ZU3 / 8YFC / 9VMX
→ DoGSite3 + fpocket + P2Rank 独立预测口袋
→ 同一结构内按空间距离和共同残基建立共识区域
→ T1 = 恰好 2 个工具支持；T2 = 3 个工具支持
→ 标准化口袋 → DrugCLIP 分子排序 → GNINA docking
→ SwissADME / ADMETlab 3.0 / ProTox 3.0 无损汇总与专家分流
→ 最终少量候选再评估膜环境 MD / MM-GBSA → 湿实验
```

> 科学边界：AI 排名、docking、MD 和 MM-GBSA 都只能生成候选与计算证据，不能证明分子已结合 PIEZO1、具有激动/抑制作用、能够治疗 COPD，或对人体安全有效。

## 现在能做什么

截至 2026-09-06，代码已经覆盖方案 2.0 的 Step 1–9 数据通路，并为 Step 8 提供可运行的导入/导出和专家复核队列；Step 10–12 仍依赖候选收敛、服务器、膜体系参数和实验条件。

| 环节 | 当前状态 | 已有可核查输出 |
|---|---|---|
| 结构与元数据 | 已跑通 | 4 个正式结构、哈希、链、分辨率和突变坐标状态 |
| 三工具口袋预测 | 已跑通 | P2Rank、DoGSite3、fpocket 原始结果与运行记录 |
| 同结构共识 | 已跑通 | 严格三工具组合 + 不可扩展双工具组合；重叠假设另存审计表 |
| 标准化口袋 | 已跑通 | 原始 PDB 坐标提取、残基/原子数、256 原子裁剪预警 |
| DrugCLIP | 已跑通 | 固定源码和权重、LMDB、GPU 排名、来源链接、输入输出哈希 |
| GNINA | 已跑通工程样例 | 固定二进制、受体/配体/搜索盒、SDF 姿势和标准化评分表 |
| ADME/毒性 | 接口已实现，真实预测待跑 | 全候选保留、原始字段留存、非口服项进入吸入/制剂专家复核 |
| 膜 MD / MM-GBSA | 待候选收敛 | 尚无可审计轨迹，不能宣称完成 |
| 湿实验 | 团队外部实验环节 | 不属于本仓库计算完成范围 |

四个结构的当前三工具结果：

| PDB | P2Rank | DoGSite3 | fpocket | 共识区域 | T1（2工具） | T2（3工具） |
|---|---:|---:|---:|---:|---:|---:|
| 8YEZ | 112 | 112 | 294 | 169 | 99 | 70 |
| 8ZU3 | 121 | 100 | 246 | 147 | 74 | 73 |
| 8YFC | 121 | 100 | 246 | 147 | 74 | 73 |
| 9VMX | 121 | 100 | 246 | 147 | 74 | 73 |

“共识区域”不是把三张表逐行硬配。程序先枚举满足中心距离 `≤12 Å` 且共同残基 `≥3` 的严格三工具组合；只有无法扩展成三工具组合的双工具配对才可成为 T1。空间上重复的假设再按固定代表分组，同时保留完整假设表，避免贪心一对一算法把三工具区域误降级。结果同时写出 `min_pairwise_jaccard_similarity`（越接近 1 越相似）和 `max_pairwise_jaccard_distance`（越接近 0 越相似），避免把相似度与距离的方向说反。

## 5 分钟本地演示

克隆后先看不需要 GPU 的代码与口袋流程：

```powershell
git clone https://github.com/PWJCSqiushan/piezo1-virtual-screening.git
cd piezo1-virtual-screening
python -m unittest discover -s tests -v
powershell -ExecutionPolicy Bypass -File scripts\demo.ps1 -PdbId 8YEZ
```

完整结果需要先生成或获得被 `.gitignore` 排除的大型运行数据。四个结构必须分别运行，禁止跨构象求共同口袋：

```powershell
python scripts\09_build_consensus.py --pdb-id 8YEZ
python scripts\08_verify_mvp.py --pdb-id 8YEZ
```

### 可复现实例：PubChem 128 分子工程压力测试

CID 1–128 是用于验证数据量、来源追踪和端到端输出的**顺序样本**，不是化学多样性库，也不是正式候选库。

```powershell
python scripts\13_fetch_pubchem_smoke_library.py --start-cid 1 --count 128 `
  --output runs\validation\pubchem_cid_1_128.csv

powershell -ExecutionPolicy Bypass -File scripts\bootstrap_drugclip_wsl.ps1
powershell -ExecutionPolicy Bypass -File scripts\drugclip_pipeline_wsl.ps1 `
  -PdbId 8YEZ -ConsensusId C005 `
  -Compounds runs\validation\pubchem_cid_1_128.csv
```

选择 C005 是因为它在当前 8YEZ 结果中由三个工具支持、属于 PIEZO1 范围，提取 18 个核心残基和 163 个原子，不触发 DrugCLIP 的 256 原子裁剪；它仍只是工程演示口袋，正式生物学优先级需要团队审核。

主要输出：

- `molecule_manifest.csv`：输入分子、规范 SMILES、PubChem CID/链接及拒绝原因；
- `standardized_pocket.pdb`：从同一个 8YEZ 原始 PDB 提取的标准化口袋；
- `input_run.json`：口袋、模型限制、随机种子和所有输入/输出哈希；
- `ranked_compounds.csv`：每个分子的口袋归属、相对名次和 DrugCLIP 分数；
- `report/index.html`：适合演示的自包含证据报告。

DrugCLIP 分数只允许在**同一口袋、同一分子库、同一模型设置**内做相对排序；它不是结合自由能。

同一 128 分子样本还可分别用 FP16/FP32 运行，再用 `scripts/19_compare_drugclip_runs.py` 输出 Spearman 排名相关系数和 Top-k 重合率；这是数值稳定性检查，不是生物学验证。

### GNINA 小规模 docking

首次运行会下载并校验固定的 GNINA 1.3.3 二进制，并在 Ubuntu 24.04 WSL 中安装官方 CUDA 12.8/cuDNN 9 运行库（总下载量约 4 GB）：

```powershell
powershell -ExecutionPolicy Bypass -File scripts\bootstrap_gnina_wsl.ps1

powershell -ExecutionPolicy Bypass -File scripts\gnina_pipeline_wsl.ps1 `
  -PdbId 8YEZ -ConsensusId C005 `
  -RankedCsv runs\drugclip\8YEZ\C005\<run-id>\ranked_compounds.csv `
  -TopN 3 -Exhaustiveness 2 -NumModes 3
```

程序会检查口袋范围、原子数、PDB/口袋归属、二进制大小与 SHA256；搜索盒超过 40 Å 时默认停止，不会静默裁剪。输出包含 `docked.sdf`、全部姿势表、每个分子的最佳姿势表、命令、版本和哈希。`CNNscore` 与 `minimizedAffinity` 是计算排序信号，不是实验结合证据。

### Step 8：ADME/毒性数据无损汇总

先把 DrugCLIP 或 GNINA 结果转换为外部工具输入：

```powershell
python scripts\17_prepare_compound_assessment.py `
  --ranked-csv runs\drugclip\8YEZ\C005\<run-id>\ranked_compounds.csv `
  --docking-csv runs\gnina\8YEZ\C005\<run-id>\best_poses.csv `
  --top-n 20 --output-dir runs\assessment\<run-id>
```

从 SwissADME、ADMETlab 3.0、ProTox 3.0 导出 CSV 后再合并。必须明确指出哪个 SwissADME 列是团队选定的口服优先启发式：

```powershell
python scripts\18_merge_compound_assessment.py `
  --master runs\assessment\<run-id>\compound_assessment_master.csv `
  --swissadme <swissadme.csv> --admetlab3 <admetlab3.csv> --protox3 <protox3.csv> `
  --oral-pass-column <exact-column-name> --output-dir runs\assessment_merged\<run-id>
```

所有工具原始列都会带前缀保留。口服启发式为 `No` 的分子不会被删除，而会进入 `inhalation_expert_review`；缺失、冲突或无法识别的结果进入 `general_expert_review`。

## 安装与数据说明

- Python 核心测试只依赖标准库；分子准备和结果解析使用隔离 WSL 环境中的 RDKit。
- P2Rank 需要 Java 17；fpocket 建议在 Linux/WSL2/HPC 使用。Windows 一键入口是 `scripts/bootstrap_fpocket_wsl.ps1`。
- DrugCLIP 环境为 PyTorch 2.7.1 + CUDA 12.8、RDKit 2022.09.5 和固定 Uni-Core 提交；入口是 `scripts/bootstrap_drugclip_wsl.ps1`。
- `data/`、`runs/`、`results/`、`tools/`、第三方源码和权重默认不进入 Git，新克隆不会自带大型结果。团队共享正式结果时应另行提供带哈希的归档。
- 版本和下载地址见 `config/tool_sources.json`，第三方许可见 [docs/THIRD_PARTY.md](docs/THIRD_PARTY.md)。
- 嘶，今天我啥也没干，哭了555

## 结构限制

8YFC 和 9VMX 分别由条目元数据标识为 A1988V 与 E756del 变体，但当前 8ZU3、8YFC、9VMX 的 ATOM/HETATM 坐标记录相同；A1988V 位点也未解析进 8YFC 原子坐标。它们按团队要求分别进入流程，但重复口袋结果**不能**解释为突变特异性效应。`results/structure_manifest.csv` 会输出坐标哈希和 `mutation_coordinate_status` 供复核。

## 质量门禁与目录

```powershell
python -m compileall -q src scripts tests
python -m unittest discover -s tests -v
python scripts\08_verify_mvp.py --pdb-id 8YEZ
```

```text
config/       权威结构、共识、评估和第三方版本配置
scripts/      可重复执行的一键入口
src/piezo_vs/ 可测试的解析、共识、报告、docking 和评估逻辑
examples/     明确标记为技术演示的输入
runs/         每次原始运行及日志（不进 Git）
results/      标准化 CSV/JSON 结果（不进 Git）
tests/        不依赖大型真实数据的快速回归测试
```

每次正式运行都应保留 UTC 时间、参数、输入/输出 SHA256、软件版本、标准输出/错误和失败记录。禁止手工修改原始结果后覆盖；算法修正应进入代码并生成新结果。

原创代码采用 [MIT License](LICENSE)。第三方软件、模型权重、结构和服务结果继续受各自来源条款约束。
