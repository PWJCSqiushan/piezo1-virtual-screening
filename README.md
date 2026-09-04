# PIEZO1 Virtual Screening

第一届全球大学生生命科学挑战赛赛道三的可复现计算项目。当前目标是先完成最小闭环：

```text
四个正式结构：8YEZ / 8ZU3 / 8YFC / 9VMX
→ 每个结构分别用 DoGSite3 + fpocket + P2Rank 独立预测
→ 只在同一个结构、同一个编号体系内统计每个区域的工具支持数
→ 禁止跨构象求共同口袋
→ 8ZU8仅保留为参考，不进入正式计算
→ DrugCLIP 小库检索
→ GNINA 小规模 docking
→ SwissADME + ADMETlab 3.0 + ProTox 3.0 分层评估（非口服导向候选不自动删除）
→ 仅对最终 1–3 个候选评估膜环境 MD / MM-GBSA
```

本仓库只产生候选分子和计算证据，不能把 AI、docking、MD 或 MM-GBSA 结果表述为已验证药物、激动剂/抑制剂或 COPD 治疗效果。

## 当前阶段（2026-09-04）

- [x] 固化4个正式结构：8YEZ、8ZU3、8YFC、9VMX；8ZU8标记为仅参考
- [x] 下载 PDB/mmCIF 和 RCSB API 元数据并记录 SHA256（RCSB 验证 PDF 链接返回 404，已作为非阻塞警告记录）
- [x] 生成 `results/structure_manifest.csv`
- [x] 本地跑通 P2Rank 2.5.1：8YEZ 112 个；8ZU3、8YFC、9VMX各121个候选口袋
- [x] 跑通 DoGSite3 REST：8YEZ 112 个；8ZU3、8YFC、9VMX各100个候选口袋
- [x] 生成4个结构的P2Rank/DoGSite3双工具预检查：18/20、19/20、19/20、19/20
- [x] 实现跨平台P2Rank安装器、fpocket运行器和三工具同结构共识生成器，并通过合成测试
- [x] 记录突变位点的坐标可见性：8YFC的A1988V位点未出现在当前原子坐标中；9VMX的删除变体不能仅凭ATOM记录直接验证
- [ ] 跑通 fpocket（当前 Windows 缺少 Linux/容器环境）
- [ ] 在Linux/WSL/HPC获得真实fpocket结果后，对四个结构分别生成最终分层口袋

分类规则以 `config/consensus_rules.json` 为准。fpocket完成前，当前P2Rank/DoGSite3的18/20匹配只能称为预检查，不能提前标记T1。
- [x] 固定 DrugCLIP 官方源码版本和许可证
- [ ] DrugCLIP 官方示例烟雾测试（需要 Linux/WSL/HPC、checkpoint 和 Uni-Core 环境）
- [ ] GNINA 小规模 docking
- [x] 固化 Step 8 数据分流政策：SwissADME/ADMETlab作综合早筛，ProTox作主要计算毒性筛查；保留所有原始输出
- [ ] 候选分子产生后运行 Step 8，并导出供药学指导老师审核的全量表
- [ ] 仅对最终 1-3 个候选评估 MD/MM-GBSA

本次组会的正式修订、科学纠错和PPT证据清单见 `docs/MEETING_DECISIONS_2026-09-01_PART2.md`。机器可读的Step 8策略见 `config/compound_assessment.json`。

## 快速开始

```bash
git clone https://github.com/PWJCSqiushan/piezo1-virtual-screening.git
cd piezo1-virtual-screening
python -m venv .venv
```

Python核心代码只使用标准库。另需Java 17运行P2Rank；fpocket建议在Linux、WSL2或实验室HPC安装。

```bash
python scripts/01_fetch_structures.py
python scripts/02_build_structure_manifest.py
python scripts/bootstrap_p2rank.py
python scripts/04_run_p2rank.py --pdb-id 8YEZ
python scripts/03_run_dogsite3.py --pdb-id 8YEZ
python scripts/05_match_pockets.py --pdb-id 8YEZ
python scripts/08_verify_mvp.py --pdb-id 8YEZ
```

fpocket完成前可生成不带Tier的双工具预检查：

```bash
python scripts/09_build_consensus.py --pdb-id 8YEZ --allow-partial
```

在Linux/WSL/HPC安装[fpocket](https://github.com/Discngine/fpocket)后生成最终共识：

```bash
python scripts/run_fpocket.py --pdb-id 8YEZ
python scripts/09_build_consensus.py --pdb-id 8YEZ
```

如果fpocket由另一台机器运行，可传入其完整输出目录：

```bash
python scripts/09_build_consensus.py --pdb-id 8YEZ --fpocket-output /path/to/8YEZ_out
```

最终CSV/JSON位于`results/consensus/`，保留`support_count`、工具列表、原始口袋编号、中心、共同残基和来源文件。对其余三个PDB逐个重复以上命令，禁止把不同PDB放进同一次共识计算。

可选择演示任一已完成DoGSite3结果的结构：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/demo.ps1 -PdbId 8YEZ
powershell -ExecutionPolicy Bypass -File scripts/demo.ps1 -PdbId 8ZU3
powershell -ExecutionPolicy Bypass -File scripts/demo.ps1 -PdbId 8YFC
powershell -ExecutionPolicy Bypass -File scripts/demo.ps1 -PdbId 9VMX
```

所有脚本默认从仓库根目录运行，也可从任意目录调用。原始下载放在 `data/raw/`，单次运行放在 `runs/`，标准化汇总放在 `results/`。

### 测试

```bash
python -m compileall -q src scripts tests
python -m unittest discover -s tests -v
```

GitHub Actions会在push和Pull Request时自动执行语法、JSON配置和共识算法测试。

## 目录

```text
config/                 权威输入配置
data/raw/structures/    PDB/mmCIF、RCSB API 和验证报告
data/processed/         预处理结构与口袋文件
docs/                   实施计划、决策和组会材料
scripts/                一键入口
src/piezo_vs/           可复用 Python 代码
tools/                  项目内便携工具
runs/                   带时间戳的原始运行记录
results/                标准化 CSV/JSON 结果
tests/                  不依赖真实数据的轻量测试
```

## 可复现要求

每次运行都应保留输入文件 SHA256、软件版本、命令参数、UTC 时间、标准输出/错误和失败记录。禁止手工修改原始结果后覆盖；修正逻辑应进入脚本并生成新运行目录。

大型结构、运行目录、结果、第三方源码和模型权重默认不会提交到Git。协作规范见[CONTRIBUTING.md](CONTRIBUTING.md)，第三方工具及许可证核对清单见[docs/THIRD_PARTY.md](docs/THIRD_PARTY.md)，当前可公开的工具状态见[docs/STATUS.md](docs/STATUS.md)。

### 重要结构限制

8YFC和9VMX分别由条目元数据标识为A1988V与E756del变体，但当前8ZU3、8YFC、9VMX的ATOM/HETATM坐标记录相同；其中A1988V位点也未解析进原子坐标。因此这些条目可以按团队要求进入独立流程，却不能被表述为已经观察到突变引起的局部口袋差异。`structure_manifest`会输出坐标哈希和`mutation_coordinate_status`供报告与人工复核。

## 许可

本仓库原创代码采用[MIT License](LICENSE)。下载的结构、远程服务结果、第三方软件和模型权重继续受各自来源的条款约束。
