# DrugCLIP 实施状态（更新至 2026-09-06）

## 给队友的准确口径

> DrugCLIP 项目脚本已经写好并完成两次端到端技术演示：可以把指定 PIEZO1 同结构共识口袋和 SMILES 分子表转换为 LMDB，调用官方权重在 RTX 5060 上检索，并导出带口袋归属的 CSV。当前只证明流程跑通，尚未选定正式分子库，演示排名不是候选药物结论。

## 已确认事项

- 上游源码位置：`third_party/DrugCLIP-main`。
- 固定提交：`7a3a3fa33673f8668c811790f2e4681c98af44ef`。
- 上游 `retrieval.sh` 只是写死路径的示例，不能直接处理本项目四个构象。
- DrugCLIP 检索需要三个外部输入：`checkpoint_best.pt`、分子库 `mols.lmdb`、口袋 `pocket.lmdb`。
- 官方 Google Drive 文件已核实：
  - `retrieval/mols.lmdb`：文件 ID `17OjIF7_HgDcC72m-nyLU26o5jcR9U801`
  - `retrieval/pocket.lmdb`：文件 ID `1CtpveHTBXCHnfcuIYoDkAAx1WXTJUJJF`
  - `checkpoint_best.pt`：文件 ID `1i87thnbNk8qeLF_tLx_BzelTukWbHaTR`，约 1.18 GB
- 官方 checkpoint 已完整下载：1,183,713,459 bytes，SHA256 为 `dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e`。
- 官方 `pocket.lmdb` 已下载至 `tools/drugclip/official/pocket.lmdb`（45,056 bytes）。官方演示 `mols.lmdb` 有约 4.68 GB，不是项目小库验证的必要条件，暂未继续下载。
- WSL `Ubuntu-fpocket` 已建立 `/opt/drugclip-venv`：Python 3.10.21、PyTorch 2.7.1+cu128、RDKit 2022.09.5、Uni-Core 0.0.1。
- RTX 5060 8 GB 的 CUDA 张量测试成功，设备能力为 12.0。
- 项目 `results/consensus/*_consensus_pockets.csv` 已包含 DrugCLIP 口袋生成所需的构象 ID、共识编号、链、残基编号和中心坐标。无需人工 Excel 转录。
- `8YEZ/C001`：19 个三工具共同残基、178 个口袋原子、5 个演示分子，端到端检索成功。
- `8YEZ/C002`：20 个三工具共同残基、180 个口袋原子，一键总流程复测成功。
- `8YEZ/C060`：3 个两工具共同残基、28 个口袋原子，T1 一键流程复测成功。

## 已实现脚本

- `scripts/10_prepare_drugclip_inputs.py`：共识口袋与分子表转 LMDB，记录输入/输出 SHA256。
- `scripts/11_run_drugclip.py`：校验单口袋约束和官方 checkpoint 哈希，运行 DrugCLIP，保留日志与环境版本。
- `scripts/12_normalize_drugclip_results.py`：把上游文本排名映射回分子 ID 并导出标准 CSV。
- `scripts/drugclip_pipeline_wsl.ps1`：把以上三步合成一个 Windows 命令。
- `scripts/bootstrap_drugclip_wsl.ps1`：建立固定版本的 WSL/GPU 运行环境。

## 下一阶段待办

1. 由医学/药学成员确定正式小分子库的来源、许可和筛选规模；演示分子不能直接沿用。
2. 由团队确定各结构先跑哪些 T2/T1 共识口袋，避免对全部 483 个区域盲目全量计算。
3. 对 8YEZ、8ZU3、8YFC、9VMX 分别独立运行正式 DrugCLIP 检索。
4. 根据同一分子在不同结构/口袋中的排名形成稳健性表，再进入 GNINA 对接。
5. 如比赛审计明确要求复现上游 retrieval 数据，再续传约 4.68 GB 的官方 `mols.lmdb`；它不阻塞项目侧脚本使用。

## 关键科学边界

- 四个结构必须独立运行，禁止把不同构象的口袋做交集。
- 当前项目配置是：8YEZ、8ZU3、8YFC、9VMX；8ZU8 仅参考。
- 当前仓库的等级名按队内最新口径为：T2 = 三工具共识、T1 = 两工具共识。
- DrugCLIP 分数只用于虚拟筛选排序，不等于真实结合、疗效或安全性证据；命中分子必须再经对接、MD/MM-GBSA、ADME/毒性和实验验证。

## 恢复工作时首先检查

```powershell
Set-Location 'D:\丘山\F_大二上学期\R_生物竞赛\piezo1-vs'
Get-Content .\docs\DRUGCLIP_HANDOFF_2026-09-05.md
powershell -ExecutionPolicy Bypass -File .\scripts\drugclip_pipeline_wsl.ps1 -PdbId 8YEZ -ConsensusId C001 -Compounds .\examples\drugclip_demo_compounds.csv
git status --short --branch
```
