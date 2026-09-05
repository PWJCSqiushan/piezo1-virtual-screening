# DrugCLIP 阶段交接（2026-09-05）

## 给队友的准确口径

> DrugCLIP 官方源码和原始检索脚本已经固定到项目中，但项目自己的批量脚本尚未完成验收。目前已经确认了 PIEZO1 共识口袋到 DrugCLIP 输入所需的数据字段，并取得官方测试文件清单；下一步是完成环境、输入 LMDB 转换、官方样例冒烟测试和 CSV 排名导出。因此现在不能说“DrugCLIP 已经跑通并产出候选药物”。

## 已确认事项

- 上游源码位置：`third_party/DrugCLIP-main`。
- 固定提交：`7a3a3fa33673f8668c811790f2e4681c98af44ef`。
- 上游 `retrieval.sh` 只是写死路径的示例，不能直接处理本项目四个构象。
- DrugCLIP 检索需要三个外部输入：`checkpoint_best.pt`、分子库 `mols.lmdb`、口袋 `pocket.lmdb`。
- 官方 Google Drive 文件已核实：
  - `retrieval/mols.lmdb`：文件 ID `17OjIF7_HgDcC72m-nyLU26o5jcR9U801`
  - `retrieval/pocket.lmdb`：文件 ID `1CtpveHTBXCHnfcuIYoDkAAx1WXTJUJJF`
  - `checkpoint_best.pt`：文件 ID `1i87thnbNk8qeLF_tLx_BzelTukWbHaTR`，约 1.18 GB
- 官方 `pocket.lmdb` 已下载至 `tools/drugclip/official/pocket.lmdb`（45,056 bytes）。
- 官方 `mols.lmdb` 因 Google Drive 下载限流尚未取得。
- checkpoint 下载已按用户关机要求中止；留下约 3.67 MB 的临时分片 `tools/drugclip/official/checkpoint_best.ptg22d9as9.part`，不是可用模型文件。
- 本机 WSL `Ubuntu-fpocket` 可识别 RTX 5060 8 GB GPU；但 DrugCLIP/Uni-Core 旧依赖是否兼容仍需以官方样例实测为准。
- 项目 `results/consensus/*_consensus_pockets.csv` 已包含 DrugCLIP 口袋生成所需的构象 ID、共识编号、链、残基编号和中心坐标。无需人工 Excel 转录。

## 尚未完成（恢复后按顺序执行）

1. 在 WSL 中建立隔离的 DrugCLIP Python 环境，安装 PyTorch、RDKit、LMDB、Uni-Core 等依赖。
2. 重新下载完整 `checkpoint_best.pt`，并取得官方 `mols.lmdb`；记录 SHA256。
3. 用官方 `mols.lmdb + pocket.lmdb + checkpoint` 跑通一次官方 retrieval 冒烟测试。
4. 编写项目输入转换脚本：
   - 从 SMILES/分子表生成 DrugCLIP `mols.lmdb`；
   - 从指定 PDB 和指定共识口袋残基生成 `pocket.lmdb`；
   - 拒绝跨构象合并口袋。
5. 编写项目运行封装：对 8YEZ、8ZU3、8YFC、9VMX 分别独立运行，保存命令、日志、版本、输入哈希和输出。
6. 将 `ranked_compounds.txt` 规范化为 CSV，至少包含 `pdb_id`、`consensus_id`、`tier`、`rank`、`smiles`、`drugclip_score`。
7. 添加小型测试和 README 使用说明；验证后再提交并推送 GitHub。

## 关键科学边界

- 四个结构必须独立运行，禁止把不同构象的口袋做交集。
- 当前项目配置是：8YEZ、8ZU3、8YFC、9VMX；8ZU8 仅参考。
- 当前仓库的等级名按队内最新口径为：T2 = 三工具共识、T1 = 两工具共识。
- DrugCLIP 分数只用于虚拟筛选排序，不等于真实结合、疗效或安全性证据；命中分子必须再经对接、MD/MM-GBSA、ADME/毒性和实验验证。

## 恢复工作时首先检查

```powershell
Set-Location 'D:\丘山\F_大二上学期\R_生物竞赛\piezo1-vs'
Get-Content .\docs\DRUGCLIP_HANDOFF_2026-09-05.md
Get-ChildItem .\tools\drugclip\official -Force
git status --short --branch
```

