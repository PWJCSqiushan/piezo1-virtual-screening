# 实施计划

## Gate 1：结构闭环

验收标准：

1. 8YEZ、8ZU3、8YFC、9VMX的PDB/mmCIF与RCSB API元数据均已下载；8ZU8只作参考。
2. 每个文件有 SHA256。
3. `structure_manifest.csv` 明确分辨率、组成、UniProt、链、建模/缺失残基比例和 MDFIC 情况。
4. 原始文件保持不变，任何预处理输出进入 `data/processed/`。
5. 对突变结构记录突变位点是否真实出现在原子坐标中；条目标题或序列标识不能替代坐标可见性检查。

## Gate 2：口袋闭环

验收标准：

1. 每个正式PDB都必须分别产生DoGSite3、fpocket、P2Rank原始输出；若fpocket被环境阻塞，只能报告双工具预检查，不能生成最终Tier。
2. 统一口袋字段：结构、工具、原口袋编号、中心、体积/分数、链、残基集合、原始文件。
3. 三聚体链和 UniProt 残基编号映射经过人工抽查。
4. 共识算法按空间中心和残基集合匹配，不只按口袋编号或单个残基计数；输出必须保留`support_count=2/3`。
5. 禁止跨PDB/跨构象求共同口袋。不同结构只能在各自完成正式口袋后进行下游结果对照。
6. 含MDFIC的结构须将结果标记为`piezo1`、`piezo1_mdfic_interface`、`mdfic`或`other_or_unknown`，由结构组决定哪些范围进入下游。

## Gate 3：AI 闭环

验收标准：

1. DrugCLIP 官方 DUD-E 或 retrieval 示例可重复运行。
2. 1 个 PIEZO1 pocket 与 100-1000 个小分子可生成 LMDB 并输出分数、排名和运行日志。
3. 若两天内无法复现，DrugCLIP 降为实验分支，主流程以 GNINA 基线继续。

截至 2026-09-06：官方 checkpoint 已固定并通过 SHA256 校验，WSL/CUDA/Uni-Core 环境已跑通；项目自己的 8YEZ/C001 与 C002 单口袋、5 分子端到端冒烟测试成功。Gate 3 尚未完全关闭，因为正式 100–1000 分子库仍需团队确定；上游 retrieval 的 4.68 GB 演示分子库不阻塞项目侧输入与检索功能。

## Gate 4：Docking 闭环

验收标准：

1. 1 个构象、1 个口袋、20-100 个分子完成 GNINA。
2. docking box、模型版本、随机种子、失败记录和输出 SDF 可追溯。
3. 解析 Vina/CNN 分数并生成可人工检查的前十名姿势。

## Gate 5：扩展与验证

只有Gate 1-4完成后，才扩大化合物库或并行处理更多正式口袋。MD/gmx_MMPBSA只对最终1-3个候选执行，并以服务器/HPC与膜蛋白参数审核为启动条件。
