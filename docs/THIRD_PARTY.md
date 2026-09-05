# 第三方工具与服务

本仓库只提供编排、解析和共识逻辑，不重新分发下列软件或权重。

| 名称 | 用途 | 官方地址 | 备注 |
|---|---|---|---|
| P2Rank | 结合位点预测 | https://github.com/rdk/p2rank | 使用前核对官方许可证与引用 |
| DoGSite3 / ProteinsPlus | 远程口袋预测 | https://proteins.plus/ | 受远程服务可用性和使用条款约束 |
| fpocket 4.2.3 (`4bb0d844…`) | 几何口袋预测 | https://github.com/Discngine/fpocket | 已用`scripts/bootstrap_fpocket_wsl.ps1`在隔离WSL2环境从官方源码编译；MIT许可证 |
| DrugCLIP | 口袋—小分子检索 | https://github.com/bowen-gao/DrugCLIP | 固定提交见配置；源码Apache-2.0，权重/输出CC-BY-NC-4.0；checkpoint仅本地保存，不进入Git |
| GNINA 1.3.3 | 分子对接 | https://github.com/gnina/gnina | 固定官方CUDA 12.8二进制、大小与SHA256；WSL需官方CUDA运行库和cuDNN 9 |
| SwissADME | 理化/药物相似性早筛 | https://www.swissadme.ch/ | 预测不是实验结论 |
| ADMETlab 3.0 | ADME与药化性质预测 | https://admetlab3.scbdd.com/ | 保存API参数、原始响应和不确定性 |
| ProTox 3.0 | 计算毒性筛查 | https://tox.charite.de/protox3/ | 免费Web/API不等于平台源码开源 |
| GROMACS / gmx_MMPBSA | MD与相对结合能分析 | https://www.gromacs.org/ / https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA | 仅在专家审核和算力就绪后使用 |

论文、PPT和最终提交物必须按照各项目官方要求引用并遵守其许可。
