# 协作指南

## 推荐流程

1. 从`main`建立短分支，例如`feature/fpocket-parser`。
2. 每个分支只解决一个清晰问题。
3. 运行测试并在Pull Request中写明输入PDB、工具版本、完整命令和实际验证范围。
4. 由另一名成员复核后合并。

## 提交前检查

```bash
python -m compileall -q src scripts tests
python -m unittest discover -s tests -v
python -m json.tool config/structures.json > /dev/null
python -m json.tool config/consensus_rules.json > /dev/null
```

PowerShell可把`> /dev/null`改成`| Out-Null`。

## 严禁提交

- API key、token、密码、Cookie或账号信息。
- `C:\Users\...`、`D:\...`等个人绝对路径。
- 大型结构、运行结果、模型权重或第三方完整源码。
- 将计算预测描述为已验证药物、临床安全性或COPD疗效。
- 将不同PDB构象的口袋交集当作正式共识。

## 数据约定

- PDB ID统一大写，每个结果必须携带`pdb_id`。
- Tier名称可能调整，任何表格都必须保留`support_count`。
- 失败记录不得覆盖；重跑应产生新运行目录。
