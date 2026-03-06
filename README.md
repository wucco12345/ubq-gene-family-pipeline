# 基因家族鉴定流程

这是一个基于蛋白组的家族鉴定流程框架。

当前已经稳定完成的家族包括：

- UBQ（Ubiquitin，PF00240）
- EF1A（eukaryotic elongation factor 1-alpha，PF00009）
- ACT（Actin，PF00022）
- PRK（Phosphoribulokinase，PF00485）

---

## 一、当前统一原则

本流程当前遵循以下原则：

1. Pfam 用于召回候选
   - Pfam 命中只用于召回可能相关的候选蛋白，不直接等同于最终家族成员。

2. BLAST 用于判定 membership
   - 最终是否属于家族，主要由 BLASTp 对 UniProt Swiss-Prot reviewed 的注释证据决定。

3. 长度用于 high_conf / rescue / excluded 分层
   - 长度用于结果分层和 QC，不再作为 BLAST 前的硬门槛。

4. 家族特异规则允许存在
   - 不同家族可根据结构域组合、repeat 数、注释特征等设置专属规则。

---

## 二、目录结构

- config/
  - 家族配置文件
- scripts/
  - 运行脚本、后处理脚本、共用函数
- inputs/
  - 输入蛋白 FASTA
- results/
  - 各家族输出结果
- db/
  - 本地数据库（如 Swiss-Prot BLAST 数据库）

---

## 三、统一流程

每个家族的一般流程如下：

1. 清洗输入文件
2. 根据 Pfam 召回候选
3. 提取候选 FASTA
4. 对候选做 BLASTp（Swiss-Prot reviewed）
5. 用家族特异规则判定 membership
6. 输出标准化结果

---

## 四、统一输出规范

每个家族的结果放在：

results/<FAMILY>/

统一输出包括：

- <FAMILY>_high_conf.id
- <FAMILY>_high_conf.fa
- <FAMILY>_rescue.id
- <FAMILY>_rescue.fa
- <FAMILY>_excluded.tsv
- <FAMILY>_nohit.id
- <FAMILY>_qc.txt

如果某个家族还保留专属文件，也允许继续保留。

---

## 五、四个已稳定家族的说明

### 1）UBQ
- Pfam：PF00240
- 先召回候选
- 对 PF00240 命中区间做 merged repeat count
- repeat >= 2 作为 polyubiquitin 高可信成员
- repeat = 1 再做 BLAST 判定
- 排除 SUMO / NEDD8 / RUB / ubiquitin-like receptor 等非 UBQ 蛋白

### 2）EF1A
- Pfam：PF00009
- 对全部候选做 BLAST
- 最终 membership 由 BLAST 注释和最小比对质量共同决定
- 保留真正 eEF1A
- 排除 EF-Tu / TUFM / Der / QQT2 / 其它非 EF1A GTPase

### 3）ACT
- Pfam：PF00022
- 对全部候选做 BLAST
- 最终只保留 true ACT
- ARP（actin-related proteins）不再作为 ACT 成员保留，而是进入 excluded

### 4）PRK
- Pfam：PF00485
- 对全部候选做 BLAST
- 最终 membership 由 phosphoribulokinase 注释决定
- 排除 uridine kinase / uridine-cytidine kinase / glycerate kinase / TTM-like 等非 PRK 蛋白

---

## 六、运行方式

运行单个家族：

bash scripts/run_family.sh ACT
bash scripts/run_family.sh UBQ
bash scripts/run_family.sh EF1A
bash scripts/run_family.sh PRK

运行全部已稳定家族：

bash scripts/run_all_families.sh

---

## 七、当前项目状态

当前已经完成并稳定下来的家族：

- UBQ
- EF1A
- ACT
- PRK

这 4 个家族已经按照第二轮原则完成整理：

- Pfam 负责 recall
- BLAST 负责 membership
- 长度 / 结构负责 high_conf、rescue 和 QC 分层

---

## 八、后续扩展方式

如果要新增一个家族，通常新增：

- config/<FAMILY>.yml
- scripts/run_<FAMILY>.sh
- scripts/post_<FAMILY>_standardize.v2.sh

然后把家族名加入：

- scripts/run_all_families.sh

推荐扩展模式：

1. Pfam 召回候选
2. 如果 membership 容易混淆，则对全部候选做 BLAST
3. 应用家族特异规则
4. 输出统一标准化结果
