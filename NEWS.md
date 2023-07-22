### 更新日志：

写于2021.1.26

#### 2.1.1
limma如果表达矩阵有重复行名，则第一列为ID，随之改了get_deg和draw_volcano

ggplot2更新，手动指定颜色必须加values=，已加

ifelse更新，允许矩阵数据，生成结果也为矩阵，为此改了表达矩阵箱线图、KM图系列

表达矩阵与分组信息的匹配

#### 2.1.2
mRNA与lncRNA都有时矩阵去重

#### 2.1.5

GEO下载的表达矩阵有异常值时不再报错，而是只给个warning啦

#### 2.1.6

解决mac tab键编码方式报错问题

#### 2.1.7

draw_heatmap 添加参数scale，F即为原矩阵画图，不scale.

#### 2.1.8

draw_heatmap 添加参数main，标题

#### 2.1.9

surv_cox 输出结果小数点位数不限制（因为有些p值太小，会变成零）

#### 2.2.0

检查了exp_boxplot之前的所有函数及其帮助文档，quickenrich和make_tcga_group不再报warning

#### 2.2.1

检查了exp_boxplot之后的所有函数及其帮助文档，调整了脚本里的函数顺序，加上编号以后方便查找。

去掉了加载R包时的提示信息，整理了所有的seealso

我要整理整理投放到cran咯，用CMD check检查发现：

帮助文档里的示例逻辑值必须写TRUE，FALSE ,不能写T和F

示例需要保证运行正确，示例代码里library的包也需要写进依赖包。

写函数的代码里的require用requireNamespace()代替。

data用use_data生成，不自己保存，避免错误的编码方式

#### 2.2.2 

geo_download支持指定下载读取的目录，支持去除重复的ch1列

添加draw_tsne、draw_KM。ggsurvplot传参问题在：https://github.com/kassambara/survminer/issues/342

删掉了split_list函数，cran不让使用全局变量

match_exp_cl以列表形式输出，也取消了全局变量。需要赋值然后取子集。

修改了默认配色。

帮助文档的文件输出路径改为临时路径。消除全部note啦！！！

#### 2.2.3 

合并了get_deg和multi_deg，用同样的代码完成二组和多组的差异分析。

cg的组织方式简化一些

#### 2.2.4和2.2.5

帮助文档所有destdir设置为临时路径，大于5s的示例加上donttest

#### 2.2.6

cran 人工返回修改意见：

1.标题 缩写

2.包的功能详细说明

3.使用的方法参考资料，doi或者链接，可以有标题

4.TF换成TRUE FALSE

5.message代替print

6.示例里面删掉rm(list= ls())

### 2.2.8

1.exp_surv添加了cut.point参数，默认值false，即以中位数为截断值画图。修改了配色。

2.解决了sur_cox的NA报错问题

3.解决了quick_enrich画图横坐标消失的问题。

4.新增函数risk_plot，画风险因子三图联动

5.修改了suggest的包安装提示

6.内置数据exprSet_hub1改为了logcpm，不再是原来的count

### 2.2.9

1.删除heat_id 和gene_number两个参数，用my_genes代替

2.随ggplot2更新了 linewidth参数和 after_stat(p.format)

3.作图函数，加上...过渡参数

4.cor.one 和cor.full添加过滤0值参数

5.更新trans_exp函数

### 2.3.1

1.差异分析、富集分析、id转换的函数支持人，小鼠，大鼠三个物种

2.修复get_deg_all不显示差异基因、logFC阈值设置无效的bug

3.trans_exp_new改到支持数据框
