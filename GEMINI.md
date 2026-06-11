Chinese character involved in this project, mind encoding!
使用myml python环境,位于D:\Software\anaconda\envs\myml
我需要你帮助我开发使用OpenVSP的分布式电推进飞机设计框架
当前终端为powershell 5.1, 不支持 && 以及||操作符

开发规则：
1. DRY 原则 (Don't Repeat Yourself)
拒绝知识的重复。系统中的每一个功能点、算法或配置，都应有且仅有一个权威定义。禁止在多个地方手动同步相同的逻辑.

平衡点： 避免为了 DRY 而引入过度复杂的泛型或多层继承。如果消除重复会导致代码可读性急剧下降，请优先选择代码的清晰度，并辅助以显式注释.

2. 单一来源(Single Source of Truth)
常量、魔术字符串, 数据库 Schema 必须定义在集中配置文件中

3. Fail-Fast 机制暴露错误
禁止过度防御性编程，不使用 config.get('max_workers', 4)的默认参数，必须让潜在的错误直接通过报错暴露出来

Draw.io画图规范:
绘制神经网络时，输入放在左侧，输出放在右侧，横向排列。神经网络用圆角矩形表示，矩形内部用箭头连接的横向排列的圆角矩形代表每一层的配置，如激活函数，神经元数量，层类型，属性等。不写concat和reshape等格式操作。

输入和输出均使用圆角矩形，格式例如"noise, 10D", "Condition, [alpha, Re, thickness, Cl]"，注明维度和属性。连线均使用折线

控制图片长宽比适合在屏幕上阅读。绘制后导出为png,draw.io.exe位置在D:\Software\draw.io-29.6.1-windows, 导出时选择高清晰度. 生成后读取图片文件确认效果,不要删除drawio文件


mistakes:
1. 当在CLI执行 Python 代码块时，必须遵循以下原则：不要在 f-string 中使用引号嵌套；禁止在 {} 内出现反斜杠；严格引号分层。禁止在 f-string 的大括号内使用任何引号或反斜杠。如果需要打印字典内容，使用 print('text', dict['key']) 这种多参数形式，不要在字符串内部嵌套转义引号

2. 不要尝试执行单行复杂命令，换行语法容易出错。2行以上命令必须写出临时代码再执行