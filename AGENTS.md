Chinese character involved in this project, mind encoding!


使用vsppytools python环境,位于C:\Users\zyx20\anaconda3\envs\vsppytools
openvsp 位于D:\3D\Projects\OpenVSP-3.51.3-win64


我需要你帮助我开发使用OpenVSP的分布式电推进飞机设计框架

开发规则：
1. DRY 原则 (Don't Repeat Yourself)
拒绝知识的重复。系统中的每一个功能点、算法或配置，都应有且仅有一个权威定义。禁止在多个地方手动同步相同的逻辑.

平衡点： 避免为了 DRY 而引入过度复杂的泛型或多层继承。如果消除重复会导致代码可读性急剧下降，请优先选择代码的清晰度，并辅助以显式注释.

2. 单一来源(Single Source of Truth)
常量、魔术字符串, 数据库 Schema 必须定义在集中配置文件中

3. Fail-Fast 机制暴露错误
禁止过度防御性编程，不使用 config.get('max_workers', 4)的默认参数，必须让潜在的错误直接通过报错暴露出来