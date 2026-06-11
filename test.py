import openvsp
import os
import math
import numpy as np
#import skopt

tess_int = 0.01

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
# 切换当前工作目录
os.chdir(output_dir)
def next_tess_value(x):
    n = max((x - 6) // 4 + 1, 0)
    return 5 + n * 4

print(np.sin(90))
prop_choice = {8:"8x6", 9:"9x6", 10:"10x7", 11:"11x7", 12:"12x8", 13:"13x8"}
for dia in prop_choice:
    print(dia)
# 或者
for dia in prop_choice.keys():
    print(dia)

# 遍历所有尺寸字符串（值）
for size in prop_choice.values():
    print(size)

min_key = min(prop_choice)  # 8
max_key = max(prop_choice)  # 13

print("最小直径：", min_key)
print("最大直径：", max_key)

#skopt.benchmarks.bench3([5.0])
for i in range(0,0):
    print("0")

testlist = np.linspace(0,10,5)
print(testlist)

gset = [0,1, None]
if gset[2]:
    print("yes")

teststring = "hello"
if teststring:
    print("it works")

sub_cfg = []
if sub_cfg:
    print("not empty")


a = np.array([1,2,3,4])
b = np.array([1,2,3,4])
print(np.sum(a*b))