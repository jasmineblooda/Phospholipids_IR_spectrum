import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt

# --- 请在这里配置您的文件和原子 ---
TPR_FILE = 'test_distance.tpr'              # <-- 修改这里
XTC_FILE = 'test_distance.xtc'              # <-- 修改这里
ATOM_1_NUMBER = 21  # 第一个原子的编号 (例如P原子)
ATOM_2_NUMBER = 34  # 第二个原子的编号 (例如C22原子)
OUTPUT_PLOT_FILE = 'distance_plot_final.png' # 输出图片文件名
# ------------------------------------

print(f"Loading simulation from {TPR_FILE} and {XTC_FILE}...")
# 加载宇宙（拓扑+轨迹）
u = mda.Universe(TPR_FILE, XTC_FILE)

# 使用原子编号精确选择原子
atom1_selection = f"bynum {ATOM_1_NUMBER}"
atom2_selection = f"bynum {ATOM_2_NUMBER}"

atom1 = u.select_atoms(atom1_selection)
atom2 = u.select_atoms(atom2_selection)

# 检查是否成功选中原子
if len(atom1) == 0 or len(atom2) == 0:
    print(f"Error: Could not select atoms {ATOM_1_NUMBER} and/or {ATOM_2_NUMBER}. Please check the numbers.")
    exit()

print(f"Successfully selected Atom 1 (ID {atom1.ids[0]}) and Atom 2 (ID {atom2.ids[0]}).")
print("Calculating distance for each frame... (This may take a moment)")

# 准备存储数据
times = []
dist_data = []

# 遍历轨迹的每一帧
for ts in u.trajectory:
    # ！！！关键修正在这里！！！
    # 直接将 atom1 和 atom2 这两个 AtomGroup 对象传递给函数，而不是它们的 .positions
    dist_value = distances.dist(atom1, atom2)[0][0] / 10.0
    
    # 记录时间和距离
    times.append(ts.time)
    dist_data.append(dist_value)

print("Calculation complete.")

# --- 绘图 ---
print("Generating plot...")
plt.figure(figsize=(10, 6))
plt.plot(times, dist_data, label=f'Distance between atom {ATOM_1_NUMBER} and {ATOM_2_NUMBER}')
plt.xlabel('Time (ps)')
plt.ylabel('Distance (nm)')
plt.title('Distance between functional group atoms over time')
plt.grid(True)
plt.legend()
plt.savefig(OUTPUT_PLOT_FILE)
print(f"Plot saved to {OUTPUT_PLOT_FILE}")
plt.show()