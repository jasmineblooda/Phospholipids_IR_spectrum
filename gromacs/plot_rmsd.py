# plot_rmsd.py

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def plot_rmsd(xvg_file, output_image):
    """
    读取 GROMACS .xvg 文件并生成一个RMSD图。
    """
    print(f"正在读取 RMSD 数据: {xvg_file}")

    try:
        # 加载数据，跳过以 '@' 或 '#' 开头的注释行
        time, rmsd = np.loadtxt(xvg_file, comments=['@', '#'], unpack=True)
    except Exception as e:
        print(f"错误: 无法读取数据文件 '{xvg_file}'.")
        print(f"详情: {e}")
        sys.exit(1)

    print("数据加载成功，正在生成图像...")

    # --- 创建图像 ---
    fig, ax = plt.subplots(figsize=(10, 6))

    # 绘制RMSD数据，并设置颜色为蓝色
    ax.plot(time, rmsd, color='blue', linewidth=1.0, label='RMSD')

    # --- 美化图像 ---
    # 从.xvg文件的@ subtitle中获取标题
    # 您可以手动修改为您想要的任何标题
    ax.set_title('RMSD: POPC_single after lsq fit to POPC_single', fontsize=16)
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('RMSD (nm)', fontsize=12)

    # 设置坐标轴范围
    ax.set_xlim(left=0, right=max(time))
    ax.set_ylim(bottom=0)

    # 添加网格
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # --- 保存并显示图像 ---
    try:
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        print(f"图像已成功保存到: {output_image}")
    except Exception as e:
        print(f"错误: 无法保存图像到 '{output_image}'.")
        print(f"详情: {e}")

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="从 GROMACS RMSD 文件生成图像。")
    parser.add_argument("input_file", help="输入的 .xvg 文件 (例如 'rmsd_popc.xvg')")
    parser.add_argument("--output", default="RMSD_Plot_blue.png", help="输出图像的文件名。")
    
    args = parser.parse_args()
    
    plot_rmsd(args.input_file, args.output)