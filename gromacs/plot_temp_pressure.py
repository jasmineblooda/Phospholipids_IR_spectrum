import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def read_xvg_data(filename):
    """
    一个健壮的函数，用于读取GROMACS的.xvg文件，跳过注释。
    """
    try:
        # np.loadtxt 会自动处理以'@'或'#'开头的注释行
        time, value = np.loadtxt(filename, comments=['@', '#'], unpack=True)
        return time, value
    except Exception as e:
        print(f"读取文件时出错 {filename}: {e}")
        return None, None

def plot_side_by_side(temp_file, press_file, output_file):
    """
    读取温度和压力文件，并将它们并排绘制在一张图上。
    """
    # --- 1. 读取数据 ---
    print(f"正在读取温度数据: {temp_file}")
    time_temp, temp = read_xvg_data(temp_file)
    
    print(f"正在读取压力数据: {press_file}")
    time_press, press = read_xvg_data(press_file)

    # 检查数据是否成功读取
    if time_temp is None or time_press is None:
        print("数据加载失败，退出脚本。")
        sys.exit(1)

    print("数据加载成功，开始绘图...")

    # --- 2. 创建一个1行2列的子图画布 ---
    # figsize=(16, 6) 创建一个宽：16英寸，高：6英寸的长方形画布
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # --- 3. 绘制左图：温度 ---
    ax1.plot(time_temp, temp, color='black', linewidth=0.5)
    ax1.set_title('Temperature', fontsize=16)
    ax1.set_xlabel('Time (ps)', fontsize=12)
    ax1.set_ylabel('Temperature (K)', fontsize=12)
    ax1.set_xlim(left=0, right=max(time_temp))
    # 根据您的数据 和之前的图，设置一个合理的Y轴范围
    ax1.set_ylim(290, 310) 
    ax1.grid(True, linestyle='--', alpha=0.6)

    # --- 4. 绘制右图：压力 ---
    ax2.plot(time_press, press, color='black', linewidth=0.5)
    ax2.set_title('Pressure', fontsize=16)
    ax2.set_xlabel('Time (ps)', fontsize=12)
    ax2.set_ylabel('Pressure (bar)', fontsize=12)
    ax2.set_xlim(left=0, right=max(time_press))
    # 压力波动很大，让matplotlib自动设置Y轴范围
    ax2.grid(True, linestyle='--', alpha=0.6)

    # --- 5. 调整布局并保存/显示 ---
    plt.tight_layout() # 自动调整子图间距，防止标签重叠
    
    try:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"图像已成功保存到: {output_file}")
    except Exception as e:
        print(f"保存图像失败: {e}")

    plt.show() # 在屏幕上显示图像

if __name__ == "__main__":
    # 设置命令行参数，使脚本更易用
    parser = argparse.ArgumentParser(description="绘制GROMACS的温度和压力平衡图。")
    parser.add_argument("temp_file", help="输入的温度.xvg文件 (例如 'temperature.xvg')")
    parser.add_argument("press_file", help="输入的压力.xvg文件 (例如 'pressure.xvg')")
    parser.add_argument("--output", default="Temp_Press_Plot.png", help="输出图像的文件名 (例如 'equilibration.png')")
    
    args = parser.parse_args()
    
    plot_side_by_side(args.temp_file, args.press_file, args.output)