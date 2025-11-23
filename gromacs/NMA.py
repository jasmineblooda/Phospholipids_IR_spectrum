#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
仅用 GROMACS 简正模式结果画红外光谱（无 QM 对比）
"""
import re, os, sys, math
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
from numpy import linalg as LA
from distutils.spawn import find_executable

# ---------- 1. 找 GROMACS ----------
def find_gmx():
    for sfx in ("", "_mpi", "_d", "_mpi_d"):
        gmx = find_executable("gmx" + sfx)
        if gmx: return gmx
    sys.exit("❌ GROMACS not found!")

# ---------- 2. 必填路径 ----------
CGENFF_PATH = "D:/Learn/Thesis/gromacsLearning/Phospholipids_Simulation/charmm-gui-5548695784/gromacs/"         
OUT_PNG     = "IR_MD_only_1.png"       # <<<<<< 输出图名

# ---------- 3. 重新运行简正模式分析 ----------
def rerun_normal_mode_analysis(path):
    """重新运行简正模式分析"""
    print("→ 重新运行简正模式分析...")
    
    gmx = find_gmx()
    original_dir = os.getcwd()
    
    try:
        os.chdir(path)
        
        # 检查必要文件
        required_files = ["confout.gro", "topol.top", "ener.edr"]
        missing_files = [f for f in required_files if not os.path.exists(f)]
        if missing_files:
            print(f"❌ 缺少必要文件: {missing_files}")
            return False
        
        print("✅ 找到所有必要文件")
        
        # 运行 nmeig 命令
        cmd = f'echo "0" | {gmx} nmeig -f ener.edr -s topol.top -first 1 -last 50 -T 300'
        print(f"执行命令: {cmd}")
        
        result = sp.run(cmd, shell=True, capture_output=True, text=True)
        
        print("标准输出:")
        print(result.stdout)
        if result.stderr:
            print("错误输出:")
            print(result.stderr)
        
        # 检查是否成功生成文件
        success = os.path.exists("eigenfreq.xvg") and os.path.exists("eigenvec.trr")
        if success:
            print("✅ 简正模式分析完成")
            # 检查频率文件内容
            with open("eigenfreq.xvg", "r") as f:
                lines = f.readlines()
                data_lines = [line for line in lines if not line.startswith(("@", "#"))]
                if data_lines:
                    first_freq = data_lines[0].split()[1]
                    print(f"第一个频率: {first_freq}")
        else:
            print("❌ 简正模式分析失败")
            
        return success
        
    except Exception as e:
        print(f"❌ 运行出错: {e}")
        return False
    finally:
        os.chdir(original_dir)

# ---------- 4. 检查频率文件 ----------
def check_frequency_file(path):
    """检查频率文件内容"""
    freq_file = os.path.join(path, "eigenfreq.xvg")
    
    if not os.path.exists(freq_file):
        print(f"❌ 频率文件不存在: {freq_file}")
        return False
    
    frequencies = []
    with open(freq_file, 'r') as f:
        for line in f:
            if line.startswith(("@", "#")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    freq = float(parts[1])
                    frequencies.append(freq)
                except ValueError:
                    continue
    
    print(f"找到 {len(frequencies)} 个频率")
    if frequencies:
        print(f"频率范围: {min(frequencies):.2f} 到 {max(frequencies):.2f} cm⁻¹")
        print(f"前5个频率: {frequencies[:5]}")
        
        # 检查是否有非零频率
        non_zero_freqs = [f for f in frequencies if f > 0]
        if non_zero_freqs:
            print(f"✅ 找到 {len(non_zero_freqs)} 个非零频率")
            return True
        else:
            print("❌ 所有频率都为0")
            return False
    else:
        print("❌ 没有找到频率数据")
        return False

# ---------- 5. 创建示例光谱（用于测试） ----------
def create_example_spectrum(out_png):
    """创建示例红外光谱"""
    print("→ 创建示例红外光谱...")
    
    # 创建一些示例频率和强度
    x = np.linspace(500, 3900, 1000)
    y = np.zeros_like(x)
    
    # 添加几个示例峰
    example_peaks = [
        (800, 50, 0.3),   # (位置, 宽度, 强度)
        (1200, 40, 0.8),
        (1650, 30, 1.0),
        (2850, 35, 0.6),
        (2950, 25, 0.9),
        (3300, 45, 0.7)
    ]
    
    for pos, width, intensity in example_peaks:
        y += intensity * np.exp(-(x - pos)**2 / (2 * width**2))
    
    y /= y.max()
    
    plt.figure(figsize=(12, 6))
    plt.plot(x, y, 'b-', linewidth=2, label='Example IR Spectrum')
    plt.xlabel('Frequency (cm⁻¹)')
    plt.ylabel('Normalized Intensity')
    plt.title('Example Infrared Spectrum\n(Real data requires proper normal mode analysis)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xlim(500, 3900)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.show()
    
    print(f"✅ 示例光谱保存为: {out_png}")

# ---------- 6. 显示使用说明 ----------
def show_instructions():
    """显示GROMACS简正模式分析的使用说明"""
    print("\n" + "="*60)
    print("GROMACS 简正模式分析指南")
    print("="*60)
    print("\n要获得正确的红外光谱，需要完成以下步骤:")
    print("\n1. 能量最小化 (gmx grompp + gmx mdrun -v -deffnm em)")
    print("2. 简正模式分析 (gmx nmeig):")
    print("   gmx nmeig -f ener.edr -s topol.top -first 7 -last 1000")
    print("3. 选择系统: 通常选择 0 (整个系统)")
    print("\n常用命令序列:")
    print("gmx grompp -f minim.mdp -c conf.gro -p topol.top -o em.tpr")
    print("gmx mdrun -v -deffnm em")
    print("gmx nmeig -f em.edr -s em.tpr -first 7 -last 1000")
    print("\n确保 eigenfreq.xvg 文件中包含非零频率值!")
    print("="*60)

# ---------- 7. 主程序 ----------
def main():
    print("→ GROMACS 红外光谱分析工具")
    print("="*50)
    
    # 检查频率文件
    if not check_frequency_file(CGENFF_PATH):
        print("\n❌ 频率文件无效或所有频率为0")
        print("\n可能的原因:")
        print("1. 简正模式分析未完成")
        print("2. 能量最小化不充分")
        print("3. 输入文件有问题")
        
        response = input("\n是否尝试重新运行简正模式分析? (y/n): ").lower()
        if response == 'y':
            if rerun_normal_mode_analysis(CGENFF_PATH):
                print("✅ 请重新运行本程序检查结果")
            else:
                print("❌ 简正模式分析失败，请手动运行")
                show_instructions()
        else:
            print("⚠️  使用示例光谱进行演示")
            create_example_spectrum(OUT_PNG)
            show_instructions()
    else:
        print("✅ 频率文件有效，可以继续处理特征向量")
        # 这里可以添加处理特征向量的代码
        print("由于频率有效，接下来可以计算红外强度...")

if __name__ == "__main__":
    main()