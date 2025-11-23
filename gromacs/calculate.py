import re
import os
import sys
import math
import gzip
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
from numpy import linalg as LA
from distutils.spawn import find_executable

# ==============================================================================
#                      用户参数 (USER PARAMETERS)
# ==============================================================================

# --- 1. 输入文件设置 ---
# GROMACS 文件路径：包含 eigenfreq.xvg, eigenvec.trr, topol.top 的文件夹
# !!! 这是您最需要修改的地方 !!!
gromacs_path = "./"

# 分子名称 (用于图表标题和图例)
molecule_name = "Phospholipid"

# --- 2. 光谱图参数 ---
# 图谱起始和结束频率 (单位: cm⁻¹)
plot_start_freq = 500
plot_stop_freq = 3900

# 绘图步长
plot_step_size = 4

# 峰宽 (半峰全宽 FWHM, 单位: cm⁻¹)
# 凝聚相体系的合理值为 10-30 之间
gamma_fwhm = 24

# 分子是否为线性？(True 或 False)
# 这会影响扣除平动和转动自由度的数量 (非线性分子扣除6个，线性分子扣除5个)
is_linear = False

# --- 3. (可选) 与量子化学(QM)结果对比 ---
# 是否进行与 QM 结果的对比？ (True 或 False)
compare_qm = False

# QM 计算的日志文件路径 (例如 Gaussian 的 .log 或 .log.gz 文件)
# 仅在 compare_qm = True 时需要修改此路径
qm_log_file = "./path/to/your/qm_log_file.log.gz"

# --- 4. 输出文件名 ---
output_filename = "infrared_spectrum_1.png"


# ==============================================================================
#                      核心功能函数 (CORE FUNCTIONS)
# ==============================================================================

def find_gmx():
    """查找 GROMACS 执行文件"""
    gmx = None
    for mpi in ["_mpi", ""]:
        for double in ["_d", ""]:
            gmx = find_executable("gmx" + mpi + double)
            if gmx:
                return gmx
    if not gmx:
        sys.exit("错误: 未找到 GROMACS 执行文件！请确保 GROMACS 已安装并位于您的 PATH 中。")
    return None

def extract_gmx_data(path):
    """从 GROMACS 文件中提取数据"""
    print("正在从 GROMACS 文件中提取数据...")
    
    # --- 1. 读取 eigenfreq.xvg ---
    eigenfreq_file = os.path.join(path, "short_IR_run.xvg")
    if not os.path.exists(eigenfreq_file):
        sys.exit(f"错误: 文件 '{eigenfreq_file}' 不存在！请检查 gromacs_path 设置。")
    
    eigenfrequencies = []
    with open(eigenfreq_file, "r") as f:
        for line in f:
            if not line.startswith("@") and not line.startswith("#"):
                split_line = re.split(r"\s+", line.strip())
                if len(split_line) >= 2:
                    eigenfrequencies.append(float(split_line[1]))

    # --- 2. 读取 eigenvec.trr ---
    eigenvec_trr_file = os.path.join(path, "short_IR_run.trr")
    if not os.path.exists(eigenvec_trr_file):
        sys.exit(f"错误: 文件 '{eigenvec_trr_file}' 不存在！")

    gmx_executable = find_gmx()
    eigenvec_txt_file = os.path.join(path, "eigenvec.txt")
    
    dump_command = f"{gmx_executable} dump -f {eigenvec_trr_file} > {eigenvec_txt_file}"
    try:
        sp.run(dump_command, shell=True, check=True, capture_output=True)
    except sp.CalledProcessError as e:
        sys.exit(f"错误: 执行 'gmx dump' 命令失败。\n{e.stderr.decode()}")

    eigenvectors = []
    current_eigenvector = []
    first_vector_found = False
    vector_found = False
    with open(eigenvec_txt_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("x (") and not first_vector_found:
                first_vector_found = True
            elif line.startswith("x (") and first_vector_found and not vector_found:
                vector_found = True
            elif line.startswith("x[") and first_vector_found and vector_found:
                values = re.split(r"[{},\s]+", line)
                current_eigenvector.append([float(values[2]), float(values[3]), float(values[4])])
            elif vector_found and first_vector_found:
                if current_eigenvector:
                    eigenvectors.append(current_eigenvector)
                current_eigenvector = []
                vector_found = False
    
    os.remove(eigenvec_txt_file)
    if current_eigenvector:
        eigenvectors.append(current_eigenvector)
    eigenvectors = np.asarray(eigenvectors, dtype=object)
    
    # --- 3. 读取 topol.top ---
    topol_file = os.path.join(path, "topol.top")
    if not os.path.exists(topol_file):
        sys.exit(f"错误: 文件 '{topol_file}' 不存在！")

    masses = []
    charge_mass_factors = []
    in_atoms_section = False
    with open(topol_file, "r") as f:
        for line in f:
            if line.strip().startswith("[ atoms ]"):
                in_atoms_section = True
                continue
            if in_atoms_section and line.strip().startswith("["):
                in_atoms_section = False
                continue
            if in_atoms_section and not line.strip().startswith(";") and line.strip():
                words = line.strip().split()
                if len(words) >= 8:
                    masses.append(float(words[7]))
                    charge_mass_factors.append(float(words[6]) / math.sqrt(float(words[7])))

    charge_mass_factors = np.asarray(charge_mass_factors)

    # --- 4. 校正特征向量 ---
    corrected_eigenvectors = []
    for eigenvector_raw in eigenvectors:
        eigenvector = np.array(eigenvector_raw, dtype=float)
        for j in range(np.size(eigenvector, 0)):
            eigenvector[j, :] = eigenvector[j, :] * math.sqrt(masses[j])
        corrected_eigenvectors.append(eigenvector / LA.norm(eigenvector))
    
    print("GROMACS 数据提取完成。")
    return eigenfrequencies, np.array(corrected_eigenvectors), charge_mass_factors


def calculate_intensities(eigenvectors, charge_mass_factors):
    """根据公式 (13) 计算红外强度"""
    print("正在计算红外强度...")
    intensities = []
    for n in range(np.size(eigenvectors, 0)):
        intensity = 0
        for k in range(3):  # 遍历 x, y, z 三个维度
            intensity += (charge_mass_factors.dot(eigenvectors[n][:, k]))**2
        intensities.append(intensity)
    print("强度计算完成。")
    return intensities

def get_QM_IR(log_file):
    """从 QM 日志文件中提取频率和强度"""
    print(f"正在从 QM 文件 '{log_file}' 中提取数据...")
    if not os.path.exists(log_file):
        print(f"警告: QM 文件 '{log_file}' 不存在！将跳过对比。")
        return None, None
        
    frequencies = []
    intensities = []
    open_func = gzip.open if log_file.endswith(".gz") else open
    
    try:
        with open_func(log_file, 'rt') as f:
            lines = f.readlines()
        
        for line in lines:
            if "Frequencies" in line:
                words = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                frequencies.extend([float(word) for word in words])
            elif "IR Inten" in line:
                words = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                intensities.extend([float(word) for word in words])
    except Exception as e:
        print(f"警告: 读取 QM 文件时出错: {e}。将跳过对比。")
        return None, None

    if frequencies and intensities:
        print("QM 数据提取完成。")
        return frequencies, intensities
    else:
        print("警告: 未在 QM 文件中找到频率和强度数据。将跳过对比。")
        return None, None

def create_cauchy_distribution(v, v0, gamma, I):
    """为单个振动模式创建柯西/洛伦兹分布"""
    return I * (1 / math.pi) * (((1/2) * gamma) / ((v - v0)**2 + ((1/2) * gamma)**2))

def plot_spectra(gmx_freqs, gmx_intens, **kwargs):
    """绘制并保存红外光谱图"""
    print("正在生成光谱图...")
    plt.figure(figsize=(18, 8))
    
    v = np.linspace(plot_start_freq, plot_stop_freq, 
                    int((plot_stop_freq - plot_start_freq) / plot_step_size) + 1)
    
    # --- 绘制 GROMACS 谱图 ---
    spectrum_total_gmx = np.zeros(len(v))
    non_vibrational = 5 if is_linear else 6
    
    if len(gmx_freqs) > non_vibrational:
        for i in range(non_vibrational, len(gmx_freqs)):
            spectrum = create_cauchy_distribution(v, gmx_freqs[i], gamma_fwhm, gmx_intens[i])
            spectrum_total_gmx += spectrum
        # 归一化处理
        if np.amax(spectrum_total_gmx) > 0:
            spectrum_total_gmx /= np.amax(spectrum_total_gmx)
    plt.plot(v, spectrum_total_gmx, label=f'GROMACS-based spectrum for {molecule_name}')

    # --- (可选) 绘制 QM 谱图 ---
    if kwargs.get('qm_freqs') and kwargs.get('qm_intens'):
        qm_freqs = kwargs['qm_freqs']
        qm_intens = kwargs['qm_intens']
        spectrum_total_qm = np.zeros(len(v))
        for i, freq in enumerate(qm_freqs):
            spectrum = create_cauchy_distribution(v, freq, gamma_fwhm, qm_intens[i])
            spectrum_total_qm += spectrum
        
        # 归一化处理
        if np.amax(spectrum_total_qm) > 0:
            spectrum_total_qm /= np.amax(spectrum_total_qm)
        plt.plot(v, spectrum_total_qm, label=f'QM-based spectrum for {molecule_name}', linestyle='--')
        
    plt.legend(loc='upper right')
    plt.ylabel('Normalized IR Intensity')
    plt.xlabel('Frequency (cm⁻¹)')
    plt.title(f'Calculated Infrared Spectrum for {molecule_name}')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.xlim(plot_start_freq, plot_stop_freq)
    plt.gca().invert_xaxis() # 红外光谱通常x轴反向
    
    # --- 保存图像 ---
    plt.savefig(output_filename, dpi=300)
    print(f"光谱图已成功保存为 '{output_filename}'")
    plt.show()

# ==============================================================================
#                      主执行模块 (MAIN EXECUTION)
# ==============================================================================

if __name__ == "__main__":
    # 从 GROMACS 获取数据并计算强度
    gmx_frequencies, gmx_eigenvectors, gmx_cm_factors = extract_gmx_data(gromacs_path)
    gmx_intensities = calculate_intensities(gmx_eigenvectors, gmx_cm_factors)
    
    plot_kwargs = {}
    
    # 如果需要，获取 QM 数据
    if compare_qm:
        qm_frequencies, qm_intensities = get_QM_IR(qm_log_file)
        if qm_frequencies and qm_intensities:
            plot_kwargs['qm_freqs'] = qm_frequencies
            plot_kwargs['qm_intens'] = qm_intensities
            
    # 绘制光谱
    plot_spectra(gmx_frequencies, gmx_intensities, **plot_kwargs)