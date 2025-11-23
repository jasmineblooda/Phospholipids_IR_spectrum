import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.constants import c, k, h

# --- 1. 用户定义的参数 ---
ACF_FILE = 'week10_dipole.xvg'  # 您的输入文件
TEMPERATURE = 310.15                 # 模拟温度 (K)
DT = 0.01                            # 轨迹帧之间的时间间隔 (ps)
                                     # 必须与您的模拟设置完全一致！
                                     # 从您的文件名看，时间间隔是0.01 ps

# --- 2. 读取.xvg文件 (修正后的函数) ---
def read_xvg(filename):
    """一个用于解析GROMACS.xvg文件的简单函数。"""
    time, acf = [], []  # <--- 这是修正后的初始化
    with open(filename, 'r') as f:
        for line in f:
            # 跳过以'@'或'#'开头的注释行
            if not line.startswith(('@', '#')):
                cols = line.split()
                if len(cols) >= 2: # 确保行中有足够的数据
                    time.append(float(cols[0])) # <--- 修正：读取第一列
                    acf.append(float(cols[1]))  # <--- 修正：读取第二列
    return np.array(time), np.array(acf)

time_ps, acf_data = read_xvg(ACF_FILE)

# --- 3. 应用窗函数以减少伪影 ---
# 有限时间的模拟会导致“频谱泄漏”。窗函数可以平滑ACF的边缘，
# 从而产生更清晰的光谱。
window = np.hanning(len(acf_data))
acf_windowed = acf_data * window

# --- 4. 执行快速傅里叶变换 (FFT) ---
# 将时域信号转换为频域信号。
power_spectrum = np.abs(fft(acf_windowed))**2

# --- 5. 计算频率轴并应用校正 ---
# 5.1 创建频率轴 (单位: cm⁻¹)
n_points = len(time_ps)
freq_ps_inv = np.fft.fftfreq(n_points, d=DT) # 频率单位: 1/ps
wavenumber = freq_ps_inv * 1e12 / (c * 100)  # 转换为波数 cm⁻¹

# 5.2 应用量子校正因子 (QCF)
# 经典MD模拟不满足量子统计，导致高频峰强度偏低。
# 这个谐振子QCF有助于校正这个问题。
# 避免在 omega = 0 时除以零
omega = 2 * np.pi * freq_ps_inv * 1e12 # 角频率 (rad/s)
qcf = np.ones_like(omega)
non_zero_mask = omega!= 0
beta_hbar_omega = (h / (2 * np.pi)) * omega[non_zero_mask] / (k * TEMPERATURE)
qcf[non_zero_mask] = beta_hbar_omega / (1 - np.exp(-beta_hbar_omega))

# 最终的红外吸收强度 I(ω) ∝ ω² * PowerSpectrum * QCF
absorption = power_spectrum * qcf * (omega**2)

# --- 6. 绘制最终的光谱图 ---
# 通常只关心正频率部分
positive_freq_mask = wavenumber > 0

plt.figure(figsize=(12, 6))
plt.plot(wavenumber[positive_freq_mask], absorption[positive_freq_mask])
plt.title('Calculated IR Spectrum')
plt.xlabel('Wavenumber (cm⁻¹)')
plt.ylabel('Intensity (arbitrary units)')
plt.xlim(0, 1800) # 设置典型的红外光谱范围
plt.grid(True)
plt.show()