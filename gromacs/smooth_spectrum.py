import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# =====================================================================
# ==                      配置参数 (请在这里修改)                      ==
# =====================================================================

# 输入和输出文件名
INPUT_FILE = 'week10_averaged_spectra.dat'
OUTPUT_FILE = 'week10_averaged_spectra_smoothed.dat'

# 高斯平滑的关键参数：半峰全宽 (Full Width at Half Maximum)
# 这个值越大, 光谱就越平滑。单位是 cm⁻¹。
# 建议从 15 或 20 开始尝试。
FWHM = 30

# 绘图的X轴范围
PLOT_X_MIN = 0
PLOT_X_MAX = 2000

# =====================================================================
# ==                          主程序 (无需修改)                        ==
# =====================================================================

print(f"正在读取文件: {INPUT_FILE}")
# 1. 读取数据文件 (忽略以 @ 和 # 开头的注释行)
try:
    data = np.loadtxt(INPUT_FILE, comments=['@', '#'])
    wavenumber = data[:, 0]
    intensity = data[:, 1]
except Exception as e:
    print(f"读取文件时出错: {e}")
    exit()

print("数据读取成功。")

# 2. 计算高斯滤波所需的 sigma 值
# FWHM 和 sigma 的转换公式为: FWHM = 2 * sqrt(2 * ln(2)) * sigma ≈ 2.355 * sigma
sigma_in_wavenumbers = FWHM / 2.35482

# 3. 将 sigma 从物理单位 (cm⁻¹) 转换为数据点的索引单位
# 首先计算数据点之间的平均间距
spacing = np.mean(np.diff(wavenumber))
sigma_in_indices = sigma_in_wavenumbers / spacing

print(f"设置的 FWHM = {FWHM} cm⁻¹")
print(f"计算出的数据点间距 (spacing) ≈ {spacing:.4f} cm⁻¹")
print(f"用于高斯滤波的 sigma (以数据点为单位) = {sigma_in_indices:.2f}")

# 4. 执行高斯平滑
# 这是最核心的一步，我们对强度数据进行一维高斯滤波
smoothed_intensity = gaussian_filter1d(intensity, sigma=sigma_in_indices, mode='constant')
print("高斯平滑处理完成。")


# 5. 绘制结果图进行对比
plt.style.use('default')
fig, ax = plt.subplots(figsize=(12, 7))

# 绘制原始数据
ax.plot(wavenumber, intensity, color='black', linewidth=0.7, label='Original Calculated Data', alpha=0.8)

# 绘制平滑后的数据
ax.plot(wavenumber, smoothed_intensity, color='red', linewidth=2.0, label=f'Smoothed Data (FWHM = {FWHM} cm⁻¹)')

# 设置图表样式
ax.set_title('Comparison of Original and Smoothed IR Spectrum', fontsize=16)
ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
ax.set_ylabel('Intensity (a.u.)', fontsize=12)
ax.set_xlim(PLOT_X_MIN, PLOT_X_MAX)
ax.set_ylim(bottom=0)
ax.grid(True, which='both', linestyle='--', linewidth=0.5)
ax.legend(fontsize=12)

print("正在显示对比图...")
plt.show()


# 6. 保存平滑后的数据到新文件
smoothed_data = np.vstack((wavenumber, smoothed_intensity)).T
header_text = f"Smoothed IR Spectrum. Original file: {INPUT_FILE}, FWHM: {FWHM} cm^-1\nColumn 1: Wavenumber (cm^-1)\nColumn 2: Intensity (a.u.)"
np.savetxt(OUTPUT_FILE, smoothed_data, fmt='%.6f', header=header_text)

print(f"平滑后的数据已成功保存到: {OUTPUT_FILE}")