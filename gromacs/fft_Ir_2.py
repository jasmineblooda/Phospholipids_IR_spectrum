import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.signal.windows import blackman

# --- Constants for Quantum Correction Factor (QCF) ---
T = 303.15  # Temperature in Kelvin, from your .mdp file
h = 6.62607015e-34  # Planck constant in J*s
k = 1.380649e-23  # Boltzmann constant in J/K

# --- Data Loading (using the robust manual method) ---
print("Reading and parsing dipcorr.xvg...")
time_list = []
acf_list = []
with open('dipcorr_OneNanoSec.xvg', 'r') as f:
    for line in f:
        if line.strip().startswith('@') or line.strip().startswith('#'):
            continue
        try:
            parts = line.strip().split()
            if len(parts) >= 2:
                time_list.append(float(parts[0]))
                acf_list.append(float(parts[1]))
        except (ValueError, IndexError):
            continue
data = np.array([time_list, acf_list]).T
print("Data loaded successfully.")

time = data[:, 0]
acf = data[:, 1]
dt = time[1] - time[0]

# --- FFT and Quantum Correction ---
print("Performing FFT and applying Quantum Correction...")
acf_windowed = acf * blackman(len(acf))
spectrum_raw = fft(acf_windowed)

# Calculate frequency axis (in Hz and cm^-1)
freq_hz = np.fft.fftfreq(len(spectrum_raw), d=dt * 1e-12)  # dt in seconds
wavenumber = freq_hz / (2.99792458e10)  # c in cm/s

# Calculate Quantum Correction Factor
qcf = np.zeros_like(freq_hz)
non_zero_freq_indices = np.where(freq_hz != 0)
hbar_w = (h / (2 * np.pi)) * (2 * np.pi * freq_hz[non_zero_freq_indices])
qcf[non_zero_freq_indices] = (hbar_w / (2 * k * T)) / np.tanh(hbar_w / (2 * k * T))
intensity = qcf * np.real(spectrum_raw) * wavenumber
print("Calculation finished.")


# ###################################################################
# ########               绘图部分修改开始                  ########
# ###################################################################
print("Generating zoomed and normalized plot...")

# 1. 定义我们感兴趣的放大区域
zoom_min_wn = 1500
zoom_max_wn = 1800

# 2. 找到这个区域内的数据并进行归一化
#   - 首先，找到所有波数在放大区域内的点的索引
zoom_indices = np.where((wavenumber > zoom_min_wn) & (wavenumber < zoom_max_wn))
#   - 然后，找到这个区域内的最大强度值
max_intensity_in_zoom = np.max(intensity[zoom_indices])
#   - 用这个最大值来归一化整个光谱的强度
if max_intensity_in_zoom != 0:
    normalized_intensity = intensity / max_intensity_in_zoom
else:
    normalized_intensity = intensity # 避免除以零

# 3. 开始绘图
plt.figure(figsize=(10, 6))
positive_freq_indices = np.where(wavenumber > 0)
plt.plot(wavenumber[positive_freq_indices], normalized_intensity[positive_freq_indices])

# 4. 设置新的坐标轴范围来实现“放大”效果
plt.gca().invert_xaxis()
plt.xlim(zoom_max_wn, zoom_min_wn)  # <-- 修改X轴范围
plt.ylim(-0.1, 1.2)                  # <-- 修改Y轴范围以适应归一化的数据

# 设置标题和标签
plt.xlabel('Wavenumber (cm⁻¹)')
plt.ylabel('Normalized Intensity (arb. units)')
plt.title('Zoomed-in Infrared Spectrum')
plt.grid(True)
plt.show()

# ###################################################################
# ########               绘图部分修改结束                  ########
# ###################################################################