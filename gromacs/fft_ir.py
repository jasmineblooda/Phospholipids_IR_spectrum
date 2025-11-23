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
# Formula: I(w) is proportional to w * tanh(hbar*w / 2kT) * FT(ACF)
# We handle the w=0 case to avoid division by zero.
qcf = np.zeros_like(freq_hz)
non_zero_freq_indices = np.where(freq_hz != 0)

hbar_w = (h / (2 * np.pi)) * (2 * np.pi * freq_hz[non_zero_freq_indices])
qcf[non_zero_freq_indices] = (hbar_w / (2 * k * T)) / np.tanh(hbar_w / (2 * k * T))

# An older, simpler QCF is sometimes used: I(w) ~ w^2 * FT(ACF), we'll use the proper one.

# Calculate final, corrected intensity
intensity = qcf * np.real(spectrum_raw) * wavenumber

print("Calculation finished.")

# --- Plotting ---
print("Generating plot...")
plt.figure(figsize=(12, 6))
positive_freq_indices = np.where(wavenumber > 0)
plt.plot(wavenumber[positive_freq_indices], intensity[positive_freq_indices])
plt.gca().invert_xaxis()
plt.xlabel('Wavenumber (cm⁻¹)')
plt.ylabel('Intensity (arbitrary units)')
plt.title('Simulated Infrared Spectrum (with Quantum Correction)')
plt.xlim(2000, 400) # We limit to 2000 since max wavenumber is ~1668
plt.grid(True)
plt.show()