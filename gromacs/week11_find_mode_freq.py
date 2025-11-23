# find_mode_freq.py (Corrected Version)

import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_xvg(filename):
    """
    A robust function to parse GROMACS .xvg files, skipping malformed lines.
    """
    time = []
    projection = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('@') or line.startswith('#'):
                continue
            
            cols = line.split()
            # CRITICAL FIX: Only process lines that have at least 2 columns
            if len(cols) >= 2:
                try:
                    time.append(float(cols[0]))
                    projection.append(float(cols[1]))
                except ValueError:
                    # Skip lines that cannot be converted to numbers
                    print(f"Skipping non-numeric line: {line.strip()}")
                    
    return np.array(time), np.array(projection)


def find_peak_frequency(xvg_file):
    """
    Reads a GROMACS projection file, performs an FFT, and finds the peak frequency.
    """
    print(f"Analyzing projection file: {xvg_file}")
    
    # --- 1. Load the projection data using the robust function ---
    time, projection = read_xvg(xvg_file)
    if len(time) == 0:
        print("Error: No valid data could be read from the file.")
        return

    # --- 2. Perform Fourier Transform ---
    fft_result = np.fft.fft(projection)
    power_spectrum = np.abs(fft_result)**2
    
    # --- 3. Calculate frequency axis in wavenumbers (cm⁻¹) ---
    dt = time[1] - time[0]  # Time step in ps
    n_points = len(time)
    freqs_thz = np.fft.fftfreq(n_points, d=dt)
    wavenumbers = freqs_thz * 33.3564

    # --- 4. Find the peak frequency in the positive range ---
    positive_mask = wavenumbers > 0
    pos_wavenumbers = wavenumbers[positive_mask]
    pos_power = power_spectrum[positive_mask]
    
    peak_index = np.argmax(pos_power)
    peak_wavenumber = pos_wavenumbers[peak_index]

    print("\n" + "="*40)
    print(f"  Characteristic Frequency of this Mode: {peak_wavenumber:.2f} cm⁻¹")
    print("="*40 + "\n")
    
    # --- 5. Plot the result for verification ---
    plt.figure(figsize=(10, 6))
    plt.plot(pos_wavenumbers, pos_power, color='blue')
    plt.title(f'Power Spectrum of {xvg_file}')
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('Power (a.u.)')
    plt.axvline(peak_wavenumber, color='red', linestyle='--', label=f'Peak at {peak_wavenumber:.2f} cm⁻¹')
    plt.legend()
    plt.xlim(0, 2000)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the characteristic frequency from a GROMACS projection file.")
    parser.add_argument("input_file", help="Path to the input projection.xvg file.")
    args = parser.parse_args()
    
    find_peak_frequency(args.input_file)