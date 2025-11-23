# calculate_ir.py (Corrected for 5-column input)

import numpy as np
import matplotlib.pyplot as plt
import argparse

def calculate_acf(series):
    """Calculates the autocorrelation function using the Wiener-Khinchin theorem."""
    n = len(series)
    # Pad the series to avoid circular correlation artifacts
    padded_series = np.append(series, np.zeros(n))
    # Fourier transform
    fft_series = np.fft.fft(padded_series)
    # Power spectral density
    psd = np.abs(fft_series)**2
    # Inverse Fourier transform gives the ACF
    acf = np.fft.ifft(psd)
    # Return the real part of the first half (positive lags)
    return acf[:n].real

def calculate_ir_spectrum(xvg_file, output_image='ir_spectrum.png'):
    """
    Calculates and plots the IR spectrum from a GROMACS dipole.xvg file.
    The method follows the formula: IR ~ FT[<dM/dt(t) * dM/dt(0)>]
    """
    print(f"Reading dipole data from {xvg_file}...")
    # Load the dipole moment data, skipping comment lines
    try:
        # --- THIS IS THE CORRECTED LINE ---
        # We specify usecols to only read the first four columns (0, 1, 2, 3)
        time, mu_x, mu_y, mu_z = np.loadtxt(
            xvg_file, 
            comments=['@', '#'], 
            usecols=(0, 1, 2, 3), 
            unpack=True
        )
    except Exception as e:
        print(f"Error reading {xvg_file}. Make sure it's a valid GROMACS xvg file.")
        print(e)
        return

    # --- 1. Calculate the time derivative of the dipole moment ---
    dt = time[1] - time[0]  # Time step in ps
    
    # Use numpy.diff for finite difference calculation
    mu_dot_x = np.diff(mu_x) / dt
    mu_dot_y = np.diff(mu_y) / dt
    mu_dot_z = np.diff(mu_z) / dt
    
    print(f"Calculated time derivatives. Time step dt = {dt:.4f} ps.")

    # --- 2. Calculate the autocorrelation function (ACF) of the derivative ---
    print("Calculating autocorrelation function (ACF)...")
    acf_x = calculate_acf(mu_dot_x)
    acf_y = calculate_acf(mu_dot_y)
    acf_z = calculate_acf(mu_dot_z)
    
    # The total ACF is the sum of the components
    total_acf = acf_x + acf_y + acf_z

    # --- 3. Fourier Transform the ACF to get the spectrum ---
    print("Performing Fourier Transform to get the spectrum...")
    # Apply a window function (e.g., Hanning) to reduce spectral leakage
    window = np.hanning(len(total_acf))
    spectrum = np.abs(np.fft.fft(total_acf * window))

    # --- 4. Convert frequency to wavenumbers (cm⁻¹) ---
    # FFT frequencies are in cycles/ps (THz)
    freqs_thz = np.fft.fftfreq(len(spectrum), d=dt)
    
    # Conversion factor: 1 THz = 33.3564 cm⁻¹
    wavenumbers = freqs_thz * 33.3564

    # We only need the positive frequency part of the spectrum
    positive_freq_indices = np.where(wavenumbers >= 0)
    wavenumbers = wavenumbers[positive_freq_indices]
    spectrum = spectrum[positive_freq_indices]
    
    print(f"Plotting spectrum and saving to {output_image}")

    # --- 5. Plot the final spectrum ---
    plt.figure(figsize=(10, 6))
    plt.plot(wavenumbers, spectrum, color='black', linewidth=1.0)
    plt.title('Calculated IR Spectrum', fontsize=16)
    plt.xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    plt.ylabel('Intensity (a.u.)', fontsize=12)
    plt.xlim(0, 2000)  # Set a typical range for mid/far-IR
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(output_image, dpi=300)
    plt.show()
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate IR spectrum from GROMACS dipole file.")
    parser.add_argument("input_file", help="Path to the input dipole.xvg file.")
    parser.add_argument("--output", default="ir_spectrum.png", help="Path to save the output image file.")
    args = parser.parse_args()
    
    calculate_ir_spectrum(args.input_file, args.output)