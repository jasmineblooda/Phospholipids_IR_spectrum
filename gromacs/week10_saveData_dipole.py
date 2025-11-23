# calculate_ir_save_data.py

import numpy as np
import argparse

def calculate_acf(series):
    """Calculates the autocorrelation function using the Wiener-Khinchin theorem."""
    n = len(series)
    padded_series = np.append(series, np.zeros(n))
    fft_series = np.fft.fft(padded_series)
    psd = np.abs(fft_series)**2
    acf = np.fft.ifft(psd)
    return acf[:n].real

def calculate_and_save_spectrum(xvg_file, output_file):
    """
    Calculates the IR spectrum from a GROMACS dipole.xvg file and saves the
    data (wavenumber, intensity) to a text file.
    """
    print(f"Reading dipole data from {xvg_file}...")
    try:
        time, mu_x, mu_y, mu_z = np.loadtxt(
            xvg_file, 
            comments=['@', '#'], 
            usecols=(0, 1, 2, 3), 
            unpack=True
        )
    except Exception as e:
        print(f"Error reading {xvg_file}: {e}")
        return

    # --- 1. Calculate the time derivative of the dipole moment ---
    dt = time[1] - time[0]
    mu_dot_x = np.diff(mu_x) / dt
    mu_dot_y = np.diff(mu_y) / dt
    mu_dot_z = np.diff(mu_z) / dt
    
    # --- 2. Calculate the ACF of the derivative ---
    acf_x = calculate_acf(mu_dot_x)
    acf_y = calculate_acf(mu_dot_y)
    acf_z = calculate_acf(mu_dot_z)
    total_acf = acf_x + acf_y + acf_z

    # --- 3. Fourier Transform the ACF to get the spectrum ---
    window = np.hanning(len(total_acf))
    spectrum = np.abs(np.fft.fft(total_acf * window))

    # --- 4. Convert frequency to wavenumbers (cm⁻¹) ---
    freqs_thz = np.fft.fftfreq(len(spectrum), d=dt)
    wavenumbers = freqs_thz * 33.3564

    # We only need the positive frequency part
    positive_freq_indices = np.where(wavenumbers >= 0)
    wavenumbers_pos = wavenumbers[positive_freq_indices]
    spectrum_pos = spectrum[positive_freq_indices]
    
    # --- 5. Save the data to a text file ---
    # Combine the two arrays into a 2-column format
    output_data = np.vstack((wavenumbers_pos, spectrum_pos)).T
    
    # Define a header for the output file
    header = "Column 1: Wavenumber (cm^-1)\nColumn 2: Intensity (a.u.)"
    
    # Save the data
    np.savetxt(output_file, output_data, header=header, fmt='%18.8f')
    
    print(f"Spectrum data successfully saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate IR spectrum from GROMACS dipole file and save data.")
    parser.add_argument("input_file", help="Path to the input dipole.xvg file.")
    parser.add_argument("--output", default="spectrum.dat", help="Path to save the output data file (e.g., spectrum_1.dat).")
    args = parser.parse_args()
    
    calculate_and_save_spectrum(args.input_file, args.output)