# average_spectra.py

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def average_spectra(file_list, output_data_file, output_image_file):
    """
    Reads multiple spectrum data files, averages them, and saves/plots the result.
    """
    if len(file_list) < 2:
        print("Error: Please provide at least two files to average.")
        sys.exit(1)

    print(f"Averaging {len(file_list)} spectrum files...")

    # --- 1. Load the first file to get the reference wavenumbers ---
    try:
        print(f"  - Reading reference file: {file_list[0]}")
        reference_wavenumbers, first_intensity = np.loadtxt(file_list[0], unpack=True)
    except Exception as e:
        print(f"Error reading the first file {file_list[0]}: {e}")
        sys.exit(1)

    # List to hold all intensity arrays
    intensities = [first_intensity]

    # --- 2. Loop through the rest of the files ---
    for filename in file_list[1:]:
        try:
            print(f"  - Reading file: {filename}")
            wavenumbers, intensity = np.loadtxt(filename, unpack=True)
            
            # --- CRITICAL CHECK: Ensure wavenumbers are consistent ---
            if not np.allclose(reference_wavenumbers, wavenumbers):
                print(f"Error: Wavenumber axes in {filename} do not match the reference file.")
                sys.exit(1)
            
            intensities.append(intensity)

        except Exception as e:
            print(f"Error reading file {filename}: {e}")
            sys.exit(1)
            
    # --- 3. Perform the averaging ---
    print("\nAll files read successfully. Performing average...")
    # Convert list of arrays into a 2D numpy array
    intensity_matrix = np.array(intensities)
    # Calculate the mean along the columns (axis=0)
    averaged_intensity = np.mean(intensity_matrix, axis=0)

    # --- 4. Save the averaged data ---
    output_data = np.vstack((reference_wavenumbers, averaged_intensity)).T
    header = "Column 1: Wavenumber (cm^-1)\nColumn 2: Averaged Intensity (a.u.)"
    np.savetxt(output_data_file, output_data, header=header, fmt='%18.8f')
    print(f"Averaged spectrum data saved to {output_data_file}")

    # --- 5. Plot the final, averaged spectrum ---
    plt.figure(figsize=(12, 7))
    plt.plot(reference_wavenumbers, averaged_intensity, color='red', linewidth=1.2)
    plt.title('POPC IR Spectrum', fontsize=16)
    plt.xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    plt.ylabel('Intensity (a.u.)', fontsize=12)
    plt.xlim(0, 2000)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(output_image_file, dpi=300, bbox_inches='tight')
    print(f"Final plot saved to {output_image_file}")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Average multiple GROMACS IR spectrum data files.")
    parser.add_argument("input_files", nargs='+', help="A list of all .dat files to be averaged.")
    parser.add_argument("--output_data", default="averaged_spectrum.dat", help="Path to save the final averaged data file.")
    parser.add_argument("--output_image", default="averaged_spectrum.png", help="Path to save the final plot image.")
    args = parser.parse_args()
    
    average_spectra(args.input_files, args.output_data, args.output_image)