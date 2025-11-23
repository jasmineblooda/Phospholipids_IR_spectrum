# plot_spectrum.py

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def plot_spectrum(data_file, output_image):
    """
    Reads a .dat file with spectrum data and generates a plot.
    
    Args:
        data_file (str): Path to the input data file (e.g., 'week10_averaged_spectra.dat').
        output_image (str): Path to save the output plot image (e.g., 'final_spectrum.png').
    """
    print(f"Reading spectrum data from: {data_file}")
    
    # --- 1. Load the data from the file ---
    # The 'comments' argument tells numpy to ignore lines starting with '#'
    try:
        wavenumber, intensity = np.loadtxt(data_file, comments='#', unpack=True)
    except Exception as e:
        print(f"Error: Could not read the data file '{data_file}'.")
        print(f"Details: {e}")
        sys.exit(1)
        
    print("Data loaded successfully. Generating plot...")

    # --- 2. Create the plot ---
    plt.style.use('default') # Use a standard plot style
    fig, ax = plt.subplots(figsize=(12, 7)) # Create a figure and an axes object

    # Plot the data
    ax.plot(wavenumber, intensity, color='red', linewidth=1.2, label='Averaged Spectrum')

    # --- 3. Customize the plot for clarity ---
    # Set titles and labels
    ax.set_title('POPC IR Spectrum', fontsize=16)
    ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Intensity (a.u.)', fontsize=12)

    # Set the viewing range for the x-axis
    ax.set_xlim(0, 2000)
    # Set the y-axis to start from zero
    ax.set_ylim(bottom=0)

    # Add a grid for easier reading of peak positions
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

    # --- 4. Save and show the plot ---
    try:
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        print(f"Plot successfully saved to: {output_image}")
    except Exception as e:
        print(f"Error: Could not save the plot to '{output_image}'.")
        print(f"Details: {e}")

    plt.show()

if __name__ == "__main__":
    # Set up command-line argument parsing to make the script easy to use
    parser = argparse.ArgumentParser(description="Generate a plot from a spectrum data file.")
    parser.add_argument("input_file", help="The input .dat file containing the spectrum data (e.g., week10_averaged_spectra.dat).")
    parser.add_argument("--output", default="Final_Averaged_Spectrum.png", help="The file name for the output plot image.")
    
    args = parser.parse_args()
    
    plot_spectrum(args.input_file, args.output)