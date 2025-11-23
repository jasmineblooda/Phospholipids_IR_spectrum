# modify_trr.py (Corrected Version)

import MDAnalysis as mda
import numpy as np
import argparse

def modify_trajectory(tpr_file, traj_file, output_file):
    """
    Reads a GROMACS trajectory, modifies velocities by multiplying with charges,
    and writes to a new trajectory file.

    Args:
        tpr_file (str): Path to the input .tpr file.
        traj_file (str): Path to the input .trr trajectory file.
        output_file (str): Path for the output modified .trr file.
    """
    print(f"Loading topology from: {tpr_file}")
    print(f"Loading trajectory from: {traj_file}")
    
    # Load the universe (topology from TPR, trajectory from TRR)
    try:
        u = mda.Universe(tpr_file, traj_file)
    except Exception as e:
        print(f"Error loading files: {e}")
        return

    n_atoms = len(u.atoms)
    n_frames = len(u.trajectory)
    
    print(f"Found {n_atoms} atoms and {n_frames} frames.")

    # Get charges and reshape for broadcasting with velocities
    # The shape becomes (n_atoms, 1) to multiply with (n_atoms, 3) velocity array
    charges = u.atoms.charges.reshape(-1, 1)

    # Prepare the output file writer using the stable, general mda.Writer
    # It automatically detects the file type from the '.trr' extension
    with mda.Writer(output_file, n_atoms=n_atoms) as writer:
        print(f"Writing modified trajectory to: {output_file}")
        
        # Loop through each frame in the trajectory
        for i, ts in enumerate(u.trajectory):
            # Get the original velocities for the current frame
            # We must check if the frame actually contains velocity information
            if not ts.has_velocities:
                print(f"Warning: Frame {i+1} does not contain velocities. Skipping.")
                continue

            velocities = ts.velocities
            
            # The core calculation: multiply velocities by charges
            modified_velocities = velocities * charges
            
            # It's good practice to work on a copy of the timestep
            new_ts = ts.copy()
            
            # Replace the velocities in the copied frame
            new_ts.velocities = modified_velocities
            
            # Write the entire modified frame to the new file
            writer.write(new_ts)
            
            # Progress indicator
            if (i + 1) % 100 == 0 or (i + 1) == n_frames:
                print(f"  Processed frame {i + 1} of {n_frames}")

    print("\nModification complete. The output file has been saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Modify a GROMACS .trr file by replacing velocities (v) with charge-weighted velocities (q*v)."
    )
    parser.add_argument("--tpr", required=True, help="Input .tpr file (e.g., md.tpr)")
    parser.add_argument("--traj", required=True, help="Input trajectory .trr file (e.g., md.trr)")
    parser.add_argument("--out", required=True, help="Output modified .trr file (e.g., md_qv.trr)")
    
    args = parser.parse_args()
    
    modify_trajectory(args.tpr, args.traj, args.out)