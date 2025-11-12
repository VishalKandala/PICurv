#!/usr/bin/env python3
"""
Convert LES-Flat grid format to PICGRID format
Input: Coordinates stored separately by dimension
Output: PICGRID format with all coordinates per point
"""

import numpy as np

def read_les_flat_grid(filename):
    """Read grid file in LES-Flat format"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # First line: number of blocks
    nblocks = int(lines[0].strip())
    
    # Second line: grid dimensions
    ni, nj, nk = map(int, lines[1].strip().split())
    
    print(f"Grid dimensions: {ni} x {nj} x {nk}")
    
    # Read X coordinates (lines 2 to 2+ni)
    x_coords = np.zeros(ni)
    for i in range(ni):
        vals = list(map(float, lines[2 + i].strip().split()))
        x_coords[i] = vals[0]  # X is in first column
    
    # Read Y coordinates (lines 2+ni to 2+ni+nj)
    y_coords = np.zeros(nj)
    for j in range(nj):
        vals = list(map(float, lines[2 + ni + j].strip().split()))
        y_coords[j] = vals[1]  # Y is in second column
    
    # Read Z coordinates (lines 2+ni+nj to 2+ni+nj+nk)
    z_coords = np.zeros(nk)
    for k in range(nk):
        vals = list(map(float, lines[2 + ni + nj + k].strip().split()))
        z_coords[k] = vals[2]  # Z is in third column
    
    return ni, nj, nk, x_coords, y_coords, z_coords

def write_picgrid_format(filename, ni, nj, nk, x_coords, y_coords, z_coords):
    """Write grid file in PICGRID format"""
    
    total_points = ni * nj * nk
    print(f"Writing {total_points} grid points to {filename}")
    
    with open(filename, 'w') as f:
        # Write header
        f.write("PICGRID\n")
        f.write("1\n")
        f.write(f"{ni} {nj} {nk}\n")
        
        # Write coordinates
        # The ordering appears to be: for i, for k, for j (based on cpipe_coarse.grid)
        # This means i varies fastest, then j, then k
        for k in range(ni):
            for j in range(nk):
                for i in range(nj):
                    x = x_coords[i] + min(x_coords)
                    y = y_coords[j] + min(y_coords)
                    z = z_coords[k] + min(z_coords)
                    f.write(f"{x:.8e} {y:.8e} {z:.8e}\n")
    
    print(f"Successfully wrote {filename}")

def main():
    input_file = "/mnt/user-data/uploads/LES-Flat_grid.dat"
    output_file = "/mnt/user-data/outputs/LES-Flat_PICGRID.grid"
    
    print(f"Reading {input_file}...")
    ni, nj, nk, x_coords, y_coords, z_coords = read_les_flat_grid(input_file)
    
    print(f"\nCoordinate ranges:")
    print(f"  X: [{x_coords.min():.6f}, {x_coords.max():.6f}]")
    print(f"  Y: [{y_coords.min():.6f}, {y_coords.max():.6f}]")
    print(f"  Z: [{z_coords.min():.6f}, {z_coords.max():.6f}]")
    
    print(f"\nWriting to {output_file}...")
    write_picgrid_format(output_file, ni, nj, nk, x_coords, y_coords, z_coords)
    
    print("\nConversion complete!")

if __name__ == "__main__":
    main()
