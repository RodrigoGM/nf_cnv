#!/usr/bin/env python3
"""
Aggregate seg.txt files into matrices for CNV analysis
Fast matrix generation using pandas for large datasets
"""

import pandas as pd
import numpy as np
import argparse
import glob
import os
import sys
from pathlib import Path

def extract_cell_id(filename):
    """Extract cell ID from filename"""
    base = os.path.basename(filename)
    # Remove _seg.txt and extract cell_id part
    cell_id = base.replace('_seg.txt', '')
    # Remove resolution and genome parts if present
    parts = cell_id.split('.')
    if len(parts) >= 4:  # e.g., LS0000_A_1201_b1.hsa38.PE_FW.5k.bwa
        return parts[0]  # Return just the cell ID part
    return cell_id

def main():
    parser = argparse.ArgumentParser(description='Aggregate seg.txt files into matrices')
    parser.add_argument('--input-dir', required=True, help='Directory with seg files')
    parser.add_argument('--resolution', required=True, help='Resolution (5, 20, 50)')
    parser.add_argument('--output-prefix', default='.', help='Output prefix/directory')
    
    args = parser.parse_args()
    
    # Find all seg files for this resolution
    pattern = f"*{args.resolution}k.bwa_seg.txt"
    seg_files = glob.glob(os.path.join(args.input_dir, pattern))
    
    if not seg_files:
        print(f"No seg files found with pattern: {pattern}")
        sys.exit(1)
    
    print(f"Processing {len(seg_files)} seg files for {args.resolution}k resolution")
    
    # Read first file to get chromosome info and establish n_bins
    first_file = seg_files[0]
    try:
        df_first = pd.read_csv(first_file, sep='\t')
        n_bins = len(df_first)
        print(f"Expected number of bins: {n_bins}")
    except Exception as e:
        print(f"Error reading first file {first_file}: {e}")
        sys.exit(1)
    
    # Create chromosome info file
    chrom_info = df_first[['chr', 'start', 'end', 'gc', 'chr.arm']].copy()
    chrom_info_file = f"{args.output_prefix}/chrom_info_{args.resolution}k.txt"
    chrom_info.to_csv(chrom_info_file, sep='\t', index=False)
    print(f"Wrote chromosome info: {chrom_info_file}")
    
    # Initialize data collectors for each column
    data_collectors = {
        'count': [],
        'ratio': [],
        'lowess.ratio': [],
        'seg.mean': [],
        'lowess.ratio.quantal': [],
        'seg.mean.quantal': []
    }
    cell_ids = []

    # Process each file and collect data
    for seg_file in seg_files:
        cell_id = extract_cell_id(seg_file)

        try:
            df = pd.read_csv(seg_file, sep='\t')

            # Verify same number of bins
            if len(df) != n_bins:
                print(f"Warning: {seg_file} has {len(df)} bins, expected {n_bins}")
                continue

            # Collect data for each metric
            cell_ids.append(cell_id)
            data_collectors['count'].append(df['count'].values)
            data_collectors['ratio'].append(df['ratio'].values)
            data_collectors['lowess.ratio'].append(df['lowess.ratio'].values)
            data_collectors['seg.mean'].append(df['seg.mean'].values)
            data_collectors['lowess.ratio.quantal'].append(df['lowess.ratio.quantal'].values)
            data_collectors['seg.mean.quantal'].append(df['seg.mean.quantal'].values)

        except Exception as e:
            print(f"Error processing {seg_file}: {e}")
            continue

    # Check if we have any valid data
    if not cell_ids:
        print("No valid seg files processed")
        sys.exit(1)

    print(f"Successfully processed {len(cell_ids)} cells")

    # Write matrices to files (memory-efficient approach)
    matrix_names = {
        'count': 'matrix_counts',
        'ratio': 'matrix_ratio', 
        'lowess.ratio': 'matrix_lowess_ratio',
        'seg.mean': 'matrix_seg_mean',
        'lowess.ratio.quantal': 'matrix_lowess_ratio_quantal',
        'seg.mean.quantal': 'matrix_seg_mean_quantal'
    }
    
    # Process each matrix type separately to reduce memory usage
    for key, data_list in data_collectors.items():
        if data_list:
            print(f"Processing {key} matrix...")
            matrix_data = np.column_stack(data_list)
            matrix = pd.DataFrame(matrix_data, columns=cell_ids)
            
            filename = f"{args.output_prefix}/{matrix_names[key]}_{args.resolution}k.txt"
            matrix.to_csv(filename, sep='\t', index=False)
            print(f"Wrote {filename}: {matrix.shape[0]} bins x {matrix.shape[1]} cells")
            
            # Free memory immediately
            del matrix_data, matrix
            import gc; gc.collect()
        else:
            print(f"No data collected for {key}")

    print(f"Matrix aggregation completed for {args.resolution}k resolution")

if __name__ == '__main__':
    main()

# __EOF__

##% ommited -- high memory usage
##%     # Create matrices using efficient array stacking
##%     matrices = {}
##%     matrix_names = {
##%         'count': 'matrix_counts',
##%         'ratio': 'matrix_ratio',
##%         'lowess.ratio': 'matrix_lowess_ratio',
##%         'seg.mean': 'matrix_seg_mean',
##%         'lowess.ratio.quantal': 'matrix_lowess_ratio_quantal',
##%         'seg.mean.quantal': 'matrix_seg_mean_quantal'
##%     }
##% 
##%     for key, data_list in data_collectors.items():
##%         if data_list:
##%             try:
##%                 # Stack arrays and create DataFrame in one step
##%                 matrix_data = np.column_stack(data_list)
##%                 matrices[key] = pd.DataFrame(matrix_data, columns=cell_ids)
##%                 
##%                 # Write matrix to file
##%                 filename = f"{args.output_prefix}/{matrix_names[key]}_{args.resolution}k.txt"
##%                 matrices[key].to_csv(filename, sep='\t', index=False)
##%                 print(f"Wrote {filename}: {matrices[key].shape[0]} bins x {matrices[key].shape[1]} cells")
##%                 
##%             except Exception as e:
##%                 print(f"Error creating matrix for {key}: {e}")
##%                 continue
