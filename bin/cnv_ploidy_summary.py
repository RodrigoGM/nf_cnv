#!/usr/bin/env python3
"""
CNV Ploidy Statistics Calculator - Fixed Version
Handles irregular tab spacing in TSV files
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
import re

def calculate_stats(data, column):
    """Calculate statistics for a given column"""
    # Convert to string first, then numeric, handling various formats
    if column not in data.columns:
        print(f"Warning: Column '{column}' not found in data", file=sys.stderr)
        return {
            'count': 0,
            'mean': 'NA',
            'median': 'NA',
            'min': 'NA',
            'max': 'NA',
            'std': 'NA'
        }
    
    # Clean the data - remove extra spaces and convert to numeric
    clean_series = data[column].astype(str).str.strip()
    clean_series = clean_series.replace('', np.nan)  # Replace empty strings with NaN
    clean_data = pd.to_numeric(clean_series, errors='coerce').dropna()
    
    if len(clean_data) == 0:
        return {
            'count': 0,
            'mean': 'NA',
            'median': 'NA',
            'min': 'NA',
            'max': 'NA',
            'std': 'NA'
        }
    
    return {
        'count': len(clean_data),
        'mean': round(clean_data.mean(), 4),
        'median': round(clean_data.median(), 4),
        'min': round(clean_data.min(), 4),
        'max': round(clean_data.max(), 4),
        'std': round(clean_data.std(), 4)
    }

def clean_dataframe(df):
    """Clean dataframe by removing extra whitespace and empty columns"""
    # Strip whitespace from all string columns
    for col in df.columns:
        if df[col].dtype == 'object':
            df[col] = df[col].astype(str).str.strip()
    
    # Remove completely empty columns
    df = df.dropna(axis=1, how='all')
    
    # Remove rows where all values are empty/NaN
    df = df.dropna(axis=0, how='all')
    
    return df

def print_summary_table(stats_dict, title, output_file=sys.stdout):
    """Print formatted summary statistics"""
    print(f"\n{title}", file=output_file)
    print("=" * len(title), file=output_file)
    print(f"{'Metric':<10} {'Count':<8} {'Mean':<10} {'Median':<10} {'Min':<10} {'Max':<10} {'Std':<10}", file=output_file)
    print("-" * 68, file=output_file)
    
    for metric, stats in stats_dict.items():
        print(f"{metric:<10} {stats['count']:<8} {stats['mean']:<10} {stats['median']:<10} "
              f"{stats['min']:<10} {stats['max']:<10} {stats['std']:<10}", file=output_file)

def main():
    parser = argparse.ArgumentParser(
        description='Calculate statistics for CNV ploidy results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cnv_stats.py ploidy_results_combined.txt
  python cnv_stats.py -i data.txt -o summary.txt --group-by resolution
  python cnv_stats.py -i data.txt --group-by resolution,genome
        """
    )
    
    parser.add_argument('input', nargs='?', 
                       help='Input file (TSV format, default: stdin)')
    parser.add_argument('-i', '--input-file', 
                       help='Input file path (alternative to positional argument)')
    parser.add_argument('-o', '--output', 
                       help='Output file (default: stdout)')
    parser.add_argument('--group-by', 
                       help='Group by columns (comma-separated): resolution,genome,read')
    parser.add_argument('--format', choices=['table', 'csv', 'json'], 
                       default='table', help='Output format')
    parser.add_argument('--debug', action='store_true',
                       help='Print debug information')
    
    args = parser.parse_args()
    
    # Determine input source
    if args.input_file:
        input_file = args.input_file
    elif args.input:
        input_file = args.input
    else:
        input_file = sys.stdin
    
    # Read data with better parsing for irregular spacing
    try:
        if input_file == sys.stdin:
            df = pd.read_csv(sys.stdin, sep='\t', skipinitialspace=True)
        else:
            # Try multiple parsing strategies
            try:
                # First try: normal tab-separated
                df = pd.read_csv(input_file, sep='\t', skipinitialspace=True)
            except:
                # Second try: handle multiple tabs/spaces
                with open(input_file, 'r') as f:
                    lines = f.readlines()
                
                # Clean up the lines - replace multiple tabs/spaces with single tab
                cleaned_lines = []
                for line in lines:
                    # Replace multiple whitespace with single tab
                    cleaned_line = re.sub(r'\s+', '\t', line.strip())
                    cleaned_lines.append(cleaned_line)
                
                # Create dataframe from cleaned lines
                from io import StringIO
                df = pd.read_csv(StringIO('\n'.join(cleaned_lines)), sep='\t')
                
    except Exception as e:
        print(f"Error reading input: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Clean the dataframe
    df = clean_dataframe(df)
    
    if args.debug:
        print(f"DataFrame shape: {df.shape}", file=sys.stderr)
        print(f"Columns: {df.columns.tolist()}", file=sys.stderr)
        print(f"First few rows:\n{df.head()}", file=sys.stderr)
        print(f"Data types:\n{df.dtypes}", file=sys.stderr)
        print(f"Sample ploidy values: {df['ploidy'].head().tolist() if 'ploidy' in df.columns else 'Column not found'}", file=sys.stderr)
        print(f"Sample error values: {df['error'].head().tolist() if 'error' in df.columns else 'Column not found'}", file=sys.stderr)
    
    # Validate required columns
    required_cols = ['cellID', 'ploidy', 'error']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}", file=sys.stderr)
        print(f"Available columns: {df.columns.tolist()}", file=sys.stderr)
        sys.exit(1)
    
    # Setup output
    if args.output:
        output_file = open(args.output, 'w')
    else:
        output_file = sys.stdout
    
    try:
        # Overall statistics
        overall_stats = {
            'Ploidy': calculate_stats(df, 'ploidy'),
            'Error': calculate_stats(df, 'error')
        }
        
        if args.format == 'table':
            print("CNV Ploidy Analysis Results", file=output_file)
            print("=" * 40, file=output_file)
            print(f"Total cells analyzed: {len(df)}", file=output_file)
            
            print_summary_table(overall_stats, "Overall Statistics", output_file)
            
            # Group by analysis
            if args.group_by:
                group_cols = [col.strip() for col in args.group_by.split(',')]
                valid_group_cols = [col for col in group_cols if col in df.columns]
                
                if valid_group_cols:
                    print(f"\nGrouped Analysis by: {', '.join(valid_group_cols)}", file=output_file)
                    print("=" * 50, file=output_file)
                    
                    # Clean group columns
                    for col in valid_group_cols:
                        df[col] = df[col].astype(str).str.strip()
                    
                    for group_name, group_data in df.groupby(valid_group_cols):
                        if len(valid_group_cols) == 1:
                            group_title = f"{valid_group_cols[0]}: {group_name}"
                        else:
                            group_title = f"{dict(zip(valid_group_cols, group_name))}"
                        
                        group_stats = {
                            'Ploidy': calculate_stats(group_data, 'ploidy'),
                            'Error': calculate_stats(group_data, 'error')
                        }
                        
                        print(f"\n{group_title} (n={len(group_data)})", file=output_file)
                        print("-" * len(group_title), file=output_file)
                        print(f"{'Metric':<8} {'Count':<6} {'Mean':<8} {'Median':<8} {'Min':<8} {'Max':<8} {'Std':<8}", file=output_file)
                        print("-" * 54, file=output_file)
                        
                        for metric, stats in group_stats.items():
                            print(f"{metric:<8} {stats['count']:<6} {stats['mean']:<8} {stats['median']:<8} "
                                  f"{stats['min']:<8} {stats['max']:<8} {stats['std']:<8}", file=output_file)
        
        elif args.format == 'csv':
            # CSV output implementation
            results = []
            for metric, stats in overall_stats.items():
                results.append({
                    'group': 'Overall',
                    'metric': metric,
                    **stats
                })
            
            if args.group_by:
                group_cols = [col.strip() for col in args.group_by.split(',')]
                valid_group_cols = [col for col in group_cols if col in df.columns]
                
                if valid_group_cols:
                    for group_name, group_data in df.groupby(valid_group_cols):
                        group_label = f"{dict(zip(valid_group_cols, group_name))}" if len(valid_group_cols) > 1 else str(group_name)
                        
                        for metric in ['ploidy', 'error']:
                            stats = calculate_stats(group_data, metric)
                            results.append({
                                'group': group_label,
                                'metric': metric.title(),
                                **stats
                            })
            
            pd.DataFrame(results).to_csv(output_file, index=False)
        
    finally:
        if args.output:
            output_file.close()

if __name__ == '__main__':
    main()
