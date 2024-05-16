#!/usr/bin/env python 

import pandas as pd
import argparse

if __name__ == '__main__':

    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--pc_file", help="Path to the file containing the first two PCs (REQUIRED)")
    parser.add_argument("--z_score_threshold", help="Z-score threshold (REQUIRED)")
    args = parser.parse_args()

    # Get input directory
    pc_file = args.pc_file
    if pc_file == None:
        raise ValueError('Please provide a valid path to file with the first two PCs')

    # Get cohort name
    z_score_threshold = int(args.z_score_threshold)
    if z_score_threshold == None:
        raise ValueError('Please provide a valid z-score')

    # Load PCs
    df = pd.read_csv('pc1_2.txt', sep='\t')

    # Calculate z-score for first column
    df['z'] = abs((df['Comp1'] - df['Comp1'].mean()) / df['Comp1'].std())

    # Filter out samples with z-score above threshold
    df_no_outliers = df[df['z'] < z_score_threshold]

    # Write filtered df to file
    df_no_outliers.to_csv('pc1_2_no_outliers.txt', sep='\t', index=False)

    # Create list of remaining samples
    remaining_samples = df_no_outliers['Sample'].tolist()
    remaining_samples = [f'{s}\n' for s in remaining_samples]

    # Write list of remaning samples out to file
    with open('passed_samples.txt', 'w') as remaining_samples_file:
        remaining_samples_file.writelines(remaining_samples)

