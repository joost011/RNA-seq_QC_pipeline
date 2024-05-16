#!/usr/bin/env python 

import sys
from scipy.special import ndtri
import numpy as np
import pandas as pd


def force_normalise(x, axis):
    '''Function that forces a normal distribution on a matrix'''
    return ndtri((pd.DataFrame(x).rank(axis=axis, ascending=True) - 0.5) / x.shape[axis])


if __name__ == '__main__':

    # Check if the expected number of input parameters is defined
    if len(sys.argv) < 3:
        print("Usage: infile outfile")
        sys.exit()

    # Assign input parameters to variables
    infile=sys.argv[1]
    outfile=sys.argv[2]

    # If parameter does not contain gzip extension, add it
    if not outfile.endswith(".gz"):
        outfile = outfile + ".gz"

    # Read in dataframe
    matrix = pd.read_csv(infile, sep='\t')

    # Set the index to the first '-' column
    matrix.set_index('-', inplace=True)

    # Normalize the data
    normalized = force_normalise(matrix.values, 1)

    # Convert matrix back to dataframe and add columns and an index
    final_matrix = pd.DataFrame(normalized)
    final_matrix.columns = matrix.columns
    final_matrix['gene'] = matrix.index
    final_matrix.set_index('gene', inplace=True)

    # Write dataframe out to file
    final_matrix.to_csv(outfile, compression='gzip', sep='\t')