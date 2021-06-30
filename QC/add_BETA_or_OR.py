import os, sys
import numpy as np
import pandas as pd

filename = sys.argv[1]

df = pd.read_csv(filename, sep='\s+', header=0)

if 'OR' in df.columns:
    df['BETA'] = np.log(df['OR'])
elif 'BETA' in df.columns:
    df['OR'] = np.exp(df['BETA'])
    df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'TEST', 'OBS_CT', 'OR', 'SE', 'T_STAT', 'P', 'BETA']]
else:
    print("No BETA or OR in the summary statistics")

df.to_csv(filename, sep='\t', index=False)