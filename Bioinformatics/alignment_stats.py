### Consolidate the alignment stat output files from samtools flagstat 

import pandas as pd
import os

# Create an empty DataFrame
df = pd.DataFrame()

# Path to the directory containing *stats.txt files
directory_path = '/g/data/ht96/McLay_UQ/inversion_paper/2_align/'

# Iterate through files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith("stats.txt"):
        file_path = os.path.join(directory_path, filename)

        # Extract header from the file name
        header = filename.split('_')[0]

        # Read contents of the file and extract numbers before the first space
        with open(file_path, 'r') as file:
            numbers = [line.split(' ')[0] for line in file.readlines()]

        # Add a new column to the DataFrame with the extracted numbers
        df[header] = numbers

# Display the resulting DataFrame
df.to_csv('/home/564/km6006/Scripts/inversion_paper/alignment_stats.csv')