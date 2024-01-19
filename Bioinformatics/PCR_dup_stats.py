### Extract Picard MarkDuplicates stats and save as table

import pandas as pd
import os
import linecache

# Create an empty DataFrame
df = pd.DataFrame()

# Path to the directory containing *stats.txt files
directory_path = '/g/data/ht96/McLay_UQ/inversion_paper/3_clean/'

# Iterate through files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith("metrics"):
        file_path = os.path.join(directory_path, filename)
        print(file_path)
        
        # Read the 8th line from each file
        eighth_line = linecache.getline(file_path, 8).strip()

        # Split the 8th line by tabs to get each column
        columns_data = eighth_line.split('\t')

        # Create a DataFrame with the 8th line data
        df = pd.concat([df, pd.DataFrame([columns_data], columns=range(len(columns_data)))], ignore_index=True)

# Save the DataFrame to a CSV file
df.to_csv('/home/564/km6006/Scripts/inversion_paper/PCRdup_stats.tsv', sep='\t', index=False, header=None)
