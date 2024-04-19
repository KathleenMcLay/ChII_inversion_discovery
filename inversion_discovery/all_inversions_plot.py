import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import numpy as np

# Read data from CSV
data = pd.read_csv("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/data/inversion_list/inversions_all.csv")

# Define unique colors for each population
population_colors = {
    "D00H00": "darkgreen",
    "D01H01": "darkseagreen",
    "D03H02": "forestgreen",
    "D04H05": "seagreen",
    "D05H06": "olive",
    "D35D09": "goldenrod",
    "D32H12": "black",
    "H03": "#55C667FF",
    "H04": "#95D840FF", 
    "H07": "#DCE319FF", 
    "H14": "darkolivegreen", 
    "D12": "burlywood",
    "A371011": "#29AF7FFF",
}

# Define the order of populations
population_order = [
    "D00H00", "D01H01", "D03H02", "D04H05", "D05H06", "D35D09", 
    "D32H12", "H03", "H04", "H07", "H14", "D12", "A371011"
]

# Group data by Scaffold
grouped_data = data.groupby("Scaffold")

# Iterate over each scaffold group
for scaffold, group in grouped_data:
    # Skip scaffold if no data exists for it
    if group.empty:
        continue
    
    # Calculate the height ratio for the current scaffold
    height_ratio = len(group['Population'].unique())
    
    # Create a new figure with a specified width and dynamic height based on the number of ticks
    fig, ax = plt.subplots(figsize=(group['chr_end'].max() / 20000000, height_ratio/2.5))
    
    # Sort data for the current scaffold
    group_sorted = group.sort_values("Start (bp)")
    
    # Filter and reorder data based on population order
    group_filtered = group_sorted[group_sorted["Population"].isin(population_order)]
    group_filtered["Population"] = pd.Categorical(group_filtered["Population"], categories=population_order, ordered=True)
    group_filtered = group_filtered.sort_values("Population")
    
    # Set start and finish based on chr_start and chr_end
    start = group_filtered['chr_start'].min()
    end = group_filtered['chr_end'].max()

    # Plot bars and horizontal grid lines for each data point
    for index, row in group_filtered.iterrows():
        # Get color for the current population
        color = population_colors.get(row["Population"], "gray")
        
        # Plot the bar
        bar = ax.barh(row["Population"], row["End (bp)"] - row["Start (bp)"], left=row["Start (bp)"], height=0.5, color=color)

        # Draw a horizontal grid line at the y-axis tick positions
        ax.axhline(bar[0].get_y(), color=color, linestyle='--', linewidth=0.5)

    # Set x-axis limit and ticks
    ax.set_xlim(start, end)
    ax.set_xticks(range(start, end + 1000000, 50000000))  # Ticks every 50 Mb
    ax.set_xticklabels([x // 1000000 for x in range(start, end + 1000000, 50000000)])  # Labels every 50 Mb, divided by 1,000,000

    # Set y-axis title without underscore and with capitalization
    scaffold_title = scaffold.replace("_", " ").capitalize()
    ax.set_ylabel(scaffold_title)
    #ax.yaxis.set_label_coords(-0.12, 0.5)

    # Set labels and title
    ax.set_xlabel('Position (Mb)')

    # Save each plot separately
    plt.savefig(f'{scaffold}_inversions.tiff', dpi=300, format='tiff')
    # Show the plots
    plt.show()
    plt.close()