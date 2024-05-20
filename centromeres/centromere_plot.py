import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Data
data = {
    "scaffold": ["1", "2", "3", "4", "5",
                 "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15"],
    "centromere": [8.23, 182.36, 56.26, 74.52, 170.72, 36.54, 1.59, 124.89, 78.57, 72.89,
                   36.54, 55.79, 52.42, 17.39, 14.82],
    "chr_start": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "chr_end": [271.4, 256.6, 224.5, 224.5, 216.8, 178.1, 176.4, 172.3, 139, 132.8,
                90.8, 64.3, 53.4, 49, 33.3]
}

# Plotting
plt.figure(figsize=(14, 10))

# Plotting bars
for i in range(len(data["scaffold"])):
    #plt.plot([i, i], [data["chr_start"][i], data["chr_end"][i]], color='#C1C4AD', linewidth=10)  # scaffold bars
    plt.bar(i, data["chr_end"][i] - data["chr_start"][i], bottom=data["chr_start"][i],
            edgecolor='black', linewidth=0.5, color='none', width=0.5) 
    plt.plot(i, data["centromere"][i], marker='o', color='#596D22', markersize=10, label='Centromere')  # centromere point

# Formatting
plt.ylabel('Position (Mb)', fontsize=18, fontname='Helvetica', labelpad=20)
plt.xlabel('Scaffold', fontsize=18, fontname='Helvetica', labelpad=20)

plt.xticks(range(len(data["scaffold"])), data["scaffold"])
plt.xlim(-0.75, len(data["scaffold"]) - 0.5)
plt.ylim(min(data["chr_start"]) - 15, max(data["chr_end"]) + 10)  # Adjusted y-axis limit

plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
#plt.gca().spines['bottom'].set_visible(False)

plt.tick_params(axis='x', which='both', bottom=False)

legend_elements = [Line2D([0], [0], marker='o', color='#596D22', markersize=10, linestyle='', label='Predicted centromere')]

# Show legend
plt.legend(handles=legend_elements, frameon=False, prop={'family': 'Helvetica', 'size': 16})

plt.savefig('scaffold_visualization.tiff', format='tiff', bbox_inches='tight', dpi=300)

# Displaying the plot
plt.show()