# Subset a fasta file to remove scaffolds < 3Mbp for RepeatOBserver

from Bio import SeqIO 

# create variablies 
input_file = "/g/data/ht96/McLay_UQ/inversion_paper/centromeres/2_D01_30.fasta"
output_file = "/g/data/ht96/McLay_UQ/inversion_paper/centromeres/2_D01_24.fasta"
scaffold_limit = 24

# open the file, and write to the output file only the scaffolds within the limit 
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    records = list(SeqIO.parse(infile, "fasta"))
    selected_records = records[:scaffold_limit]
    SeqIO.write(selected_records, outfile, "fasta")

