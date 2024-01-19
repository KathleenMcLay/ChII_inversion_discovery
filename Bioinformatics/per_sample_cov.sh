cd /g/data/ht96/McLay_UQ/inversion_paper/3_clean
# # Loop through the directory and process text files with "coverage" in the name
for file in *coverage*.txt; do
    if [[ -f $file ]]; then
        # Extract the first part of the file name before the first _
        file_name=$(echo "$file" | cut -d'_' -f1)
        # Calculate the average of the values in the meandepth column (exclude Chr 13)
        average_mean_depth=$(awk '{ if(NR > 1 && $1 != "" && $1 != "SLv141Ch13" && $7 != "" && $7 != 0) { sum+=$7; count++ } } END { if(count > 0) print sum/count; else print 0 }' "$file")

        # Add data to the output file
        echo "$file_name $average_mean_depth" >> "$output_file"
    fi
done

### Consolidate per chromosome coverage 

chrom_file="/home/564/km6006/Scripts/inversion_paper/per_chrom_cov.txt"

# Create a temporary file to store intermediate results
temp_file=$(mktemp)

# Loop through coverage files
for file in *coverage*.txt; do
    if [[ -f $file ]]; then
        # Extract the file name
        file_name=$(echo "$file" | cut -d'_' -f1)
        # Extract the whole 7th column without the header and store it in the temporary file
        awk '{print $7}' "$file" > "${temp_file}" 
        # Replace the header with "file_name"
        sed -i "1s/.*/$file_name/" "${temp_file}"

        # Use paste to merge the temporary file with the output file
        paste -d "\t" "${chrom_file}" "${temp_file}" > "${chrom_file}.new"

        # Move the new file to the original output file
        mv "${chrom_file}.new" "${chrom_file}"
    fi
done

# Add a final column with the average of each row
awk '{ sum = 0; for(i=2; i<=NF; i++) sum += $i; avg = sum / (NF-1); print $0 "\t" avg }' "${chrom_file}" > "${chrom_file}.new"

# Move the new file with the final column to the original output file
mv "${chrom_file}.new" "${chrom_file}"

# Remove the temporary file
rm "${temp_file}"

