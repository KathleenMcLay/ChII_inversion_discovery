### Merge the samtools per sample coverage files into a single file

output_file="/g/data/ht96/McLay_UQ/inversion_paper/per_sample_coverage.txt"

cd /g/data/ht96/McLay_UQ/inversion_paper/3_clean
# Loop through the directory and process text files with "table" in the name
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
