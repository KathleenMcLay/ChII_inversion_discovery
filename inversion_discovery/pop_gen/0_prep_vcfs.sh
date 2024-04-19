module load bcftools

# convert bcf to vcf and index 
directory="/QRISdata/Q6656/new_lpca/"


for file in $directory/*.bcf; do 
    filename=$(basename -- "$file")
    filename="${filename%.*}"
    bcftools view -O z -o "$directory/$filename.vcf.gz" "$file"
    bcftools index -t "$directory/$filename.vcf.gz"
done


# create a vcf for each populations inversions 
pops="/home/uqkmcla4/scripts/current/pops.txt"
region_dir="/QRISdata/Q6656/new_lpca/regions"

while IFS= read -r pop; do
    vcf="/QRISdata/Q6656/new_lpca/population_vcfs/$pop.vcf.gz" 
    for file in "$region_dir"/"${pop}"*_regions.txt; do 
        if [ -f "$file" ]; then
            while IFS= read -r row || [ -n "$row" ]; do
                inversion=$(echo "$row" | cut -d ' ' -f1) 
                fi
            done < "$file"
        else 
            echo "Error: $file does not exist"
        fi
    done
done < "$pops"