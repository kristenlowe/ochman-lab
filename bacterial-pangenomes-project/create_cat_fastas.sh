# Define the base directory containing all bacterial species folders
base_dir="/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species"

# Loop through each bacterial species directory
for species_dir in ${base_dir}/*; do
    # Extract the species name from the directory path
    species=$(basename "$species_dir")
    
    # Define the three subdirectories: CDSs, genomes, and proteins
    for subdir in "${species_dir}/${species}_CDSs" "${species_dir}/${species}_genomes" "${species_dir}/${species}_proteins"; do
        # Check if the subdirectory exists to avoid errors
        if [[ -d "$subdir" ]]; then
            # Determine the suffix based on the subdirectory name
            suffix=""
            case "$subdir" in
                *CDSs) suffix="cds.fna" ;;
                *genomes) suffix="genomes.fna" ;;
                *proteins) suffix="proteins.faa" ;;
            esac
            
            # Define the output concatenated file name
            output_file="${subdir}/all_${species}_${suffix}"
            
            # Concatenate all files in the subdirectory into the output file
            cat ${subdir}/*.f* > "$output_file"
        fi
    done
done

