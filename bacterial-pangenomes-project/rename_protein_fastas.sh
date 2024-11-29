# Define the source and destination directories for protein files
src_dir="/stor/scratch/Ochman/kristen/pangenome/pangenome_project_proteins"
dest_dir="/stor/scratch/Ochman/kristen/bacterial_species_proteins"

# Create the destination directory if it doesn't already exist
mkdir -p "$dest_dir"

# Initialize an associative array to track file counts for each species
declare -A species_count

# Loop through each line in the accession-species mapping file
while IFS=$'\t' read -r accession species; do
    # Find the corresponding protein file using the accession number
    file=$(ls ${src_dir}/${accession}_*.faa 2>/dev/null)

    # If a matching file is found, proceed with copying and renaming
    if [[ -n "$file" ]]; then
        # Increment the species count to number the files
        ((species_count["$species"]++))
        
        # Format the new file name with species name, count, and suffix
        new_file="${species}_$(printf "%02d" ${species_count["$species"]})_protein.faa"
        
        # Copy the file to the new directory with the new name
        cp "$file" "${dest_dir}/${new_file}"
    fi
done < /stor/scratch/Ochman/hassan/083024_RefSeq_Genbank_genomes/accession_speciesname.tsv

