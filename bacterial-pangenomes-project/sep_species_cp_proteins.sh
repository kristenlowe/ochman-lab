# Define the source and base destination directories
src_dir="/stor/scratch/Ochman/kristen/bacterial_species_proteins"
base_dest_dir="/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species"

# Loop through each protein fasta file in the source directory
for file in ${src_dir}/*_protein.faa; do
    # Extract the bacterial species name from the file name
    species=$(basename "$file" | cut -d'_' -f1-2)
    
    # Define the destination directory for this species
    species_dest_dir="${base_dest_dir}/${species}/${species}_proteins"
    
    # Create the species-specific directory if it doesn't already exist
    mkdir -p "$species_dest_dir"
    
    # Copy the file to the species-specific directory
    cp "$file" "$species_dest_dir"
done

