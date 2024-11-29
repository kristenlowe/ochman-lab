# Define the base directory containing all bacterial species folders
base_dir="/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species"

# Loop through each bacterial species directory
for species_dir in ${base_dir}/*; do
    # Extract the species name from the directory path
    species=$(basename "$species_dir")
    
    # Define the path to the proteins directory
    proteins_dir="${species_dir}/${species}_proteins"
    concatenated_protein_file="${proteins_dir}/all_${species}_proteins.faa"
    
    # Check if the concatenated protein file exists
    if [[ -f "$concatenated_protein_file" ]]; then
        # Create the clustering directory within the proteins directory
        clustering_dir="${proteins_dir}/clustering"
        mkdir -p "$clustering_dir"
        
        # Define the database and output files within the clustering directory
        db="${clustering_dir}/${species}_protein_db"
        result_db="${clustering_dir}/resultDB"
        cluster_db="${clustering_dir}/clusterDB"
        tmp="${clustering_dir}/tmp"
        output_tsv="${clustering_dir}/clusters.tsv"
        
        # Run MMSeqs2 commands
        # Step 1: Create the database
        mmseqs createdb "$concatenated_protein_file" "$db"
        
        # Step 2: Search for homologs
        mmseqs search "$db" "$db" "$result_db" "$tmp" --min-seq-id 0.8 -c 0.8 --cov-mode 1
        
        # Step 3: Convert results to human-readable format
        mmseqs convertalis "$db" "$db" "$result_db" "${clustering_dir}/resultDB.m8"
        
        # Step 4: Cluster homologous proteins using Linclust
        mmseqs linclust "$db" "$cluster_db" "$tmp" --min-seq-id 0.8 -c 0.8 --cov-mode 1
        
        # Step 5: Create the clusters.tsv file
        mmseqs createtsv "$db" "$db" "$cluster_db" "$output_tsv"
        
        # Clean up temporary files
        rm -rf "$tmp"
    fi
done

