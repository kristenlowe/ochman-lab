#!/bin/bash

# Usage: ./seqprop_pipeline_each_species.sh

BASE_DIR="/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species/"

# Loop through each species directory
for species_dir in "$BASE_DIR"/*; do
    species=$(basename "$species_dir")
    cds_fasta="$species_dir/rep_${species}_cds.fna"
    protein_fasta="$species_dir/rep_${species}_proteins.faa"
    alias="$species"

    # Check if both the CDS and protein fasta files exist
    if [[ -f "$cds_fasta" && -f "$protein_fasta" ]]; then
        # Linearize the CDS fasta
        linear_cds_fasta="$species_dir/linear_rep_${species}_cds.fna"
        awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq$0} END {if (seq) print seq}' "$cds_fasta" > "$linear_cds_fasta"

        # Linearize the protein fasta
        linear_protein_fasta="$species_dir/linear_rep_${species}_proteins.faa"
        awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq$0} END {if (seq) print seq}' "$protein_fasta" > "$linear_protein_fasta"

        # Define the output directory for sequence properties
        output_dir="$species_dir/rep_seq_properties"
        mkdir -p "$output_dir"

        # Calculate sequence properties
        # Length calculation
        awk -F '\n' '{OFS=","; getline seq; match($1, /^>(.*)/, arr); print arr[1], "length", length(seq)}' "$linear_cds_fasta" > "${output_dir}/${alias}_length.csv"

        # GC content
        awk 'BEGIN {OFS=","} /^>/ {if (seq) {gc_content = (gc_count / total_length) * 100; if (header_id != "") print header_id, "GC", gc_content}; seq = ""; gc_count = 0; total_length = 0; match($0, /^>(.*)/, arr); header_id = arr[1]; next} {seq = seq $0; total_length += length($0); gc_count += gsub(/[GCgc]/, "")} END {if (seq && header_id != "") {gc_content = (gc_count / total_length) * 100; print header_id, "GC", gc_content}}' "$linear_cds_fasta" > "${output_dir}/${alias}_gc.csv"

        # GC content at 3rd codon position
        awk 'BEGIN {OFS=","} /^>/ {if (NR > 1 && seq) {gc_3rd_content = (third_gc_count / (total_length / 3)) * 100; if (header_id != "") print header_id, "GC_3rd", gc_3rd_content}; seq = ""; third_gc_count = 0; total_length = 0; match($0, /^>(.*)/, arr); header_id = arr[1]; next} {seq = seq $0; total_length += length($0); for (i = 3; i <= length($0); i += 3) {if (substr($0, i, 1) ~ /[GCgc]/) third_gc_count++}} END {if (seq && header_id != "") {gc_3rd_content = (third_gc_count / (total_length / 3)) * 100; print header_id, "GC_3rd", gc_3rd_content}}' "$linear_cds_fasta" > "${output_dir}/${alias}_gc_3rd.csv"

        # Polar amino acid content
        awk 'BEGIN {OFS=","} /^>/ {if (NR > 1) {polar_content = (polar_count / total_length) * 100; print identifier, "polarAA", polar_content} identifier = substr($1, 2); total_length = 0; polar_count = 0; next} {total_length += length($0); polar_count += gsub(/[STYECQNHKR]/, "", $0)} END {if (total_length > 0) {polar_content = (polar_count / total_length) * 100; print identifier, "polarAA", polar_content}}' "$linear_protein_fasta" > "${output_dir}/${alias}_polar_aa.csv"

        # Hydrophobic amino acid content
        awk 'BEGIN {OFS=","} /^>/ {if (NR > 1) {hydrophobic_content = (hydrophobic_count / total_length) * 100; print identifier, "hydrophobicAA", hydrophobic_content} identifier = substr($1, 2); total_length = 0; hydrophobic_count = 0; next} {total_length += length($0); hydrophobic_count += gsub(/[ACFGILMPVWY]/, "", $0)} END {if (total_length > 0) {hydrophobic_content = (hydrophobic_count / total_length) * 100; print identifier, "hydrophobicAA", hydrophobic_content}}' "$linear_protein_fasta" > "${output_dir}/${alias}_hydrophobic_aa.csv"

        # Metabolic costs
        awk 'BEGIN {OFS=","} NR==FNR {costs[$1]=$2; next} /^>/ {if (NR > FNR + 1 && total_length > 0) {print identifier, "metabol", total_cost / total_length} identifier = substr($1, 2); total_length = 0; total_cost = 0; next} {for (i = 1; i <= length($0); i++) {aa = substr($0, i, 1); if (aa in costs) {total_cost += costs[aa]}} total_length += length($0)} END {if (total_length > 0) {print identifier, "metabol", total_cost / total_length}}' /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/AA_metabolic_costs.txt "$linear_protein_fasta" > "${output_dir}/${alias}_metabol.csv"

        # CAI calculation
        /stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/cai -seqall "$linear_cds_fasta" -cfile /stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/data/CODONS/Eecoli.cut -outfile "${output_dir}/${alias}.cai"
        awk '{print $2, "cai", $4}' OFS="," "${output_dir}/${alias}.cai" > "${output_dir}/${alias}_cai.csv"

        # Instability calculation
        python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/calculate_instability.py "$linear_protein_fasta" "${output_dir}/${alias}_instability.csv"

        # Intrinsic structural disorder calculation
        python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/iupred3/iupred3.py "$linear_protein_fasta" long > "${output_dir}/${alias}_iupred3_output.txt"
        python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/calculate_disorder.py "$linear_protein_fasta" "${output_dir}/${alias}_iupred3_output.txt" "${output_dir}/${alias}_disorder.csv"

        # Transmembrane domain calculation
        /stor/work/Ochman/hassan/tools/phobius/phobius.pl -short "$linear_protein_fasta" > "${output_dir}/${alias}_tm.phobius"
        python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/tm_to_csv.py "${output_dir}/${alias}_tm.phobius" "${output_dir}/${alias}_tm.csv"

        # Transmembrane coverage calculation
        python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/calculate_tm_coverage.py "${output_dir}/${alias}_tm.phobius" "${output_dir}/${alias}_length.csv" "${output_dir}/${alias}_tm_coverage.csv"

        # Signal peptide presence/absence
        python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/sp_to_csv.py "${output_dir}/${alias}_tm.phobius" "${output_dir}/${alias}_sp.csv"

        # Amino acid composition bias calculation
        Rscript /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/calc_aa_comp_bias.R "$linear_protein_fasta" "${output_dir}/${alias}_aa_comp_bias.csv"

        # Concatenate all CSV files into one
        cat "${output_dir}/${alias}"_*.csv > "${output_dir}/${alias}_seq_properties.csv"

    fi
done

echo "Sequence properties calculation completed for all species."
