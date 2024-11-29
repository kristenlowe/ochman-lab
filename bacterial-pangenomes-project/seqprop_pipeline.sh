#!/bin/bash

# Usage: ./seqprop_pipeline.sh <DNA_FASTA> <PROTEIN_FASTA> <ALIAS>

# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <DNA_FASTA> <PROTEIN_FASTA> <ALIAS>"
    exit 1
fi

# Assign input arguments to variables
DNA_FASTA=$1
PROTEIN_FASTA=$2
ALIAS=$3

# Calculate Length
awk -F '\n' '{OFS=","; getline seq; match($1, /^>(.*)/, arr); print arr[1], "length", length(seq)}' $DNA_FASTA > "${ALIAS}_length.csv"

# Calculate GC Content
awk 'BEGIN {OFS=","} /^>/ {if (seq) {gc_content = (gc_count / total_length) * 100; if (header_id != "") print header_id, "GC", gc_content}; seq = ""; gc_count = 0; total_length = 0; match($0, /^>(.*)/, arr); header_id = arr[1]; next} {seq = seq $0; total_length += length($0); gc_count += gsub(/[GCgc]/, "")} END {if (seq && header_id != "") {gc_content = (gc_count / total_length) * 100; print header_id, "GC", gc_content}}' $DNA_FASTA > "${ALIAS}_gc.csv"

# Calculate GC at 3rd Codon Position
awk 'BEGIN {OFS=","} /^>/ {if (NR > 1 && seq) {gc_3rd_content = (third_gc_count / (total_length / 3)) * 100; if (header_id != "") print header_id, "GC_3rd", gc_3rd_content}; seq = ""; third_gc_count = 0; total_length = 0; match($0, /^>(.*)/, arr); header_id = arr[1]; next} {seq = seq $0; total_length += length($0); for (i = 3; i <= length($0); i += 3) {if (substr($0, i, 1) ~ /[GCgc]/) third_gc_count++}} END {if (seq && header_id != "") {gc_3rd_content = (third_gc_count / (total_length / 3)) * 100; print header_id, "GC_3rd", gc_3rd_content}}' $DNA_FASTA > "${ALIAS}_gc_3rd.csv"

# Calculate Polar Amino Acid Content
awk 'BEGIN {OFS=","} /^>/ {if (NR > 1) {polar_content = (polar_count / total_length) * 100; print identifier, "polarAA", polar_content} identifier = substr($1, 2); total_length = 0; polar_count = 0; next} {total_length += length($0); polar_count += gsub(/[STYECQNHKR]/, "", $0)} END {if (total_length > 0) {polar_content = (polar_count / total_length) * 100; print identifier, "polarAA", polar_content}}' $PROTEIN_FASTA > "${ALIAS}_polar_aa.csv"

# Calculate Hydrophobic Amino Acid Content
awk 'BEGIN {OFS=","} /^>/ {if (NR > 1) {hydrophobic_content = (hydrophobic_count / total_length) * 100; print identifier, "hydrophobicAA", hydrophobic_content} identifier = substr($1, 2); total_length = 0; hydrophobic_count = 0; next} {total_length += length($0); hydrophobic_count += gsub(/[ACFGILMPVWY]/, "", $0)} END {if (total_length > 0) {hydrophobic_content = (hydrophobic_count / total_length) * 100; print identifier, "hydrophobicAA", hydrophobic_content}}' $PROTEIN_FASTA > "${ALIAS}_hydrophobic_aa.csv"

# Calculate Metabolic Costs
awk 'BEGIN {OFS=","} NR==FNR {costs[$1]=$2; next} /^>/ {if (NR > FNR + 1 && total_length > 0) {print identifier, "metabol", total_cost / total_length} identifier = substr($1, 2); total_length = 0; total_cost = 0; next} {for (i = 1; i <= length($0); i++) {aa = substr($0, i, 1); if (aa in costs) {total_cost += costs[aa]}} total_length += length($0)} END {if (total_length > 0) {print identifier, "metabol", total_cost / total_length}}' AA_metabolic_costs.txt $PROTEIN_FASTA > "${ALIAS}_metabol.csv"

# Calculate CAI
/stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/cai -seqall $DNA_FASTA -cfile /stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/data/CODONS/Eecoli.cut -outfile "${ALIAS}.cai"
awk '{print $2, "cai", $4}' OFS="," "${ALIAS}.cai" > "${ALIAS}_cai.csv"

# Calculate Instability
python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/calculate_instability.py $PROTEIN_FASTA "${ALIAS}_instability.csv"

# Calculate Intrinsic Structural Disorder
python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/iupred3/iupred3.py $PROTEIN_FASTA long > "${ALIAS}_iupred3_output.txt"
python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/calculate_disorder.py $PROTEIN_FASTA "${ALIAS}_iupred3_output.txt" "${ALIAS}_disorder.csv"

# Calculate Transmembrane Domain
/stor/work/Ochman/hassan/tools/phobius/phobius.pl -short $PROTEIN_FASTA > "${ALIAS}_tm.phobius"
python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/tm_to_csv.py "${ALIAS}_tm.phobius" "${ALIAS}_tm.csv"

# Calculate Transmembrane Coverage
python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/calculate_tm_coverage.py "${ALIAS}_tm.phobius" "${ALIAS}_length.csv" "${ALIAS}_tm_coverage.csv"

# Calculate Signal Peptide Presence/Absence
python /stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/sp_to_csv.py "${ALIAS}_tm.phobius" "${ALIAS}_sp.csv"

# Concatenate All CSV Files into One
cat "${ALIAS}"_*.csv > "${ALIAS}_seq_properties.csv"

echo "Pipeline completed"

