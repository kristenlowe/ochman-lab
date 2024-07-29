cp /stor/work/Ochman/hassan/newgene_workspace/Ecoli/ECOR2018/ECOR2018_proteins/ECOR_*_protein.faa /stor/scratch/Ochman/kristen/pangenome/ECOR_proteins/
cp /stor/work/Ochman/hassan/newgene_workspace/Ecoli/ECOR2018/ECOR2018_genomes/ECOR_*_genome.faa /stor/scratch/Ochman/kristen/pangenome/ECOR_genomes/
cp /stor/work/Ochman/hassan/newgene_workspace/Ecoli/ECOR2018/ECOR2018_cds/ECOR_*_cds.faa /stor/scratch/Ochman/kristen/pangenome/ECOR_cds/

# concatenate all cds files
# add ECOR_strain_phylogroup to each protein id
# remove pseudogenes
# convert multiline FASTA to single-line FASTA
python cat_cds_files.py
/stor/scratch/Ochman/kristen/pangenome/ECOR_cds/all_ECOR_cds.faa

# LENGTH
# process all_ECOR_cds.faa to extract ids and sequence length, add to csv file
awk -F '\n' '{OFS=","; getline seq; match($1, /protein_id=([^ ]+)/, arr); gsub(/\]/, "", arr[1]); print arr[1], length(seq), "length"}' all_ECOR_cds.faa > ECOR_sequenceproperties.csv

# GC
awk 'BEGIN {OFS=","} /^>/ {if (seq) {gc_content = (gc_count / total_length) * 100; if (header_id != "") print header_id, gc_content, "GC"}; seq = ""; gc_count = 0; total_length = 0; if (match($0, /protein_id=([^ ]+)/, arr)) {gsub(/\]/, "", arr[1]); header_id = arr[1]} else {header_id = ""} next} {seq = seq $0; total_length += length($0); gc_count += gsub(/[GCgc]/, "")} END {if (seq && header_id != "") {gc_content = (gc_count / total_length) * 100; print header_id, gc_content, "GC"}}' all_ECOR_cds.faa >> ECOR_sequenceproperties.csv

# GC 3RD
awk 'BEGIN {OFS=","} /^>/ {if (NR > 1 && seq) {gc_3rd_content = (third_gc_count / (total_length / 3)) * 100; if (header_id != "") print header_id, gc_3rd_content, "GC_3rd"}; seq = ""; third_gc_count = 0; total_length = 0; if (match($0, /protein_id=([^ ]+)/, arr)) {gsub(/\]/, "", arr[1]); header_id = arr[1]} else {header_id = ""} next} {seq = seq $0; total_length += length($0); for (i = 3; i <= length($0); i += 3) {if (substr($0, i, 1) ~ /[GCgc]/) third_gc_count++}} END {if (seq && header_id != "") {gc_3rd_content = (third_gc_count / (total_length / 3)) * 100; print header_id, gc_3rd_content, "GC_3rd"}}' all_ECOR_cds.faa >> ECOR_sequenceproperties.csv

# convert protein multiline FASTA to single-line FASTA
awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq$0} END {if (seq) print seq}' /stor/scratch/Ochman/kristen/pangenome/ECOR_proteins/ECOR_*_protein.faa > all_ECOR_proteins.faa

# POLAR AA CONTENT
awk 'BEGIN {OFS=","} /^>/ {if (NR > 1) {polar_content = (polar_count / total_length) * 100; print identifier, polar_content, "polarAA"} identifier = substr($1, 2); total_length = 0; polar_count = 0; next} {total_length += length($0); polar_count += gsub(/[STYECQNHKR]/, "", $0)} END {if (total_length > 0) {polar_content = (polar_count / total_length) * 100; print identifier, polar_content, "polarAA"}}' all_ECOR_proteins.faa >> ECOR_sequenceproperties.csv

# HYDROPHOBIC AA CONTENT
awk 'BEGIN {OFS=","} /^>/ {if (NR > 1) {hydrophobic_content = (hydrophobic_count / total_length) * 100; print identifier, hydrophobic_content, "hydrophobicAA"} identifier = substr($1, 2); total_length = 0; hydrophobic_count = 0; next} {total_length += length($0); hydrophobic_count += gsub(/[ACFGILMPVWY]/, "", $0)} END {if (total_length > 0) {hydrophobic_content = (hydrophobic_count / total_length) * 100; print identifier, hydrophobic_content, "hydrophobicAA"}}' all_ECOR_proteins.faa >> ECOR_sequenceproperties.csv

# Biosynthetic cost - Akashi and Gojobari lists
awk 'BEGIN {OFS=","} NR==FNR {costs[$1]=$2; next} /^>/ {if (NR > FNR + 1 && total_length > 0) {print identifier, total_cost / total_length, "metabol"} identifier = substr($1, 2); total_length = 0; total_cost = 0; next} {for (i = 1; i <= length($0); i++) {aa = substr($0, i, 1); if (aa in costs) {total_cost += costs[aa]}} total_length += length($0)} END {if (total_length > 0) {print identifier, total_cost / total_length, "metabol"}}' AA_metabolic_costs.txt all_ECOR_proteins.faa >> ECOR_sequenceproperties.csv

# CAI
/stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/cai -seqall all_ECOR_cds_renamed.faa -cfile /stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/data/CODONS/Eecoli.cut -outfile ECOR_sequenceproperties.cai

cut -f2,4 -d " " ECOR_sequenceproperties.cai | sed "s/ /,/g" | sed "s/$/,cai/g" >> ECOR_sequenceproperties.csv

# add gene status to ECOR_sequenceproperties.csv
python add_gene_status.py

# remove genes with incorrect number of duplicates
sort -u updated_ECOR_sequenceproperties_with_status.csv | awk -F, '{count[$1]++; lines[$1] = lines[$1] "\n" $0} END {for (gene in count) if (count[gene] == 7) printf "%s", lines[gene]}' > 7_properties.csv

# add gene frequency to ECOR_sequenceproperties.csv
python add_gene_frequency.py

# remove genes with incorrect number of duplicates
sort -u updated_ECOR_sequenceproperties_with_freq.csv | awk -F, '{count[$1]++; lines[$1] = lines[$1] "\n" $0} END {for (gene in count) if (count[gene] == 7) printf "%s", lines[gene]}' > 7_properties_freq.csv

# Hassan made updated_ECOR_sequenceproperties_with_freq.csv tidy

# Intrinsic Structural Disorder (adding to the tidy csv)
python iupred3.py all_ECOR_proteins_renamed.faa long > test_output.txt
python calculate_disorder.py
cat updated_ECOR_sequenceproperties_with_freq.tidy.csv disorder_scores.csv > 9_properties.csv

# Transmembrane Domain
/stor/work/Ochman/hassan/tools/phobius/phobius.pl -short all_ECOR_proteins_renamed.faa > ECOR_tm.phobius
python tm_to_csv.py
cat 9_properties.csv tm_scores.csv > 10_properties.csv

# Transmembrane Coverage
python calculate_tm_coverage.py
cat 10_properties.csv transmembrane_coverage.csv > 11_properties.csv

# Signal Peptide presence/absence
python sp_to_csv.py
cat 11_properties.csv sp_scores.csv > 12_properties.csv

# remove genes with incorrect number of duplicates
sort -u 12_properties.csv | awk -F, '{count[$1]++; lines[$1] = lines[$1] "\n" $0} END {for (gene in count) if (count[gene] == 12) printf "%s", lines[gene]}' > 12_properties_final.csv

# scramble proteins and run phobius on scrambled proteins
python scramble_seqs.py
/stor/work/Ochman/hassan/tools/phobius/phobius.pl -short scrambled_ECOR_proteins.faa > scrambled_ECOR_tm.phobius
python tm_to_csv.py # tm domain counts written to scrambled_tm_scores.csv
python calculate_tm_coverage.py # Coverage data written to scrambled_transmembrane_coverage.csv

