# sequence properties for panaroo datasets

# protein sequences
cp /stor/work/Ochman/hassan/0528_geneage/061924_prokka_panaroo_pipeline/results_SeT/pan_genome_reference.prot.fa /stor/scratch/Ochman/kristen/pangenome/panaroo_seq_properties

# DNA sequences
cp /stor/work/Ochman/hassan/0528_geneage/061924_prokka_panaroo_pipeline/results_SeT/pan_genome_reference.fa /stor/scratch/Ochman/kristen/pangenome/panaroo_seq_properties

# conservation level
cp /stor/work/Ochman/hassan/0528_geneage/061924_prokka_panaroo_pipeline/results_SeT/gene_conservation_level.csv /stor/scratch/Ochman/kristen/pangenome/panaroo_seq_properties



# LENGTH
awk -F '\n' '{OFS=","; getline seq; header_id = substr($1, 2); print header_id, length(seq), "length"}' pan_genome_reference.fa > panaroo_lengths.csv

# GC
awk 'BEGIN {OFS=","} /^>/ {if (seq) {gc_content = (gc_count / total_length) * 100; if (header_id != "") print header_id, "GC", gc_content}; seq = ""; gc_count = 0; total_length = 0; header_id = substr($0, 2); next} {seq = seq $0; total_length += length($0); gc_count += gsub(/[GCgc]/, "")} END {if (seq && header_id != "") {gc_content = (gc_count / total_length) * 100; print header_id, "GC", gc_content}}' pan_genome_reference.fa > panaroo_seq_properties.csv

# Intrinsic Structural Disorder
python /stor/scratch/Ochman/kristen/pangenome/seq_properties/iupred3/iupred3.py pan_genome_reference.prot.fa long > isd_output.txt
python calculate_disorder.py
cat panaroo_seq_properties.csv disorder_scores.csv > 2_properties.csv

# Transmembrane Domain
/stor/work/Ochman/hassan/tools/phobius/phobius.pl -short pan_genome_reference.prot.fa > panaroo_tm.phobius
python tm_to_csv.py
cat 2_properties.csv tm_scores.csv > 3_properties.csv

# Transmembrane Coverage
python calculate_tm_coverage.py
cat 3_properties.csv transmembrane_coverage.csv > 4_properties.csv

# scramble proteins and run phobius on scrambled proteins
python scramble_seqs.py
/stor/work/Ochman/hassan/tools/phobius/phobius.pl -short scrambled_panaroo_proteins.faa > scrambled_panaroo_tm.phobius
python tm_to_csv.py # tm domain counts written to scrambled_tm_scores.csv
python calculate_tm_coverage.py # Coverage data written to scrambled_transmembrane_coverage.csv

