# Goal: get the distribution of proteins across the E. coli pangenome

#First, copy ECOR genomes and proteins in scratch directory: /stor/scratch/Ochman/kristen/pangenome/prot_dist

cp /stor/work/Ochman/hassan/newgene_workspace/Ecoli/ECOR2018/ECOR2018_proteins/ECOR_*_protein.faa .
cp /stor/work/Ochman/hassan/newgene_workspace/Ecoli/ECOR2018/ECOR2018_genomes/ECOR_*_genome.faa .

#Now make the diamond database

cat ECOR_*_protein.faa > ECOR_protein_db.faa
diamond makedb --in ECOR_protein_db.faa --db ECOR_protein_db -p 12

#You can now search proteins from (say) one strain in ECOR against everything in ECOR with a diamond code like this:

diamond blastp -q ECOR_1_protein.faa -d ECOR_protein_db --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --out ECOR_1_vs_ECOR.tsv -k 0 -b8 -c1 -p 12

#You can look up what those options are doing in the Diamond documentation pasted above
#This returns a protein blast results file, where the first column shows the name of the query protein (the one you're searching)
#And the second column shows the name of target protein

#By looking at which strains each protein have a hit against, you can figure out how many strains that protein is found in

#E.g.:

awk -F '\t' '($3>70&&$5>70&&$16<0.001)' ECOR_1_vs_ECOR.tsv | grep "ECOR_01_A_WP_001308243.1" | cut -f2 | cut -f1-3 -d "_" | sort -u

#The list shows the strains ECOR_01_A_WP_001308243.1 is present in
#67 such cases

#the ($3>70&&$5>70&&$16<0.001) part only selects those matches that had
#A percent identity above 70% (why this value?)
#A query cover above 70% (why this value?)
#And an e-value below 0.001 (why this value?)

#Standard blast metrics, you can look these up (You probably know about these already)

#In the (cut -f1-3 -d "_") part, I'm only getting the first three parts of the target protein name
#Which tells me what ECOR strain they're from

#So the task is to make a presence/absence table for all of the proteins in ECOR01

#Which we can then parse to figure out which ones are core/cloud/shell/lineage-specific/all that jazz

#Eventually we have to figure this out for all proteins in all 70 ECOR genomes