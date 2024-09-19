#!/bin/bash


TransDecoder.LongOrfs -t ./output/trinity_assembled/Trinity.fasta --gene_trans_map ./output/trinity_assembled/Trinity.fasta
.gene_trans_map -S --complete_orfs_only --output_dir ./output/transdecoder


# hmmpress uma única vez
#hmmpress /data/db/Pfam-A.hmm
hmmscan --cpu 8 --domtblout ./output/trinity_assembled/pfam.domtblout /data/db/Pfam-A.hmm ./output/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep

# makeblastdb uma única vez
#makeblastdb -dbtype prot -in ./ref/protein.fa -out ./ref/protein 
blastp -query ./output/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep  -db ./ref/protein  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > ./output/transdecoder/blastp.outfmt6

TransDecoder.Predict -t ./output/trinity_assembled/Trinity.fasta --retain_pfam_hits ./output/trinity_assembled/pfam.domtblout --retain_blastp_hits ./output/transdecoder/blastp.outfmt6 --single_best_only --no_refine_starts --output_dir ./output/transdecoder

grep '^>' ./output/transdecoder/Trinity.fasta.transdecoder.pep | cut -d ' ' -f 1 | sed 's/^>//' > peptides.txt
grep '^>' ./output/transdecoder/Trinity.fasta.transdecoder.pep | cut -d ' ' -f 1 | sed 's/^>//' | sed 's/\.p.*$//' > transcripts.txt

paste transcripts.txt peptides.txt > transcripts_peptides.txt
