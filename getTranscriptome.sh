#!/bin/bash

rm -f transcriptoma.fa

for acc in $(cut -f 1 ACCS.txt); do
        echo "Pegando FASTA para ${acc} ..."

        esearch -db nucleotide -query ${acc} | efetch \
        -format fasta >> transcriptoma.fa

done
