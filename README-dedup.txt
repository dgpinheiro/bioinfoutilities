#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2010  Fundação Hemocentro de Ribeirão Preto
#
#  Laboratório de Bioinformática
#  Núcleo de Bioinformática
#  BiT -  Bioinformatics Team
#  Fundação Hemocentro de Ribeirão Preto
#  Rua Tenente Catão Roxo, 2501
#  Ribeirão Preto - São Paulo
#  Brasil
#  CEP 14051-140
#  Fone: 55 16 39639300 Ramal 9603
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://lgmb.fmrp.usp.br
#

###
# Using perl scripts and sort/nsort
# This alternative uses less memory but takes long time with unix sort
###

perl merge-ends.pl -i1 input/P_1.fastq -i2 input/P_2.fastq -s '|' -o input/P.txt
# sort - ~30 minutes
sort -t '|' -k3d,3 -T /work/tmp/ input/P.txt > input/P_sorted.txt

# nsort - 30 seconds (http://www.ordinal.com/)
##nsort -t '|' -k3d,3 -T /work/tmp/ input/P.txt > input/P_sorted.txt

perl dedup-ends.pl -i input/P_sorted.txt -o input/P_sorted_unique.txt -s '|'

perl unmerge-ends.pl -i input/P_sorted_unique.txt -o1 input/P_sorted_unique_1.fastq -o2 input/P_sorted_unique_2.fastq -s '|'

###
# using PRINSEQ
# This alternative uses much more memory
###

#perl merge-ends-fastq.pl -i1 input/P_1.fastq -i2 input/P_2.fastq -o input/P.fastq
#prinseq-lite.pl -verbose -fastq  input/P.fastq -out_format 3 -derep 1 -out_good input/P_unique.fastq
#perl unmerge-ends-fastq.pl -i input/P_sunique.fastq -o1 input/P_unique_1.fastq -o2 input/P_unique_2.fastq

###
# using FASTX-Toolkit
# This alternative generates a fasta file without quality
###

#perl merge-ends-fastq.pl -i1 input/P_1.fastq -i2 input/P_2.fastq -o input/P.fastq
#fastx_collapser -Q 33 -i input/P.fastq -o input/P_unique.fasta
#perl unmerge-ends-fasta.pl -i input/P_sunique.fasta -o1 input/P_unique_1.fasta -o2 input/P_unique_2.fasta


###
# using Goby
# 
###

#java -Xmx3g -jar goby.jar -m fasta-to-compact --paired-end -d -x input/P_1.fastq -o input/P.compact-reads
#java -Xmx3g -jar goby.jar -m tally-reads -i input/P.compact-reads -o input/myfilter
#java -Xmx3g -jar goby.jar -m compact-to-fasta -t fastq -f myfilter-keep.filter -i HCC1954.compact-reads -o P_collapsed_1.fastq -p P_collapsed_.fastq
#sed 's/#\([0-9]\)\/1$/#\1\/2/' P_collapsed_.fastq > P_collapsed_2.fastq
#rm -f P_collapsed_.fastq	

