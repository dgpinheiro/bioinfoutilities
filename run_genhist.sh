#!/bin/bash
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
#  Copyright (C) 2012  Universidade de São Paulo
#
#  Universidade de São Paulo
#  Laboratório de Biologia do Desenvolvimento de Abelhas
#  Núcleo de Bioinformática (LBDA-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://zulu.fmrp.usp.br/bioinfo 
#

fastq_dir=$1

if [ ! ${fastq_dir} ]; then
	echo "Missing fastq directory containing *.fastqsanger files"
	exit
fi

accepted_hits_dir=$2

if [ ! ${accepted_hits_dir} ]; then
	echo "Missing accepted_hits directory containing base name directories of TopHat results for each library"
	exit
fi


cufflinks_dir=$3

if [ ! ${cufflinks_dir} ]; then
	echo "Missing cufflinks directory containing base name directories of cufflinks results for each library"
	exit
fi

output_dir=$4

if [ ! ${output_dir} ]; then
	echo "Missing output directory"
	exit
fi

div=1000000
minreads=1
minrpkm=5

gff_reference="/data/tmp/genhist/LBDAv0.4.gff"

for fastq in `ls ${fastq_dir}/*.fastqsanger`; do
	bn=`basename $fastq .fastqsanger`
	echo "Processing ${bn} ..."
	
	if [ ! -e "${accepted_hits_dir}/${bn}/accepted_hits.bam" ]; then
		echo "File not found: ${accepted_hits_dir}/${bn}/accepted_hits.bam"
		exit
	fi
	
	if [ ! -e "${cufflinks_dir}/${bn}/transcripts.gtf" ]; then
		echo "File not found: ${cufflinks_dir}/${bn}/transcripts.gtf"
		exit
	fi

	./genhist.pl -f ${fastq_dir}/${bn}.fastqsanger -a ${accepted_hits_dir}/${bn}/accepted_hits.bam -t ${cufflinks_dir}/${bn}/transcripts.gtf -g ${output_dir}/${bn} -d ${div} -m ${minreads} -k ${minrpkm} -n 5 -r ${gff_reference} &> ${output_dir}/${bn}.txt
done


