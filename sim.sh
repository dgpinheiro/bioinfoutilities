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

SIMNGS_RUNFILE="/usr/local/bioinfo/simNGS/data/s_4_0099.runfile"
SIMNGS_CYCLE=151

genome_file=$1
gene_file=$2
outdir=$3

abf_dir=$4
reps=$5

number=$6
sd=$7

if [ ! ${genome_file} ]; then
	echo "Missing full genome reference file"
	exit
fi

if [ ! ${gene_file} ]; then
	echo "Missing full gene reference file"
	exit
fi

if [ ! ${outdir} ]; then
	echo "Missing output directory"
	exit
else
	if [ ! -d ${outdir} ]; then
		echo "Wrong output directory (${outdir})"
		exit
	fi
fi

if [ ! ${number} ]; then
	number=1000
fi
if [ ! ${sd} ]; then
	sd=100
fi

if [ ! ${abf_dir} ]; then
	abf_dir='.'
fi

if [ ! ${reps} ]; then
	reps=3
fi



for abf in "${abf_dir}/abundances_A.txt" "${abf_dir}/abundances_B.txt" "${abf_dir}/genes_abundances_A.txt" "${abf_dir}/genes_abundances_B.txt"; do
	if [ ! -e ${abf} ]; then
		echo "Not found ${abf}"
		exit
	fi
done

echo "Pull sequences from genomes ..."
cut -f 1 ${abf_dir}/abundances_A.txt ${abf_dir}/abundances_B.txt | pullseq -N -i ${genome_file} > ./genome_selection.fa
echo "Pull sequences from genes ..."
cut -f 1 ${abf_dir}/genes_abundances_A.txt ${abf_dir}/genes_abundances_B.txt | pullseq -N -i ${gene_file} > ./gene_selection.fa


for r in `seq 1 ${reps}`; do
	
	echo "Generate replicate A ${r} ..."
	
	generate_fragments.py -r ./genome_selection.fa -a ./abundances_A.txt -o ./tmp.frags_A${r} -t ${number} -i 300 -s ${sd}
	
	cat ./tmp.frags_A${r}.1.fasta | renameSeqs.pl -if FASTA -of FASTA -p SAMPLEA${r} -w 1000 > frags_A${r}.fa
	
	cat ./frags_A${r}.fa | simNGS \
		-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-p paired ${SIMNGS_RUNFILE} -n ${SIMNGS_CYCLE} > SAMPLEA_${r}.fastq 2> SAMPLEA_${r}.err.txt
	
	deinterleave_pairs SAMPLEA_${r}.fastq -o ${outdir}/tmp.SAMPLEA_${r}_R1.fastq ${outdir}/tmp.SAMPLEA_${r}_R2.fastq
	# o script interleave-reads.py necessita de uma nomenclatura idêntica para o par de leituras, com apenas uma diferença:
	# o /1 e /2 ou 1: e 2:
	sed 's/\/1\s\+\(.\+\)/\t1:\1/' ${outdir}/tmp.SAMPLEA_${r}_R1.fastq > ${outdir}/SAMPLEA_${r}_R1.fastq
	sed 's/\/2\s\+\(.\+\)/\t2:\1/' ${outdir}/tmp.SAMPLEA_${r}_R2.fastq > ${outdir}/SAMPLEA_${r}_R2.fastq
	rm -f ${outdir}/tmp.SAMPLEA_${r}_R1.fastq ${outdir}/tmp.SAMPLEA_${r}_R2.fastq
	
	echo "Generate replicate B ${r} ..."
	
	generate_fragments.py -r ./genome_selection.fa -a ./abundances_B.txt -o ./tmp.frags_B${r} -t ${number} -i 300 -s ${sd}
	
	cat ./tmp.frags_B${r}.1.fasta | renameSeqs.pl -if FASTA -of FASTA -p SAMPLEB${r} -w 1000 > frags_B${r}.fa
	
	cat ./frags_B${r}.fa | simNGS \
		-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-p paired ${SIMNGS_RUNFILE} -n ${SIMNGS_CYCLE} > SAMPLEB_${r}.fastq 2> SAMPLEB_${r}.err.txt
	
	# o script interleave-reads.py necessita de uma nomenclatura idêntica para o par de leituras, com apenas uma diferença:
	# o /1 e /2 ou 1: e 2:
	deinterleave_pairs SAMPLEB_${r}.fastq -o ${outdir}/tmp.SAMPLEB_${r}_R1.fastq ${outdir}/tmp.SAMPLEB_${r}_R2.fastq
	sed 's/\/1\s\+\(.\+\)/\t1:\1/' ${outdir}/tmp.SAMPLEB_${r}_R1.fastq > ${outdir}/SAMPLEB_${r}_R1.fastq
	sed 's/\/2\s\+\(.\+\)/\t2:\1/' ${outdir}/tmp.SAMPLEB_${r}_R2.fastq > ${outdir}/SAMPLEB_${r}_R2.fastq
	rm -f ${outdir}/tmp.SAMPLEB_${r}_R1.fastq ${outdir}/tmp.SAMPLEB_${r}_R2.fastq
done

