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

refgff=$1

if [ ! ${refgff} ]; then
	echo "Missing reference GFF file"
	exit
else
	if [ ! -e ${refgff} ]; then
		echo "Wrong reference GFF file (${refgff})"
		exit
	else
		if [[ ! ${refgff} =~ .g[tf]f*$ ]]; then
			echo "Wrong name of reference GFF or GTF file"
			exit
		fi
	fi
fi

basedir_out=$2

if [ ! ${basedir_out} ]; then
	echo "Missing output directory"
	exit
else
	if [ ! -d ${basedir_out} ]; then
		echo "Wrong output directory"
		exit
	fi
fi

min_perc=$3
max_perc=$4

if [ ! ${min_perc} ]; then
	min_perc=0
else
	if [[ ${min_perc} < 0 ]] || [[ ${min_perc} > 1 ]]; then
		echo "Minimum percentile (${min_perc}) must be between 0 and 1"
		exit
	fi
fi
if [ ! ${max_perc} ]; then
	max_perc=1
else
	if [[ ${max_perc} < 0 ]] || [[ ${max_perc} > 1 ]]; then
		echo "Maximum percentile (${max_perc}) must be between 0 and 1"
		exit
	fi
fi

refgff_format=`echo ${refgff} | perl -ne 'chomp; ~/\.([^\.]+)$/; print $1;'`

echo "* Gene structure statistics from current genome reference annotation (introntab.pl)"
if [ ! -e "${basedir_out}/genome_annotation_stats.txt" ]; then
	echo "   Running introntab.pl based on ${refgff_format}"

	if [ ${refgff_format} == "gtf" ]; then
		gtfbn=`basename ${refgff} .gtf`
		if [ ! -e "${basedir_out}/introntab/${gtfbn}.gff" ]; then
			echo "      Converting gtf to gff for introntab.pl"
			gtf2gff3 ${refgff} 2> /dev/null > ${basedir_out}/${gtfbn}.gff
		fi
		grep -v '^#' ${basedir_out}/${gtfbn}.gff | introntab.pl --format gff  > ${basedir_out}/genome_annotation_stats.txt
	else
		grep -v '^#' ${refgff} | introntab.pl --format ${refgff_format} > ${basedir_out}/genome_annotation_stats.txt
	fi

fi

max_intron_size=`cut -f 10 ${basedir_out}/genome_annotation_stats.txt | perl -lane 'next if (($.<=3)||($_=~/^#/)); my @isize=split(/,/, $_);  foreach my $s ( @isize ) { print $s; } ' | nsort -n | awk -v perc="${max_perc}" '{all[NR] = $0} END{print all[int(NR*perc-1)]}' `
min_intron_size=`cut -f 10 ${basedir_out}/genome_annotation_stats.txt | perl -lane 'next if (($.<=3)||($_=~/^#/)); my @isize=split(/,/, $_);  foreach my $s ( @isize ) { print $s; } ' | nsort -n | awk -v perc="${min_perc}" '{all[NR] = $0} END{print all[int(NR*perc+1)]}' `

echo "${min_intron_size},${max_intron_size}"

