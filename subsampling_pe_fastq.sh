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
#  Based on a little BASH script by AARON
#  http://obx.cphg.virginia.edu/quinlan/?p=204

# Staring FASTQ files

FQ=( ${1} ${2})


x=0;
while [ $x != 2 ]
do
	n=$((x+1))
	if [ ! ${FQ[$x]} ]
	then
		echo "Missing pair" $n "fastq file" 1>&2
		exit
	else
		if [ ! -f ${FQ1} ]
		then
			echo "Wrong pair" $n "fastq file" 1>&2
			exit
		fi		
	fi
	x=$((x+1))
done

## How many random pairs do we want?

N=${3}

if [ ! ${N} ]
then
	echo "Missing number of random pairs" 1>&2
	exit
else
	if [ ! ${N} -gt 0 ]
	then
		echo "Number of random pairs must be greater than 0"
		exit
	fi
fi		

# The names of the random subsets
export FQ1SUBSET=${FQ[0]}.random
export FQ2SUBSET=${FQ[1]}.random

# Remove old files

rm -f ${FQ1SUBSET}
rm -f ${FQ22SUBSET}

# (1) Use paste to merge the two into a single line
# (2) Use awk to "linearize" the two mates into a single record.
#     Add a random number to the front of each line
# (3) Sort by the random number
# (4) Using head, grab the first N (randomly chosen, thanks to the awk above) records
# (5) Use awk to convert back to 2 separate FASTQ files.
paste ${FQ[0]} ${FQ[1]} | \
awk 'BEGIN{srand()}; {OFS="\t"; \
                      getline seqs; getline sep; getline quals; \
                      print rand(),$0,seqs,sep,quals}' | \
sort -k1,1 | \
head -n ${N} | \
awk 'BEGIN{FS="\t"}; {OFS="\n"; \
      print $2,$4,$6,$8 >> ENVIRON["FQ1SUBSET"]; \
      print $3,$5,$7,$9 >> ENVIRON["FQ2SUBSET"]}'
