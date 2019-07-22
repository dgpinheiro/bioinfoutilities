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

if [ ! -e "./list.csv" ]; then
	echo "Missing list.csv"
	exit
fi

rm -f ./mock.fa 
rm -f ./description.txt

IFS=$'\n'; for i in $(grep -v -P '^\s' list.csv); do 
	acc=`echo ${i} | cut -f 1`; 
	spc=`echo ${i} | cut -f 2 | sed 's/ /_/g'`; 
	if [ ! -e ${spc}_${acc}.fa ]; then
		echo "${acc}-${spc}"; 
		esearch -db nucleotide -query ${acc} | efetch -format fasta > ${spc}_${acc}.fa; 
	fi		
	grep -h '^>' ${spc}_${acc}.fa | sed 's/^>//' >> ./description.txt
done 
 
IFS=$'\n'; for i in $(grep -v -P '^\s' list.csv); do 
	acc=`echo ${i} | cut -f 1`; 
	spc=`echo ${i} | cut -f 2 | sed 's/ /_/g'`; 
	echo "${acc}-${spc}"; 
	cat ${spc}_${acc}.fa | perl -lane 'if ($_=~/^>(\S+)/) { my $acc=$1; $acc=~s/\.\d+$//; print ">",$acc; } else { print $_; }' >> mock.fa; 
done

rm -f abundances_A.txt
IFS=$'\n'; for i in $(grep -v -P '^\s' list.csv); do acc=`echo ${i} | cut -f 1`; percent=`echo ${i} | cut -f 4`; echo -e ${acc}"\t"$( echo "scale=5; $percent/100" | bc -l | sed 's/^\./0\./') >> abundances_A.txt; done

perl -F"\t" -lane ' INIT { my $sum=0;} $sum+=$F[1]; END{ print "SUM A: ",sprintf("%.2f",$sum); }' abundances_A.txt

rm -f abundances_B.txt
IFS=$'\n'; for i in $(grep -v -P '^\s' list.csv); do acc=`echo ${i} | cut -f 1`; percent=`echo ${i} | cut -f 6`; echo -e ${acc}"\t"$( echo "scale=5; $percent/100" | bc -l | sed 's/^\./0\./') >> abundances_B.txt; done

perl -F"\t" -lane ' INIT { my $sum=0;} $sum+=$F[1]; END{ print "SUM B: ",sprintf("%.2f",$sum); }' abundances_B.txt
