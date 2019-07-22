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


rm -f mock_genes.fa 
rm -f genes_description.txt

if [ ! -e "list_for_genes.csv" ]; then
	echo "Missing list_for_genes.csv"
	exit
fi

IFS=$'\n'; for i in $(cut -f 1 list_for_genes.csv | grep -v -P '^\s*$' | grep ':' | cut -d ':' -f 1 | sort -u); do 
	acc=`echo ${i} | cut -f 1`; 
	echo ${acc}

	if [ ! -e ${acc}.gb ]; then
		esearch -db nucleotide -query ${acc} | efetch -format gbwithparts > ${acc}.gb; 
	fi

	id=()
	IFS=$'\n'; for j in $(grep "${acc}" list_for_genes.csv | cut -f 1 | grep -v -P '^\s*$' | grep ':' | cut -d ':' -f 3 | sort -u); do 
		id=(${id[@]} ${j})
		symbol=`esearch -db gene -query ${j} | efetch -format tabular | cut -f 6 | tail -1`
		if [ ! -e ${acc}_${j}_${symbol}_500.fa ]; then
			echo "Get ${acc} ${j} ${symbol} ..."
			gbGeneSearch.pl -g ${acc}.gb -i $j -e 500 > ${acc}_${j}_${symbol}_500.fa
		fi
	done
	grep -h '^>' ${acc}_${j}_${symbol}_500.fa | sed 's/^>//' >> genes_description.txt
	
	cat ${acc}_${j}_${symbol}_500.fa | perl -lane 'if ($_=~/^>(\S+)/) { my $acc=$1; $acc=~s/\.\d+$//; print ">",$acc; } else { print $_; }' >> mock_genes.fa

done

rm -f genes_abundances_A.txt
IFS=$'\n'; for i in $(grep -v -P '^\s+' list_for_genes.csv); do acc=`echo ${i} | cut -f 1`; percent=`echo ${i} | cut -f 4`; echo -e ${acc}"\t"$( echo "scale=5; $percent/100" | bc -l | sed 's/^\./0\./') >> genes_abundances_A.txt; done

perl -F"\t" -lane ' INIT { my $sum=0;} $sum+=$F[1]; END{ print "SUM A: ",sprintf("%.2f",$sum); }' genes_abundances_A.txt

rm -f genes_abundances_B.txt
IFS=$'\n'; for i in $(grep -v -P '^\s+' list_for_genes.csv); do acc=`echo ${i} | cut -f 1`; percent=`echo ${i} | cut -f 7`; echo -e ${acc}"\t"$( echo "scale=5; $percent/100" | bc -l | sed 's/^\./0\./') >> genes_abundances_B.txt; done

perl -F"\t" -lane ' INIT { my $sum=0;} $sum+=$F[1]; END{ print "SUM B: ",sprintf("%.2f",$sum); }' genes_abundances_B.txt
