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
#  Copyright (C) 2012  Universidade de Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade de Estadual Paulista "Júlio de Mesquita Filho"
#  Faculdade de Ciências Agrárias e Veterinárias
#  Laboratório de Bioinformática (FCAV-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#

input="$1"
asmname="$2"

# Quantidade de processadores disponíveis na máquina
proc_number=$(grep -c '^processor' /proc/cpuinfo)
# Utilizar 90% dos processadores 
proc_percent=90

# Número de threads disponíveis para os processos

# Seria mais fácil determinar que o número vai de 0-1 para proc_percent e usar a linha a seguir
# bc sem o l ao dividir por 1 ele trunca o valor decimal, deixando-o somente a parte inteira
#threads=$( echo "(${proc_number} * 0.9)/1" | bc )

# Mas só para deixar um pouco mais complicado:
threads=$( echo "(${proc_number} * (${proc_percent}/100))" | bc -l | xargs printf "%.0f" )


# as linhas que iniciam com cerquilha são comentários
        
if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit 
else    
        if [ ! -d ${input} ]
        then
                echo "Wrong input directory ${input}"
                exit
        fi
fi

if [ ! ${asmname} ]
then   
	asmname="testasm"
fi

echo "Mapping data for ${input}..."

echo "   Indexing contig \"bowtie2-build\" ..."

curdir=`pwd`

cd ${input}

#mkdir -p ./mapped
#
#cd ./mapped
#
#### Criando atalho
#
#if [ ! -e final-assembly.fa ]; then
#	ln -s ../assembled/final-assembly.fa
#fi
#
#### Indexação da montagem
#
#bowtie2-build -f final-assembly.fa metagenome > bowtie2-build.out.txt 2> bowtie2-build.err.txt
#
#echo "   Alignment reads X contigs \"bowtie2\" ..."
#
#samples=()
#contig_counts=()
#
#rm -f contig_names.txt
#rm -f header.txt
#rm -f count_matrix.tab
#
#for i in ../processed/prinseq/*.scythe.cutadapt5p.filtered.prinseq.fastq; do
#	name=`basename ${i} .scythe.cutadapt5p.filtered.prinseq.fastq`
#	samples=($(printf "%s\n" ${samples[@]} `echo ${name}` | sort -u ))
#	echo "      ${name}"
#	bowtie2 --very-sensitive-local -p ${threads} -q -x metagenome -U ${i} 2> ${name}.bowtie2.err.txt | samtools view -Su - 2> /dev/null |  samtools sort - ${name}.metagenome 2> /dev/null
#	samtools index ${name}.metagenome.bam
#	if [ ! -e contig_names.txt ]; then
#		samtools idxstats ${name}.metagenome.bam | grep -v "*" | cut -f 1 > contig_names.txt
#	fi
#	samtools idxstats ${name}.metagenome.bam | grep -v "*" | cut -f 3 > ${name}.counts.txt
#	contig_counts=($(printf "%s\n" ${contig_counts[@]} `echo ${name}.counts.txt` | sort -u ))
#done
#
#echo "   Accounting contigs abundance in samples ..."
#
#header=("contig_name" ${samples[@]})
#
## IFS - Internal Field Separator (https://en.wikipedia.org/wiki/Internal_field_separator)
## define um separador distinto do espaço " ", tab "\t" e nova linha "\n" e depois retorna ao padrão
#IFS=$'\t';echo "${header[*]}" > header.txt;IFS=$' \t\n'
#paste contig_names.txt ${contig_counts[*]} | cat header.txt - > count_matrix.tab
#
#echo "   Annotating OTUs in contig data \"kraken\" ..."
#
#kraken --preload ./final-assembly.fa > ./sequences.kraken.out.txt 2> ./sequences.kraken.err.txt
#kraken-translate --mpa-format ./sequences.kraken.out.txt > ./sequences.kraken.out.labels.txt 2> ./sequences.kraken.err.labels.txt
#
#
#mergeSimpleAnnot.pl -t count_matrix.tab -k contig_name -i sequences.kraken.out.labels.txt -c 1 -a 2 -n taxonomy > count_matrix_tax_kraken.tab 2> count_matrix_tax_kraken.err.txt
#countOTUs.pl -i count_matrix_tax_kraken.tab -c taxonomy -n 100 > count_matrix_OTU.tab 2> count_matrix_OTU.err.txt
#cd ../
#
#echo "   Annotating genes in contig data \"prokka\" ..."
#
#mkdir -p ./annotated
#
#prokka ./assembled/final-assembly.fa --mincontiglen 500 --force --cpus ${threads} --outdir ./annotated --prefix ${asmname} --metagenome > ./annotated/prokka.out.txt 2> ./annotated/prokka.err.txt
#
cd ./annotated
#
#echo "   Indexing gene predictions \"bowtie2-build\" ..."
#
#mkdir -p ./align
#
#cd ./align
#
#if [ ! -e ${asmname}.ffn ]; then
#	ln -s ../${asmname}.ffn
#fi
#
## Agrupamento de genes
#
#cd-hit-est -i ${asmname}.ffn -o ${asmname}.cdhit.ffn -c 0.95 -n 3 -l 10 -aS 0.9 -d 0 -B 0 -p 1 -g 1 -T 23 -M 100000  > cd-hit-est.out.txt 2> cd-hit-est.err.txt
#
#### Indexação da montagem
#
#bowtie2-build -f ${asmname}.cdhit.ffn  ${asmname} > bowtie2-build.out.txt 2> bowtie2-build.err.txt
#
#echo "   Alignment reads X gene predictions \"bowtie2\" ..."
#
#samples=()
#gene_counts=()
#
#rm -f gene_names.txt
#rm -f header.txt
#rm -f count_matrix.tab
#
#for i in ../../processed/prinseq/*.scythe.cutadapt5p.filtered.prinseq.fastq; do
#	name=`basename ${i} .scythe.cutadapt5p.filtered.prinseq.fastq`
#	samples=($(printf "%s\n" ${samples[@]} `echo ${name}` | sort -u ))
#	echo "      ${name}"
#	bowtie2 --very-sensitive-local -p ${threads} -q -x ${asmname} -U ${i} 2> ${name}.bowtie2.err.txt | samtools view -Su - 2> /dev/null |  samtools sort - ${name}.${asmname} 2> /dev/null
#	samtools index ${name}.${asmname}.bam
#	if [ ! -e gene_names.txt ]; then
#		samtools idxstats ${name}.${asmname}.bam | grep -v "*" | cut -f 1 > gene_names.txt
#	fi
#	samtools idxstats ${name}.${asmname}.bam | grep -v "*" | cut -f 3 > ${name}.counts.txt
#	gene_counts=($(printf "%s\n" ${gene_counts[@]} `echo ${name}.counts.txt` | sort -u ))
#done
#
#echo "   Accounting gene abundance in samples ..."
#
#header=("gene_name" ${samples[@]})
#
## IFS - Internal Field Separator (https://en.wikipedia.org/wiki/Internal_field_separator)
## define um separador distinto do espaço " ", tab "\t" e nova linha "\n" e depois retorna ao padrão
#IFS=$'\t';echo "${header[*]}" > header.txt;IFS=$' \t\n'
#paste gene_names.txt ${gene_counts[*]} | cat header.txt - > count_matrix.tab
#
#mergefaAnnot.pl -f ${asmname}.cdhit.ffn -t count_matrix.tab -c gene_name -n description > ../count_matrix_desc.tab
#
#cd ../
#mkdir -p ./iprscan/
#
# Eliminando *s traduzidos de codons terminadores no interior de predições gênicas
#sed 's/\*//g' ${asmname}.faa > ${asmname}.cleaned.faa
interproscan.sh --seqtype p --tempdir /dev/shm --pathways --goterms --iprlookup --output-dir ./iprscan/ --input ${asmname}.cleaned.faa > interproscan.out.txt 2> interproscan.err.txt

mergeiprAnnot.pl -t ./count_matrix_desc.tab -i ./iprscan/${asmname}.cleaned.faa.tsv -c gene_name -a IPR > count_matrix_desc_ipr.tab

## humann2
#
#humann2 -i ./align/${asmname}.cdhit.ffn -o ./humann2_out --threads ${threads} > humann2.out.txt 2> humann2.err.txt
#


cd ${curdir}
