#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

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

# num_threads - inteiro com o número de processadores para serem utilizados nos parâmetros dos programas que serão executados a seguir
num_threads=2

###
# Arquivos e diretórios de entrada (input)
#

# scythe_adapters_path - caminho para o arquivo contendo os adaptadores para o programa scythe
scythe_adapters_p1_path="refs/truseq_3p_adapt_R1.fa"
scythe_adapters_p2_path="refs/truseq_3p_adapt_R2.fa"

# refgff - caminho para o arquivo contendo as coordenadas gênicas
# Dbuz
#refgff="refs/D.buzzatii_annotations.gff3"
# Dmoj
#refgff="refs/dmoj-all-filtered-r1.3.gff"
# Dmel
refgff="refs/dmel-all-r6.05.gtf"

# refbt2idx - caminho para a base com o índice do genoma para o bowtie2
# Dbuz
#refbt2idx="refs/Dbuz-index/D.buzzatii_Freeze1_Scaffolds"
# Dmoj
#refbt2idx="refs/Dmoj-index/dmoj-all-chromosome-r1.3"
# Dmel
refbt2idx="refs/Dmel-index/dmel-all-chromosome-r6.05"

# regenomefa - caminho para o arquivo com as sequências do genoma no formato multifasta
# Dbuz
#refgenomefa="refs/Dbuz_Freeze1.fa"
# Dmoj
#refgenomefa="refs/dmoj-all-chromosome-r1.3.fasta"
# Dmel
refgenomefa="refs/dmel-all-chromosome-r6.05.fasta"

# tophat_transcriptome_index - bowtie2 transcriptome index relacionado ao genes.gtf selecionado
# Dbuz
#tophat_transcriptome_index="./refs/transcriptome_index/Dbuz-transcriptome_index/"
# Dmouj
#tophat_transcriptome_index="./refs/transcriptome_index/Dmoj-transcriptome_index/"
# Dmel
tophat_transcriptome_index="./refs/transcriptome_index/Dmel-transcriptome_index/"

# cutadapt_5p_adapter - adaptador 3' para ser usado com cutadapt
cutadapt_adapters_p1_path="refs/truseq_5p_adapt_R1.fa"
cutadapt_adapters_p2_path="refs/truseq_5p_adapt_R2.fa"

###
# Arquivos e diretórios de saída (output) 
#

basedir_out="."

# fastx_out - caminho para o diretório de saída do fastx_trimmer
fastx_out="${basedir_out}/raw2"

# cutadapt_out - caminho para o diretório de saída do cutadapt
cutadapt_out="${basedir_out}/processed/cutadapt"

# scythe_out - caminho para o diretório de saída do scythe
scythe_out="${basedir_out}/processed/scythe"
	
# tophat_out - caminho para o diretório de saída do tophat
tophat_out="${basedir_out}/tophat"

# prinseq_out - caminho para o diretório de saída do prinseq
prinseq_out="${basedir_out}/processed/prinseq"

# cufflinks_out - caminho para o diretório de saída do cufflinks
cufflinks_out="${basedir_out}/cufflinks"

# cuffmerge_out - caminho para o diretório de saída do cuffmerge
cuffmerge_out="${basedir_out}/cuffmerge"

# cuffquant_out - caminho para o diretório de saída do cuffquant
cuffquant_out="${basedir_out}/cuffquant"

# cuffnorm_out - caminho para o diretório de saída do cuffnorm
cuffnorm_out="${basedir_out}/cuffnorm"

# cuffdiff_out - caminho para o diretório de saída do cuffdiff
cuffdiff_out="${basedir_out}/cuffdiff"


echo "* Creating directories"

# Criando diretórios para as saídas dos programas que serão utilizados a seguir
mkdir -p ${fastx_out}
mkdir -p ${tophat_out}
mkdir -p ${cufflinks_out}
mkdir -p ${cuffmerge_out}
mkdir -p ${cuffquant_out}
mkdir -p ${cuffnorm_out}
mkdir -p ${cuffdiff_out}
mkdir -p ${scythe_out}
mkdir -p ${cutadapt_out}
mkdir -p ${prinseq_out}


# biosamps - array com os nomes dos arquivos .fastq de cada réplica (técnica/biológica)
biosamps=()
# biogroups - array com os nomes dos grupos de amostras obtidos a partir do nome do arquivo fastq inicial
biogroups=()

# menor tamanho de sequência
min=1000

echo "* Identifying sequence size: "

# percorrendo múltiplos valores retornados por um comando no shell (o comando irá listar todos os arquivos com extensão .fastq do diretório ${input})
for fastq in `ls ${input}/*.fastq`; do

	# obtendo nome do arquivo 
	fastqbn=`basename ${fastq} .fastq`;
	

	# carregando valores de biosamps com o nome de cada réplica biológica. Por ex.: FLN_B1 obtido de FLN_B1_T1
	# Exemplo na linha de comando para remover valores duplicados * pois é possível que haja mais de 1 réplica técnica: 
	# LIST=("A1" "A1" "A1" "A2" "C2" "C1" "B1" "B2") && echo $(printf "%s\n" ${LIST[@]} | sort -u)
	biosamps=($( printf "%s\n" ${biosamps[@]} `echo ${fastqbn}  | cut -d '_' -f 1,2` | sort -u ))
	
	# carregando valores de biogroups com o nome de cada grupo de amostras. Por ex.: FLN obtido de FLN_B1_T1
	biogroups=($(printf "%s\n" ${biogroups[@]} `echo ${fastqbn}  | cut -d '_' -f 1` | sort -u ))
	
	# obtendo a primeira sequência do arquivo ${fastq}
	firstSeq=`head -2 ${fastq} |  tail -1`

	echo "   ${fastqbn}: ${#firstSeq} bases"

	# testando se o tamanho de ${firstSeq} é menor que mínimo
	if [ ${#firstSeq} -le ${min} ]
	then
		min=${#firstSeq}
	fi
done

echo "   Minimum sequence size: ${min}"

for bs in ${biosamps[@]}; do
	# Concatena os dados das réplicas técnicas em uma única réplica biológica e podas as sequências para deixá-las do mesmo tamanho
	# caso ainda não tenham sido concatenados e podados (não existe -e)
	if [ ! -e "${fastx_out}/${bs}_1.fastq" ]  
	then
		
		techreps=(`ls ${input}/${bs}_T*_P1.fastq`)
		if ((${#techreps[@]} > 1))
		then
			echo "* Trimming ${bs} by minimum length (${min}) and merging ${#techreps[@]} technical replicates of ${bs} [P1] ..."
			cat ${input}/${bs}_T*_P1.fastq | fastx_trimmer -l ${min} -Q 33 > ${fastx_out}/${bs}_1.fastq
		else
			echo "* Trimming ${bs} by minimum length (${min}) [P1] ..."
			cat ${input}/${bs}_T*_P1.fastq | fastx_trimmer -l ${min} -Q 33 > ${fastx_out}/${bs}_1.fastq
		fi
	fi
	if [ ! -e "${fastx_out}/${bs}_2.fastq" ]  
	then
		
		techreps=(`ls ${input}/${bs}_T*_P2.fastq`)
		if ((${#techreps[@]} > 1))
		then
			echo "* Trimming ${bs} by minimum length (${min}) and merging ${#techreps[@]} technical replicates of ${bs} [P2] ..."
			cat ${input}/${bs}_T*_P2.fastq | fastx_trimmer -l ${min} -Q 33 > ${fastx_out}/${bs}_2.fastq
		else
			echo "* Trimming ${bs} by minimum length (${min}) [P2] ..."
			cat ${input}/${bs}_T*_P2.fastq | fastx_trimmer -l ${min} -Q 33 > ${fastx_out}/${bs}_2.fastq
		fi
	fi
done


# removendo arquivo com o caminho de cada transcriptoma (.gtf)
rm -f ${cuffmerge_out}/assembly_GTF_list.txt

echo "* Checking pre-processing:"

for fastq in `ls ${fastx_out}/*_1.fastq`; do
	
	# fqname - string com o nome do arquivo fastq

	fqname=`basename ${fastq} _1.fastq`

	echo "   ${fqname} [P1] ..."	

        if [ ! -e "${scythe_out}/${fqname}_1.scythe.fastq" ]
        then   

                mkdir -p ${scythe_out}
                
		echo "      Adapter trimming - step 1/2 (scythe) [P1] ..."

                scythe -p 0.1 -m ${scythe_out}/${fqname}_1.scythe.matches.txt -o ${scythe_out}/${fqname}_1.scythe.fastq -n 5 -M 15  -q sanger  -a ${scythe_adapters_p1_path} -i ${fastq} 1> ${scythe_out}/${fqname}_1.scythe.log.txt 2> ${scythe_out}/${fqname}_1.scythe.err.txt

        fi


        if [ ! -e "${cutadapt_out}/${fqname}_1.scythe.cutadapt.fastq" ]
        then
                
		mkdir -p ${cutadapt_out}

                echo "      Adapter trimming - step 2/2 (cutadapt) [P1] ..."
				
		cutadapt_5p_adapter_p1=`perl -ane 'chomp; if ($_=~/^>(\S+)/) { } else { print "--front=\"$_\""," "; }' ${cutadapt_adapters_p1_path}`;
                cutadapt --format=fastq ${cutadapt_5p_adapter_p1} --error-rate=0.1 --times=2 --overlap=5 --minimum-length=15 --output=${cutadapt_out}/${fqname}_1.scythe.cutadapt.fastq ${scythe_out}/${fqname}_1.scythe.fastq > ${cutadapt_out}/${fqname}_1.scythe.cutadapt.report.txt 2> ${cutadapt_out}/${fqname}_1.scythe.cutadapt.log.txt

        fi
done

for fastq in `ls ${fastx_out}/*_2.fastq`; do
	
	# fqname - string com o nome do arquivo fastq

	fqname=`basename ${fastq} _2.fastq`

	echo "   ${fqname} [P2] ..."	

        if [ ! -e "${scythe_out}/${fqname}_2.scythe.fastq" ]
        then   

                mkdir -p ${scythe_out}
                
		echo "      Adapter trimming - step 1/2 (scythe) [P2] ..."

                scythe -p 0.1 -m ${scythe_out}/${fqname}_2.scythe.matches.txt -o ${scythe_out}/${fqname}_2.scythe.fastq -n 5 -M 15  -q sanger  -a ${scythe_adapters_p2_path} -i ${fastq} 1> ${scythe_out}/${fqname}_2.scythe.log.txt 2> ${scythe_out}/${fqname}_2.scythe.err.txt

        fi


        if [ ! -e "${cutadapt_out}/${fqname}_2.scythe.cutadapt.fastq" ]
        then
                
		mkdir -p ${cutadapt_out}

                echo "      Adapter trimming - step 2/2 (cutadapt) [P2] ..."
		
		cutadapt_5p_adapter_p2=`perl -ane 'chomp; if ($_=~/^>(\S+)/) { } else { print "--front=\"$_\""," "; }' ${cutadapt_adapters_p2_path}`;

                cutadapt --format=fastq ${cutadapt_5p_adapter_p2} --error-rate=0.1 --times=2 --overlap=5 --minimum-length=15 --output=${cutadapt_out}/${fqname}_2.scythe.cutadapt.fastq ${scythe_out}/${fqname}_2.scythe.fastq > ${cutadapt_out}/${fqname}_2.scythe.cutadapt.report.txt 2> ${cutadapt_out}/${fqname}_2.scythe.cutadapt.log.txt

        fi
done


for fastq in `ls ${fastx_out}/*_1.fastq`; do

	# fqname - string com o nome do arquivo fastq

	fqname=`basename ${fastq} _1.fastq`

	echo "   ${fqname} ..."	

        if [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ] && [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq" ];
        then

                echo "      Quality and poly-A/T trimming (prinseq) ..."

                mkdir -p ${prinseq_out}
		
		# -graph_data ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.gd -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn
                # trimming reads by quality value
                prinseq-lite.pl -fastq ${cutadapt_out}/${fqname}_1.scythe.cutadapt.fastq -fastq2 ${cutadapt_out}/${fqname}_2.scythe.cutadapt.fastq -out_format 3 -out_good ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq -out_bad ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq-bad -qual_noscale -no_qual_header -min_len 15 -ns_max_p 80 -noniupac -trim_tail_left 5 -trim_tail_right 5 -trim_qual_right 25 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 3 -trim_qual_step 1 &> ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.log.txt

        fi
	
	
	# Execução do tophat2 para obter o alinhamento de cada réplica biológica em relação ao genoma
	# caso ainda não tenha sido executado
	tophat_sample_out="${tophat_out}/${fqname}"
	if [ !  -e "${tophat_sample_out}/accepted_hits.bam" ]
	then

		mkdir -p ${tophat_sample_out}
		
		echo "      Aligning ${fqname} against genome (tophat2) ..."

		tophat2 --num-threads=${num_threads} --library-type=fr-unstranded --prefilter-multihits --coverage-search --max-multihits 10 --b2-very-sensitive --min-coverage-intron 20 --max-coverage-intron 20000 --microexon-search  --max-deletion-length 3 --max-insertion-length 3 --read-mismatches 2 --min-anchor 8 --splice-mismatches 0 --min-intron-length 20 --max-intron-length 200000 --min-segment-intron 20 --max-segment-intron 200000 --segment-mismatches=2 --segment-length=20 --output-dir ${tophat_sample_out}/ --transcriptome-index ${tophat_transcriptome_index} -G ${refgff} ${refbt2idx} ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq,${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq,${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq &> ${tophat_sample_out}/${fqname}.scythe.cutadapt.prinseq.tophat.log.txt

	fi
	
	# Execução do cufflinks para obter a montagem para cada réplica biológica
	# caso ainda não tenha sido executado
	cufflinks_sample_out="${cufflinks_out}/${fqname}"
	if [ !  -e "${cufflinks_sample_out}/transcripts.gtf" ]
	then

		mkdir -p ${cufflinks_sample_out}

		echo "      Assembling ${fqname} (cufflinks) ..."

		cufflinks --num-threads ${num_threads} --library-type fr-unstranded --multi-read-correct --GTF-guide ${refgff} --frag-bias-correct ${refgenomefa} --total-hits-norm --output-dir ${cufflinks_sample_out} ${tophat_sample_out}/accepted_hits.bam --max-bundle-frags 2000000 &> ${cufflinks_sample_out}/${fqname}.cufflinks.log.txt

	fi
	
	# gravando no arquivo ${cuffmerge_out}/assembly_GTF_list.txt o caminho de cada transcriptoma para cada réplica biológica. Esse arquivo é requerido para a execução do cuffmerge
	echo "${cufflinks_sample_out}/transcripts.gtf" >> ${cuffmerge_out}/assembly_GTF_list.txt
done

# Execução do cuffmerge para fusionar os transcriptomas gerados em um transcriptoma referência
# Será executado caso ainda não tenha o arquivo de saída esperado (merged.gtf)

if [ ! -e "${cuffmerge_out}/merged.gtf" ]
then
	rowcount=`wc -l ${cuffmerge_out}/assembly_GTF_list.txt | sed 's/ .*//'`;
	echo "* Merging ${rowcount} transcriptomes (cuffmerge) ..."

	cuffmerge --num-threads ${num_threads} --ref-sequence ${refgenomefa} -o ${cuffmerge_out} --ref-gtf ${refgff} ${cuffmerge_out}/assembly_GTF_list.txt &> ${cuffmerge_out}/cuffmerge.log.txt

fi

echo "Creating quantification files (.cxb):"
for bs in ${biosamps[@]}; do
	
       	echo "	${bs} ..."
	
	mkdir -p ${cuffquant_out}/${bs}
	
	if [ ! -e "${cuffquant_out}/${bs}/abundances.cxb" ]
	then
		cuffquant --library-type fr-unstranded --multi-read-correct --frag-bias-correct ${refgenomefa} ${cuffmerge_out}/merged.gtf ${tophat_out}/${bs}/accepted_hits.bam -o ${cuffquant_out}/${bs} &> ${cuffquant_out}/${bs}.log.txt
	fi
	
done;

biogroupscxb=()

echo "Checking cxb files:"
for bg in ${biogroups[@]}; do
       	echo "	${bg} ..."
	cxblist=(`find ${cuffquant_out}/ -name abundances.cxb | grep "${bg}" | sort`)
	# aqui temos que usar o asterisco ("*") para expandir o conteúdo de cxblist pois queremos apenas uma palavra, separada pela variável especial IFS e não múltiplas palavras como quando usamos o arroba ("@")
	
	biogroupscxb=(${biogroupscxb[@]} $(IFS=, ; echo "${cxblist[*]}"))
done;

# Execução do cuffnorm para obter as estimativas de expressão gênica
# caso não tenha sido executado
if [ !  -e "${cuffnorm_out}/run.info" ]
then
	echo "* Estimating Gene Expression (cuffnorm) from groups:" "${biogroups[*]}" 
	cuffnorm --num-threads ${num_threads} --library-type fr-unstranded --labels $(IFS=, ; echo "${biogroups[*]}")  -o ${cuffnorm_out} --library-norm-method classic-fpkm ${cuffmerge_out}/merged.gtf ${biogroupscxb[@]} &> ${cuffnorm_out}/cuffnorm.log
fi

# Execução do cuffdiff para obter as estimativas de expressão gênica diferencial
# caso não tenha sido executado
if [ !  -e "${cuffdiff_out}/run.info" ]
then
	echo "* Estimating Differential Expression (cuffdiff) among groups:" "${biogroups[*]}" 
	cuffdiff --num-threads ${num_threads} --library-type fr-unstranded --multi-read-correct --frag-bias-correct ${refgenomefa} --upper-quartile-norm --compatible-hits-norm --output-dir ${cuffdiff_out} ${cuffmerge_out}/merged.gtf ${biogroupscxb[@]} --FDR  0.05 --labels $(IFS=, ; echo "${biogroups[*]}") --max-bundle-frags 999999999 &> ${cuffdiff_out}/cuffdiff.log
fi

# Documentação do bash
#	http://www.gnu.org/software/bash/manual/bashref.html#Arrays
# Discussão sobre * ou @ ao expandir valores do array
#	Suppose we have perls=(perl-one perl-two) and use ls "${perls[*]}" -- that's equivalent to ls "perl-one perl-two", which will look for single file named perl-one perl-two, which is probably not what you wanted. ls "${perls[@]}" is equivalent to ls "perl-one" "perl-two", which is much more likely to do something useful.
# 	http://www.coderexception.com/0166bbz36XXQJQQX/a-confusion-about-array-versus-array-in-the-context-of-a-bash-completion
