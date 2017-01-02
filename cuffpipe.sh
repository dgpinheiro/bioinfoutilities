#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

aligner=$2

if [ ! ${aligner} ]
then
	aligner="tophat2"
else
	if [ ! ${aligner} == "tophat2" ] && [ ! ${aligner} == "star" ] && [ ! ${aligner} == "bbmap" ]
	then
		echo "Wrong aligner (${aligner})"
		exit
	fi
fi


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

####
### Loading config file
####

configfile='./cuffpipe.cfg'
configfile_secured=$(mktemp)

# checagem se o arquivo contem algo indesejavel (que não segue o padrao esperado para arquivos de configuracao chave=valor)
if egrep -q -v '^#|^[^ ]*=[^;]*' "${configfile}"; then
  echo "Config file is unclean, cleaning it..." >&2
  # filtragem do arquivo original para um novo arquivo
  egrep '^#|^[^ ]*=[^;&]*'  "${configfile}" > "${configfile_secured}"
  configfile="${configfile_secured}"
fi

# carrega o arquivo original ou sua versao filtrada
source "${configfile}"

# remove arquivo de configuração temporario
rm ${configfile_secured}

###
# Arquivos e diretórios de entrada (input)
#

# scythe_adapters_path - caminho para o arquivo contendo os adaptadores para o programa scythe
scythe_adapters_r1_path="/share/tank/db/adap_barc_primer/truseq_illumina_3p_adapt_R1.fa"
scythe_adapters_r2_path="/share/tank/db/adap_barc_primer/truseq_illumina_3p_adapt_R2.fa"

# refstaridx - caminho para a base com o índice do genoma para o STAR
refstaridx="./refs/STAR"

# refbbmap - caminho para a base com o índice do genoma para o bbmap
refbbmapidx="./refs/bbmap"

# num_threads - inteiro com o número de processadores para serem utilizados nos parâmetros dos programas que serão executados a seguir
# Checando se num_threads foi definida
if [ ! ${num_threads} ]; then
	num_threads=1
	echo "Using default number of threads (${num_threads})";
fi

# refgff - caminho para o arquivo contendo as coordenadas gênicas
# Checando se refgff foi definida
if [ ! ${refgff} ]; then
	echo "Reference gff file (refgff) was not defined";
	exit;
fi

# refbt2idx - caminho para a base com o índice do genoma para o bowtie2
# Checando se refbt2idx foi definida
if [ ! ${refbt2idx} ]; then
	echo "Reference Bowtie2 index (refbt2idx) was not defined";
	exit;
fi

# refgenomefa - caminho para o arquivo com as sequências do genoma no formato multifasta
# Checando se refgenomefa foi definida
if [ ! ${refgenomefa} ]; then
	echo "Reference genome fasta file (refgenomefa) was not defined";
	exit;
fi

# tophat_transcriptome_index - bowtie2 transcriptome index relacionado ao genes.gtf selecionado
# Checando se tophat_transcriptome_index foi definida
if [ ! ${num_threads} ]; then
	tophat_transcriptome_index="./refs/transcriptome_index/transcriptome"
	echo "Using default tophat_transcriptome_index (${tophat_transcriptome_index})";
fi

mkdir -p $(dirname ${tophat_transcriptome_index})

# cutadapt_5p_adapters_r1_path - adaptador 3' para ser usado com cutadapt read R1
# cutadapt_5p_adapters_r2_path - adaptador 3' para ser usado com cutadapt read R2
cutadapt_5p_adapters_r1_path="/share/tank/db/adap_barc_primer/truseq_illumina_5p_adapt_R1.fa"
cutadapt_5p_adapters_r2_path="/share/tank/db/adap_barc_primer/truseq_illumina_5p_adapt_R2.fa"

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

# align0_out - caminho para o diretório de saída do alinhamento (alignment 0)
align0_out="${basedir_out}/align0"
	
# align1_out - caminho para o diretório de saída do alinhamento (alignment 1)
align1_out="${basedir_out}/align1"

# align2_out - caminho para o diretório de saída do alinhamento (alignment 2)
align2_out="${basedir_out}/align2"

# align3_out - caminho para o diretório de saída do alinhamento (alignment 3)
align3_out="${basedir_out}/align3"

# finalalign_out - caminho para o diretório de saída do alinhamento final
finalalign_out="${basedir_out}/finalalign"

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


refgff_format=`echo ${refgff} | perl -ne 'chomp; ~/\.([^\.]+)$/; print $1;'`

mkdir -p ${basedir_out}/introntab

echo "* Gene structure statistics from current genome reference annotation (introntab.pl)"
if [ ! -e "${basedir_out}/introntab/genome_annotation_stats.txt" ]; then
	echo "   Running introntab.pl"
	introntab.pl --format ${refgff_format} < ${refgff} > ${basedir_out}/introntab/genome_annotation_stats.txt
fi

max_intron_size=`cut -f 10 ${basedir_out}/introntab/genome_annotation_stats.txt | perl -lane 'next if (($.<=3)||($_=~/^#/)); my @isize=split(/,/, $_);  foreach my $s ( @isize ) { print $s; } ' | nsort -n | awk '{all[NR] = $0} END{print all[int(NR*0.75 - 0.5)]}' `
min_intron_size=`cut -f 10 ${basedir_out}/introntab/genome_annotation_stats.txt | perl -lane 'next if (($.<=3)||($_=~/^#/)); my @isize=split(/,/, $_);  foreach my $s ( @isize ) { print $s; } ' | nsort -n | awk '{all[NR] = $0} END{print all[int(NR*0.25 - 0.5)]}' `

if [ ! ${max_intron_size} ] || [ ! ${min_intron_size} ] ; then
	max_intron_size=`grep 'Genome_size:' /tmp/out | perl -lane '$_=~/(\d+)/; print $1;' `;
	min_intron_size=`grep 'Genome_size:' /tmp/out | perl -lane '$_=~/(\d+)/; print $1;' `;
fi

echo "   Max intron size: ${max_intron_size}"
echo "   Min intron size: ${min_intron_size}"

echo "* Creating directories"

# Criando diretórios para as saídas dos programas que serão utilizados a seguir
mkdir -p ${fastx_out}
mkdir -p ${align0_out}
mkdir -p ${align1_out}
mkdir -p ${align2_out}
mkdir -p ${align3_out}
mkdir -p ${finalalign_out}
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
for fastq in `ls ${input}/*R?.fastq`; do

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
	if [ ! -e "${fastx_out}/${bs}_R1.fastq" ]  
	then
		
		techreps=(`ls ${input}/${bs}_T*_R1.fastq`)
		if ((${#techreps[@]} > 1))
		then
			echo "* Trimming ${bs} by minimum length (${min}) and merging ${#techreps[@]} technical replicates of ${bs} [R1] ..."
		else
			echo "* Trimming ${bs} by minimum length (${min}) [R1] ..."
		fi
		cat ${input}/${bs}_T*_R1.fastq | fastx_trimmer -l ${min} -Q 33 > ${fastx_out}/${bs}_R1.fastq
	fi
	if [ ! -e "${fastx_out}/${bs}_R2.fastq" ]  
	then
		# Erro esta sendo redirecionado para /dev/null pois pode não haver arquivos com o padrão desejado, o que nesse caso seriam 
		# bibliotecas single-end, ou seja, somente possuem "*_T*_R1.fastq"
		techreps=(`ls ${input}/${bs}_T*_R2.fastq 2>/dev/null`)
		if ((${#techreps[@]} > 0))
		then
			if ((${#techreps[@]} > 1))
			then
				echo "* Trimming ${bs} by minimum length (${min}) and merging ${#techreps[@]} technical replicates of ${bs} [R2] ..."
			else
				echo "* Trimming ${bs} by minimum length (${min}) [R2] ..."
			fi
			cat ${input}/${bs}_T*_R2.fastq | fastx_trimmer -l ${min} -Q 33 > ${fastx_out}/${bs}_R2.fastq
		fi			
	fi
done


function tophat_alignment_step {
	local samplename=$1
	local step=$2
	local alignX_sample_out=$3
	local alignL_sample_out=$4
	local fastqin=$5
	local fastq2in=$6
	
	local nonovel=" --no-novel-juncs "
	local nomixed=" --no-mixed "
	if [ ${step} == 3 ]
	then
		nomixed=""
		nonovel=""
	fi
	
	local cmdtophat2="tophat2 --keep-tmp ${nomixed} ${nonovel} --num-threads=${num_threads} --library-type=fr-unstranded --prefilter-multihits --coverage-search --max-multihits 10 --b2-very-sensitive --min-coverage-intron ${min_intron_size} --max-coverage-intron ${max_intron_size} --microexon-search  --max-deletion-length 3 --max-insertion-length 3 --read-mismatches 2 --min-anchor $((min/3)) --splice-mismatches 0 --min-intron-length ${min_intron_size} --max-intron-length ${max_intron_size} --min-segment-intron ${min_intron_size} --max-segment-intron ${max_intron_size} --segment-mismatches=2 --segment-length=$((min/2)) --output-dir ${alignX_sample_out}/ --transcriptome-index ${tophat_transcriptome_index} -G ${refgff}"
	
	if [ -s "${alignL_sample_out}/all_junctions.junc" ]
	then   
		cmdtophat2="${cmdtophat2} --raw-juncs ${alignL_sample_out}/all_junctions.junc "
	fi
	
	cmdtophat2="${cmdtophat2} ${refbt2idx} ${fastqin}"

	local unmapped_second_end_fastq=""
	if [ -e ${fastq2in} ] 
	then
		cmdtophat2="${cmdtophat2} ${fastq2in}"
		unmapped_second_end_fastq="SECOND_END_FASTQ=${alignX_sample_out}/${samplename}.unmapped_R2.fastq"
	fi

	cmdtophat2="${cmdtophat2} &> ${alignX_sample_out}/${samplename}.tophat.log.txt"
	
	eval ${cmdtophat2}
			
	if [ -e "${alignX_sample_out}/junctions.bed" ]
	then
       	        bed_to_juncs < ${alignX_sample_out}/junctions.bed > ${alignX_sample_out}/junctions.junc 2> /dev/null
	else
		touch ${alignX_sample_out}/junctions.junc				
	fi
	
	
	
	cat ${alignX_sample_out}/junctions.junc | nsort -u > ${alignX_sample_out}/all_junctions.junc
			
	local picard_input="INPUT=${alignX_sample_out}/unmapped.bam"
	if [ -e ${fastq2in} ]
	then
		fixMissingMate.pl ${fastqin} ${alignX_sample_out}/tmp/left_multimapped.bam ${alignX_sample_out}/tmp/right_multimapped.bam | samtools view -b -S - > ${alignX_sample_out}/unmapped_multi.bam
		
		samtools view -H ${alignX_sample_out}/tmp/left_multimapped.bam > ${alignX_sample_out}/tmp/header.sam
		mv ${alignX_sample_out}/unmapped.bam ${alignX_sample_out}/unmapped.tmp

		bam_merge -Q --sam-header ${alignX_sample_out}/tmp/header.sam ${alignX_sample_out}/unmapped.bam ${alignX_sample_out}/unmapped.tmp ${alignX_sample_out}/unmapped_multi.bam 2> /dev/null
		
		rm -fr ${alignX_sample_out}/unmapped.tmp
			
		picard 10g FixMateInformation VALIDATION_STRINGENCY=LENIENT I=${alignX_sample_out}/unmapped.bam O=${alignX_sample_out}/unmapped_fixed.bam &> ${alignX_sample_out}/unmapped_fixed.log
		picard_input="INPUT=${alignX_sample_out}/unmapped_fixed.bam"
	fi
	
        picard 10g SamToFastq FASTQ=${alignX_sample_out}/${samplename}.unmapped_R1.fastq ${unmapped_second_end_fastq} ${picard_input} VALIDATION_STRINGENCY=LENIENT INCLUDE_NON_PF_READS=TRUE &> ${alignX_sample_out}/${samplename}.unmapped.SamToFastq.log.txt
	
	#rm -fr ${alignX_sample_out}/tmp/
}

function bbmap_alignment_step {
	samplename=$1
	step=$2
	alignX_sample_out=$3
	alignL_sample_out=$4
	fastqin=$5
	fastq2in=$6
	
	bbmap_in2=""
	unmapped_second_end_fastq=""
	# (4) read unmapped
	samtools_accepted=" -F 4 " 
	samtools_unmapped=" -f 4 "

	if [ -e ${fastq2in} ] 
	then
		bbmap_fastqin2="in2=${fastq2in}"
		c=0
		# os singletons encontrados devem conter os correspondentes termos "null"
		for i in `echo ${fastqin} | tr "," "\n"`; do 
			if [ ${c} -gt 0 ]; then
				bbmap_fastqin2="${bbmap_fastqin2},null"
			fi
			c=$((c+1))
		done
		
		unmapped_second_end_fastq="SECOND_END_FASTQ=${alignX_sample_out}/${samplename}.unmapped_R2.fastq"

		# (2) read mapped in proper pair 
		samtools_accepted=" -f 2 "
		samtools_unmapped=" -F 2 "
	fi

	cmdbbmap="bbwrap.sh -Xmx22g path=${refbbmapidx} build=1 in=${fastqin} ${bbmap_fastqin2} out=${alignX_sample_out}/all_hits.sam ambiguous=all trimq=0 mintrimlength=100 vslow=t strictmaxindel=t maxindel=3 minid=0.95 k=9 xstag=t mdtag=t nhtag=t xmtag=t amtag=t nmtag=t stoptag=t lengthtag=t idtag=t inserttag=t scoretag=t timetag=t boundstag=t xs=us sam=1.3 &> ${alignX_sample_out}/${samplename}.bbmap.log append=t"
	
	eval ${cmdbbmap}
	
	if [ ${step} == 3 ]
	then
		# (4) read unmapped
		samtools_accepted=" -F 4 "
		samtools_unmapped=" -f 4 "
	fi
		
	samtools view -b -S ${samtools_accepted} ${alignX_sample_out}/all_hits.sam > ${alignX_sample_out}/accepted_hits.bam 2> /dev/null
	samtools view -b -S ${samtools_unmapped} ${alignX_sample_out}/all_hits.sam > ${alignX_sample_out}/unmapped.bam 2> /dev/null

	rm -f ${alignX_sample_out}/all_hits.sam 
	
	picard_input="INPUT=${alignX_sample_out}/unmapped.bam"
	if [ -e ${fastq2in} ]
	then
		picard 10g FixMateInformation VALIDATION_STRINGENCY=LENIENT I=${alignX_sample_out}/unmapped.bam O=${alignX_sample_out}/unmapped_fixed.bam &> ${alignX_sample_out}/unmapped_fixed.log
		picard_input="INPUT=${alignX_sample_out}/unmapped_fixed.bam"
	fi

        picard 10g SamToFastq FASTQ=${alignX_sample_out}/${samplename}.unmapped_R1.fastq ${unmapped_second_end_fastq} ${picard_input} VALIDATION_STRINGENCY=LENIENT INCLUDE_NON_PF_READS=TRUE &> ${alignX_sample_out}/${samplename}.unmapped.SamToFastq.log.txt

}

function star_alignment_step {
	samplename=$1
	step=$2
	alignX_sample_out=$3
	alignL_sample_out=$4
	fastqin=$5
	fastq2in=$6
	
	star_fastq2in=""
	unmapped_second_end_fastq=""
	# (4) read unmapped
	samtools_accepted=" -F 4 " 
	samtools_unmapped=" -f 4 "
	
	if [ -e "$fastq2in" ]
	then
		star_fastq2in="$fastq2in"
		# (2) read mapped in proper pair 
		samtools_accepted=" -f 2 "
		samtools_unmapped=" -F 2 "
		unmapped_second_end_fastq="SECOND_END_FASTQ=${alignX_sample_out}/${samplename}.unmapped_R2.fastq"
	fi
	
	cmdstar="STAR --limitBAMsortRAM 20000000000 --runThreadN ${num_threads} --outFilterMultimapNmax 10 --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 0.1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique  --alignSJoverhangMin 35 --alignEndsType EndToEnd --alignIntronMin 5000000 --outSAMstrandField intronMotif --outSAMattributes All --outSAMprimaryFlag AllBestScore --bamRemoveDuplicatesType UniqueIdentical --outSAMtype BAM Unsorted --sjdbOverhang 100 --sjdbGTFfile ${refgff} --genomeDir ${refstaridx} --readFilesIn ${fastqin} ${star_fastq2in} --outFileNamePrefix ${alignX_sample_out}/ --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5 --alignMatesGapMax 600 --outSAMunmapped Within --chimOutType WithinBAM --chimSegmentMin 50"
	
	if [ -s "${alignL_sample_out}/all_junctions.junc" ]
	then   
		cmdstar="${cmdstar} --sjdbFileChrStartEnd ${alignL_sample_out}/all_junctions.junc "
	fi
	
	cmdstar="${cmdstar} &> ${alignX_sample_out}/${samplename}.STAR.log.txt"
	
	
	eval ${cmdstar}
	
	if [ -e "${alignX_sample_out}/SJ.out.tab" ]
	then
		perl -F"\t" -lane ' if (  ($F[4]>0) && ($F[5]>=0) && ($F[6]>=5) && ($F[7]>=0) && ($F[8]>10)  ) { print join("\t", @F[0,1,2],(($F[3]==1) ? "+" : (($F[3]==2) ? "-" : ".")) ); }' ${alignX_sample_out}/SJ.out.tab > ${alignX_sample_out}/junctions.junc
	else
		touch ${alignX_sample_out}/junctions.junc
	fi
	
	cat ${alignL_sample_out}/all_junctions.junc ${alignX_sample_out}/junctions.junc | nsort -u > ${alignX_sample_out}/all_junctions.junc
		
	
	if [ ${step} == 3 ]
	then
		# (4) read unmapped
		samtools_accepted=" -F 4 "
		samtools_unmapped=" -f 4 "
	fi
	
	samtools view -b ${samtools_accepted} ${alignX_sample_out}/Aligned.out.bam > ${alignX_sample_out}/accepted_hits.bam 2> /dev/null
	samtools view -h ${samtools_unmapped} ${alignX_sample_out}/Aligned.out.bam | fixPrimaryInformation.pl 2> ${alignX_sample_out}/fixPrimaryInformation.log.txt | samtools view -b - > ${alignX_sample_out}/unmapped.bam 2> /dev/null

	picard_input="INPUT=${alignX_sample_out}/unmapped.bam"
	if [ -e ${fastq2in} ]
	then
		picard 10g FixMateInformation VALIDATION_STRINGENCY=LENIENT I=${alignX_sample_out}/unmapped.bam O=${alignX_sample_out}/unmapped_fixed.bam &> ${alignX_sample_out}/unmapped_fixed.log
		picard_input="INPUT=${alignX_sample_out}/unmapped_fixed.bam"
	fi

        picard 10g SamToFastq FASTQ=${alignX_sample_out}/${samplename}.unmapped_R1.fastq ${unmapped_second_end_fastq} ${picard_input} VALIDATION_STRINGENCY=LENIENT INCLUDE_NON_PF_READS=TRUE &> ${alignX_sample_out}/${samplename}.unmapped.SamToFastq.log.txt

}

# removendo arquivo com o caminho de cada transcriptoma (.gtf)
rm -f ${cuffmerge_out}/assembly_GTF_list.txt

echo "* Pre-processing:"

for fastq in `ls ${fastx_out}/*_R1.fastq`; do
	
	# fqname - string com o nome do arquivo fastq
	
	fqname=`basename ${fastq} _R1.fastq`
	
	mate_fastq=`echo ${fastq} | sed 's/_R1.fastq$/_R2.fastq/'`

	echo "   ${fqname} ..."	
	
        # Execução do align2 para obter o alinhamento de cada réplica biológica em relação ao genoma
        # caso ainda não tenha sido executado
	
	align0_sample_out="${align0_out}/${fqname}"
	mkdir -p ${align0_sample_out}
	
        if [ !  -e "${align0_sample_out}/accepted_hits.bam" ]
        then
                echo "      Aligning ${fqname} against genome - step 0/3 (${aligner}) ..."
		
		if [ ${aligner} == "tophat2" ]; then
			
			tophat_alignment_step ${fqname} 0 ${align0_sample_out} ${align0_sample_out} ${fastq} ${mate_fastq}
			
		elif [ ${aligner} == "bbmap" ]; then
			
			bbmap_alignment_step ${fqname} 0 ${align0_sample_out} ${align0_sample_out} ${fastq} ${mate_fastq}

		elif [ ${aligner} == "star" ]; then

			star_alignment_step ${fqname} 0 ${align0_sample_out} ${align0_sample_out} ${fastq} ${mate_fastq}
		else
			echo "Not found aligner (${aligner})"
			exit
		fi

        fi
	
        if [ ! -e "${scythe_out}/${fqname}.scythe_R1.fastq" ]
        then   

                mkdir -p ${scythe_out}
                
		echo "         Adapter trimming - step 1/2 (scythe) [R1] ..."

                scythe -p 0.1 -m ${scythe_out}/${fqname}.scythe.matches_R1.txt -o ${scythe_out}/${fqname}.scythe_R1.fastq -n 5 -M 0  -q sanger  -a ${scythe_adapters_r1_path} -i ${align0_sample_out}/${fqname}.unmapped_R1.fastq  1> ${scythe_out}/${fqname}.scythe_R1.log.txt 2> ${scythe_out}/${fqname}.scythe_R1.err.txt

		cat ${scythe_out}/${fqname}.scythe_R1.fastq | sed 's/ *$//' > ${scythe_out}/${fqname}.scythe_R1.tmp
		
		mv ${scythe_out}/${fqname}.scythe_R1.tmp ${scythe_out}/${fqname}.scythe_R1.fastq
        fi
        
	if [ -e "${align0_sample_out}/${fqname}.unmapped_R2.fastq" ]
	then
		if [ ! -e "${scythe_out}/${fqname}.scythe_R2.fastq" ]
        	then   

			echo "         Adapter trimming - step 1/2 (scythe) [R2] ..."

        	        scythe -p 0.1 -m ${scythe_out}/${fqname}.scythe.matches_R2.txt -o ${scythe_out}/${fqname}.scythe_R2.fastq -n 5 -M 0  -q sanger  -a ${scythe_adapters_r2_path} -i ${align0_sample_out}/${fqname}.unmapped_R2.fastq  1> ${scythe_out}/${fqname}.scythe_R2.log.txt 2> ${scythe_out}/${fqname}.scythe_R2.err.txt

			cat ${scythe_out}/${fqname}.scythe_R2.fastq | sed 's/ *$//' > ${scythe_out}/${fqname}.scythe_R2.tmp
		
			mv ${scythe_out}/${fqname}.scythe_R2.tmp ${scythe_out}/${fqname}.scythe_R2.fastq
        	fi
	fi

        # Execução do align2 para obter o alinhamento de cada réplica biológica em relação ao genoma
        # caso ainda não tenha sido executado
	
	align1_sample_out="${align1_out}/${fqname}"
	mkdir -p ${align1_sample_out}

        if [ !  -e "${align1_sample_out}/accepted_hits.bam" ]
        then
                echo "      Aligning ${fqname} against genome - step 1/3 (${aligner}) ..."
		
		if [ ${aligner} == "tophat2" ]; then

			tophat_alignment_step ${fqname} 1 ${align1_sample_out} ${align0_sample_out} "${scythe_out}/${fqname}.scythe_R1.fastq" "${scythe_out}/${fqname}.scythe_R2.fastq"

		elif [ ${aligner} == "bbmap" ]; then
			
			bbmap_alignment_step ${fqname} 1 ${align1_sample_out} ${align0_sample_out} "${scythe_out}/${fqname}.scythe_R1.fastq" "${scythe_out}/${fqname}.scythe_R2.fastq"

		elif [ ${aligner} == "star" ]; then
	
			star_alignment_step ${fqname} 1 ${align1_sample_out} ${align0_sample_out} "${scythe_out}/${fqname}.scythe_R1.fastq" "${scythe_out}/${fqname}.scythe_R2.fastq"		
		else 
			echo "Not found aligner (${aligner})"
			exit
		fi
        fi

        if [ ! -e "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" ]
        then
                
		mkdir -p ${cutadapt_out}

                echo "         Adapter trimming - step 2/2 (cutadapt) [R1] ..."
		
		cutadapt_5p_adapter_r1=`perl -ane 'chomp; if ($_=~/^>(\S+)/) { } else { print "--front=\"$_\""," "; }' ${cutadapt_5p_adapters_r1_path}`;

                cutadapt --format=fastq ${cutadapt_5p_adapter_r1} --error-rate=0.1 --times=2 --overlap=5 --minimum-length=0 --output=${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq  ${align1_sample_out}/${fqname}.unmapped_R1.fastq > ${cutadapt_out}/${fqname}.scythe.cutadapt_R1.report.txt 2> ${cutadapt_out}/${fqname}.scythe.cutadapt_R1.log.txt

        fi
	if [ -e "${align1_sample_out}/${fqname}.unmapped_R2.fastq" ]
	then
		if [ ! -e "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq" ]
        	then
                	echo "         Adapter trimming - step 2/2 (cutadapt) [R2] ..."
		
			cutadapt_5p_adapter_r2=`perl -ane 'chomp; if ($_=~/^>(\S+)/) { } else { print "--front=\"$_\""," "; }' ${cutadapt_5p_adapters_r2_path}`;

	                cutadapt --format=fastq ${cutadapt_5p_adapter_r2} --error-rate=0.1 --times=2 --overlap=5 --minimum-length=0 --output=${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq  ${align1_sample_out}/${fqname}.unmapped_R2.fastq > ${cutadapt_out}/${fqname}.scythe.cutadapt_R2.report.txt 2> ${cutadapt_out}/${fqname}.scythe.cutadapt_R2.log.txt

		fi
	fi
	
	
        # Execução do align2 para obter o alinhamento de cada réplica biológica em relação ao genoma
        # caso ainda não tenha sido executado
	
	align2_sample_out="${align2_out}/${fqname}"
	mkdir -p ${align2_sample_out}
	
	if [ !  -e "${align2_sample_out}/accepted_hits.bam" ]
	then
		echo "      Aligning ${fqname} against genome - step 2/3 (${aligner}) ..."
		
		if [ ${aligner} == "tophat2" ]; then
			
			tophat_alignment_step ${fqname} 2 ${align2_sample_out} ${align1_sample_out} "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq"
			
		elif [ ${aligner} == "bbmap" ]; then
			
			bbmap_alignment_step ${fqname} 2 ${align2_sample_out} ${align1_sample_out} "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq"
			
		elif [ ${aligner} == "star" ]; then
			
			star_alignment_step ${fqname} 2 ${align2_sample_out} ${align1_sample_out} "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq"			
	
		else 
			echo "Not found aligner (${aligner})"
			exit
		fi
		
        fi

	if [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ] && [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
        then

                echo "         Quality and poly-A/T trimming (prinseq) ..."

                mkdir -p ${prinseq_out}
		

		fastq2in=""
		
		if [ -e "${align2_sample_out}/${fqname}.unmapped_R2.fastq" ]
		then
			fastq2in="-fastq2 ${align2_sample_out}/${fqname}.unmapped_R2.fastq"
		fi

                # trimming reads by quality value
                prinseq-lite.pl -fastq ${align2_sample_out}/${fqname}.unmapped_R1.fastq ${fastq2in} -out_format 3 -out_good ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq -out_bad ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq-bad -qual_noscale -no_qual_header -min_len 15 -ns_max_p 80 -noniupac -trim_qual_right 25 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 3 -trim_qual_step 1 &> ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.log.txt

        fi
	
        # Execução do align3 para obter o alinhamento de cada réplica biológica em relação ao genoma
        # caso ainda não tenha sido executado
	
	align3_sample_out="${align3_out}/${fqname}"
	mkdir -p ${align3_sample_out}

        if [ !  -e "${align3_sample_out}/accepted_hits.bam" ]
        then
                echo "      Aligning ${fqname} against genome - step 3/3 (${aligner}) ..."

		if [ ${aligner} == "tophat2" ]
		then
			
			tophat_fastq_in=""
			if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
			then
				tophat_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq"
			elif [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]
			then
				tophat_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq"
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ]; then
					tophat_fastq_in="${tophat_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq"
				fi
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ]; then
					tophat_fastq_in="${tophat_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
				fi
			else
				echo "Missing input file for tophat2 alignment step 3"
			fi

			tophat_alignment_step ${fqname} 3 ${align3_sample_out} ${align2_sample_out} ${tophat_fastq_in} "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq"

		
		elif [ ${aligner} == "bbmap" ]; then
			
			bbmap_fastq_in=""
			if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
			then
				bbmap_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq"
			elif [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]
			then
				bbmap_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq"
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ]; then
					bbmap_fastq_in="${bbmap_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq"
				fi
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ]; then
					bbmap_fastq_in="${bbmap_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
				fi
			else
				echo "Missing input file for bbmap alignment step 3"
			fi
				
			
			bbmap_alignment_step ${fqname} 3 ${align3_sample_out} ${align2_sample_out} ${bbmap_fastq_in} "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq"
			
		else
			star_fastq_in=""
			if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
			then
				star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq"
			elif [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]
			then
				star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq"
			else
				echo "Missing input file for star alignment step 3"
			fi
			
			star_alignment_step ${fqname} 3 ${align3_sample_out} ${align2_sample_out} ${star_fastq_in} "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq"
				
			if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ]; then
				star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq"
				mkdir -p ${align3_sample_out}/singleton1
				star_alignment_step ${fqname} 3 ${align3_sample_out}/singleton1 ${align3_sample_out} ${star_fastq_in} ""
				
			fi
			if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ]; then
				star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
				mkdir -p ${align3_sample_out}/singleton2
				star_alignment_step ${fqname} 3 ${align3_sample_out}/singleton2 ${align3_sample_out} ${star_fastq_in} ""
			fi
		
			if [ ! -e "${align3_sample_out}/header.sam" ]
			then 
				samtools view -H ${align3_sample_out}/accepted_hits.bam > ${align3_sample_out}/header.sam 2> /dev/null
			fi
			
			if [ -e "${align3_sample_out}/singleton1/accepted_hits.bam" ]; then
				mv ${align3_sample_out}/accepted_hits.bam ${align3_sample_out}/accepted_hits.tmp

		                bam_merge -Q --sam-header ${align3_sample_out}/header.sam ${align3_sample_out}/accepted_hits.bam ${align3_sample_out}/accepted_hits.tmp ${align3_sample_out}/singleton1/accepted_hits.bam 2> /dev/null
				
				rm -f ${align3_sample_out}/accepted_hits.tmp
			fi
			
			if [ -e "${align3_sample_out}/singleton2/accepted_hits.bam" ]; then
				mv ${align3_sample_out}/accepted_hits.bam ${align3_sample_out}/accepted_hits.tmp

		                bam_merge -Q --sam-header ${align3_sample_out}/header.sam ${align3_sample_out}/accepted_hits.bam ${align3_sample_out}/accepted_hits.tmp ${align3_sample_out}/singleton2/accepted_hits.bam 2> /dev/null
				
				rm -f ${align3_sample_out}/accepted_hits.tmp
			fi

		fi
		
        fi


        finalalign_sample_out="${finalalign_out}/${fqname}"
        mkdir -p ${finalalign_sample_out}

        if [ ! -e "${finalalign_sample_out}/accepted_hits.bam" ]
        then
                echo "      Merging alignments of ${fqname} against genome (bam_merge) ..."

		if [ ! -e "${finalalign_sample_out}/header.sam" ]
		then 
			samtools view -H ${align0_sample_out}/accepted_hits.bam > ${finalalign_sample_out}/header.sam 2> /dev/null
		fi

                bam_merge -Q --sam-header ${finalalign_sample_out}/header.sam ${finalalign_sample_out}/accepted_hits.bam ${align0_sample_out}/accepted_hits.bam ${align1_sample_out}/accepted_hits.bam ${align2_sample_out}/accepted_hits.bam ${align3_sample_out}/accepted_hits.bam 2> /dev/null

                samtools sort ${finalalign_sample_out}/accepted_hits.bam ${finalalign_sample_out}/accepted_hits_sorted

		mv ${finalalign_sample_out}/accepted_hits_sorted.bam ${finalalign_sample_out}/accepted_hits.bam

        fi
	
	
	# Execução do cufflinks para obter a montagem para cada réplica biológica
	# caso ainda não tenha sido executado
	cufflinks_sample_out="${cufflinks_out}/${fqname}"
	if [ !  -e "${cufflinks_sample_out}/transcripts.gtf" ]
	then

		mkdir -p ${cufflinks_sample_out}

		echo "      Assembling ${fqname} (cufflinks) ..."

		cufflinks --num-threads ${num_threads} --library-type fr-unstranded --multi-read-correct --GTF ${refgff} --frag-bias-correct ${refgenomefa} --total-hits-norm --output-dir ${cufflinks_sample_out} ${finalalign_sample_out}/accepted_hits.bam --max-bundle-frags 2000000 &> ${cufflinks_sample_out}/${fqname}.cufflinks.log.txt

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
		cuffquant --library-type fr-unstranded --multi-read-correct --frag-bias-correct ${refgenomefa} ${cuffmerge_out}/merged.gtf ${finalalign_out}/${bs}/accepted_hits.bam -o ${cuffquant_out}/${bs} --max-bundle-frags 999999999 &> ${cuffquant_out}/${bs}.log.txt
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

if [ ! -e "${cuffnorm_out}/run.info" ]
then
	echo "* Scaling transcriptomes (cuffnorm) ..."
	
	cuffnorm --num-threads ${num_threads} --library-type fr-unstranded --library-norm-method geometric --output-format simple-table --compatible-hits-norm --output-dir ${cuffnorm_out} --labels $(IFS=, ; echo "${biogroups[*]}") ${cuffmerge_out}/merged.gtf ${biogroupscxb[@]} &> ${cuffnorm_out}/cuffnorm.log
fi 

# Execução do cuffdiff para obter as estimativas de expressão gênica diferencial
# caso não tenha sido executado
if [ !  -e "${cuffdiff_out}/run.info" ]
then
	echo "* Estimating Differential Expression (cuffdiff) among groups:" "${biogroups[*]}" 
	cuffdiff --num-threads ${num_threads} --library-type fr-unstranded --multi-read-correct --frag-bias-correct ${refgenomefa} --library-norm-method geometric --compatible-hits-norm --output-dir ${cuffdiff_out} ${cuffmerge_out}/merged.gtf ${biogroupscxb[@]} --FDR  0.05 --labels $(IFS=, ; echo "${biogroups[*]}") --max-bundle-frags 999999999 &> ${cuffdiff_out}/cuffdiff.log
fi

# Documentação do bash
#	http://www.gnu.org/software/bash/manual/bashref.html#Arrays
# Discussão sobre * ou @ ao expandir valores do array
#	Suppose we have perls=(perl-one perl-two) and use ls "${perls[*]}" -- that's equivalent to ls "perl-one perl-two", which will look for single file named perl-one perl-two, which is probably not what you wanted. ls "${perls[@]}" is equivalent to ls "perl-one" "perl-two", which is much more likely to do something useful.
# 	http://www.coderexception.com/0166bbz36XXQJQQX/a-confusion-about-array-versus-array-in-the-context-of-a-bash-completion


