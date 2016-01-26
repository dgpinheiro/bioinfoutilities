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


# Diretório onde está o diretório ./processed/prinseq com os arquivos .fastq processados
input="$1"

# Nome para os contigs/scaffolds montados
asmname="$2"

# Se plataforma Ion Torrent (exigido para spades.py) inclui parâmetro --iontorrent na linha de comando
# habilita correção de erros próprios do Ion Torrent
iontorrent=""
if [ ! ${iontorrent} ]
then   
	iontorrent="--iontorrent"
fi

# Quantidade de processadores disponíveis na máquina
proc_number=$(grep -c '^processor' /proc/cpuinfo)
# Utilizar 90% dos processadores 
proc_percent=90

# Número de threads disponíveis para os processos

# Seria mais fácil determinar que o número vai de 0-1 para proc_percent e usar a linha a seguir
# bc sem o l ao dividir por 1 ele trunca o valor decimal, deixando-o somente a parte inteira
#threads=$( echo "(${proc_number} * 0.9)/1" | bc )

# Mas só para deixar um pouco mais complicado, além disso, arredondando o valor e não truncando-o.
threads=$( echo "(${proc_number} * (${proc_percent}/100))" | bc -l | xargs printf "%.0f" )

# Número de threads para o processamento dos arquivos pmap utilizando parallel
# Obter o número de 10% dos processadores para lançar os processos com o parallel para o processamento dos arquivo *.pmap,
# os quais contêm as reads particionadas
pmap_proc_percent=10
pmap_threads=$( echo "(${proc_number} * (${pmap_proc_percent}/100))" | bc -l | xargs printf "%.0f" )

# http://khmer.readthedocs.org/en/v1.0/choosing-table-sizes.html
# Memória = -N (n_tables) * -x (tablesize)
# 20Gb
# khmer parâmetro -N
khmer_N="4"
# khmer parâmetro -x
# khmer_byte_x = 5Gb
# 4 * khmer_byte_x = 20Gb (total disponível de memória)
# e9 = Gb
khmer_byte_x="5e9"
# avail_mem_gb = 4 * khmer_byte_x
# comando interno pega o inteiro antes do e o múltiplica por khmer_N
avail_mem="$(($(echo ${khmer_byte_x} | cut -d e -f 1) * ${khmer_N}))"

# khmer parâmetro -k
khmer_k=32

# Tamanho mínimo para a avaliação das montagens, realizada pela soma dos tamanhos de todos os contigs maiors que este valor
cutoff=300

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

echo "Assembling data from ${input}..."

echo "   Digital normalization - \"diginorm\" ..."

curdir=`pwd`

cd ${input}

mkdir -p diginorm

cd ./diginorm

echo "      First normalization C 20 ..."
### Normaliza tudo para uma cobertura de 20, considerando um k-mer de 20
normalize-by-median.py -k ${khmer_k} -C 20 -N ${khmer_N} -x ${khmer_byte_x} --savetable normC20k20.kh ../processed/prinseq/*.fastq > normC20k20.out.txt 2> normC20k20.err.txt

echo "      Filtering low abundance k-mers ..."
### Poda leituras em k-mers que são pouco abundantes em leituras de alta cobertura. A opção -V é usada para datasets com cobertura variável.
filter-abund.py --threads ${threads} -V normC20k20.kh *.keep > filter-abund.out.txt 2> filter-abund.err.txt

echo "      Second normalization C 5 ..."
### Após eliminar k-mers errôneos, vamos abandonar mais alguns dados de alta cobertura.
normalize-by-median.py -C 5 -k ${khmer_k} -N ${khmer_N} -x ${khmer_byte_x} --savetable normC5k20.kh *.keep.abundfilt > normC5k20.out.txt 2> normC5k20.err.txt

echo "      Compacting resulting files ..."
# Compactar cada amostra (gzip)
for i in `ls *.scythe.cutadapt5p.filtered.prinseq.fastq.keep.abundfilt.keep`; do bn=`basename ${i} .scythe.cutadapt5p.filtered.prinseq.fastq.keep.abundfilt.keep`; gzip -9c ${i} > ${bn}.scfp.kak.fastq.gz; done

# Remover arquivos desnecessários
rm -f normC20k20.kh *.keep *.abundfilt

echo "      Generating reads statistics ..."
readstats.py *.kak.fastq.gz ../processed/prinseq/*.fastq > diginorm.out.txt 2> diginorm.err.txt

cd ../

echo "   Partitioning ..."

mkdir -p ./partitioned

cd ./partitioned

### Eliminação dos k-mers altamente repetitivos que podem juntar múltiplas espécies (quimeras) e renomear de forma apropriada:

filter-below-abund.py ../diginorm/normC5k20.kh ../diginorm/*.fastq.gz > filter-below-abund.out.txt 2> filter-below-abund.err.txt

## Carrega grafo de k-mers
load-graph.py -k ${khmer_k} -T ${threads} -N ${khmer_N} -x ${khmer_byte_x} lump *.kak.fastq.gz.below > load-graph.out.txt 2> load-graph.err.txt

## Encontrando k-mers altamente conectados iniciais (possíveis artefatos)
make-initial-stoptags.py -k ${khmer_k} -N ${khmer_N} -x ${khmer_byte_x} lump > make-initial-stoptags.out.txt 2> make-initial-stoptags.err.txt

## Particiona o grafo de acordo com a sua conectividade
partition-graph.py --threads ${threads}  --stoptags lump.stoptags lump > partition-graph.out.txt 2> partition-graph.err.txt

## Encontrando k-mers altamente conectados (possíveis artefatos)
pmap_count=`ls -l lump.*.pmap | wc -l`

pmap_limit=$((pmap_count / pmap_threads))

echo "      PMAP count: ${pmap_count}"

pmap_c=0
pmap_dir=0
rm -f ./run-find-knots.sh
for i in lump.*.pmap; do
        pmap_c=$((pmap_c + 1))
        if [ $((pmap_c % pmap_limit)) == 1 ]; then
                pmap_dir=$((pmap_dir + 1))
		mkdir -p ./fk.${pmap_dir}
		ln -s $(readlink -f lump.pt) ./fk.${pmap_dir}/lump.pt
		ln -s $(readlink -f lump.tagset) ./fk.${pmap_dir}/lump.tagset
		echo "(cd ./fk.${pmap_dir} && find-knots.py lump > find-knots.out.txt 2> find-knots.err.txt)" >> ./run-find-knots.sh
        fi
        echo "         Move ${i} to fk.${pmap_dir}"
	mv ${i} ./fk.${pmap_dir}/
done

echo "      Running find-knots.sh using ${pmap_threads} parallel processes ..."

parallel --gnu -j ${pmap_threads} < ./run-find-knots.sh

echo "      Merging all .stoptags files to merge.stoptags ..."

merge-stoptags.py -k ${khmer_k} ./fk > merge-stoptags.out.txt 2> merge-stoptags.err.txt

echo "      Filtering stoptags ..."

### Poda sequências nos k-mers considerados artefatos
filter-stoptags.py -k ${khmer_k} merge.stoptags *.kak.fastq.gz.below > filter-stoptags.out.txt 2> filter-stoptags.err.txt

echo "      do-partitions ..."

### Particionamento, gera arquivos que contêm as anotações das partições
do-partition.py -N ${khmer_N} -k ${khmer_k} -x ${khmer_byte_x} -T ${threads} kak *.kak.fastq.gz.below.stopfilt > do-partition.out.txt 2> do-partition.err.txt

echo "      extract-partitions ..."

### Extraindo partições em grupos
extract-partitions.py -m 0 -X 1000000 kak *.part > extract-partitions.out.txt 2> extract-partitions.err.txt

echo "      sweep-reads ..."

### Alocar os dados processados (sem normalização) em partições - Recupera reads baseadas no compartilhamento de k-mers
cat ../processed/prinseq/*.fastq > all.fq

sweep-reads3.py -k ${khmer_k} -N ${khmer_N} -x ${khmer_byte_x} kak.group*.fq all.fq > sweep-reads3.out.txt 2> sweep-reads3.err.txt


# Renomeando

for i in kak*.fq; do name=$(basename $i .fq); mv ${i} ${name}.dn.fq; done
for i in kak*.sweep3; do name=$(basename $i .fq.sweep3); mv ${i} ${name}.nodn.fq; done

# Remover arquivos desnecessários
rm all.fq
rm -f *.part *.sweep3 *.below.fq
rm -f ../diginorm/*.kh

cd ../

echo "   Assembling ..."

mkdir -p ./assembled

cd ./assembled

groups=()
for i in ../partitioned/*.fq; do
   name=$(basename $i .fq);

   echo "   [${name}] ..."
   
   echo "      IDBA-UD ..."
   # scaffold.fa
   
   # O arquivo de entrada para o idba_ud deve ser um arquivo do tipo paired-end, neste caso o script fastq-se2ipe.pl vai
   # gerar uma sequência complementar reversa para cada sequência e posicioná-las logo após cada sequência de origem, de
   # modo intercalado, e também deve estar no formato fasta, aqui utilizando fastq2fasta.pl
   mkdir -p ${name}.idba_ud.0.d
   fastq-se2ipe.pl ${i} | /usr/local/bioinfo/bioinfoutilities/fastq2fasta.pl -n 100 > ${name}.idba_ud.0.d/${name}.fa

   # --min_pairs zero pois são sequências de origem single-end transformadas em paired-end, neste caso não vamos exigir restrições com relação
   # a pares para determinar contigs/scaffolds
   idba_ud --num_threads ${threads} -r ${name}.idba_ud.0.d/${name}.fa --min_contig 100 --maxk 51 --step 8 --mink 27 --min_pairs 0 -o ${name}.idba_ud.0.d > ${name}.idba_ud.0.log.out.txt 2> ${name}.idba_ud.0.log.err.txt
   
   echo "      Newbler ..."
   # 454AllContigs.fna

   # -mi - mínimo de idêntidade em sobreposições
   # -ml - tamanho mínimo para as sobreposições
   # -force - sobrepor diretório já gravado (caso exista)
   runAssembly -o ${name}.newbler.0.d -mi 95 -ml 20 -force -cpu ${threads} ${i} > ${name}.newbler.0.log.out.txt 2> ${name}.newbler.0.log.err.txt

   echo "      MEGAHIT ..."
   # final.contigs.fa

   # desabilitado pois não há memória suficiente, mínimo ~30Gb
   # seria interessante utilizar --presets (pré-configurações) ex.: meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'       (large & complex metagenomes, like soil)
   # --min-count - múltiplicidade mínima para eliminar nós do grafo de k-mers
   megahit --mem-flag 2  -t $((threads/2)) --min-contig-len 100 -o ${name}.megahit.0.d -m 0.8 --min-count 2 --k-list 19,35,51 -r ${i} > ${name}.megahit.0.log.out.txt 2> ${name}.megahit.0.log.err.txt

   echo "      SPAdes ..."
   # scaffolds.fasta
   # desabilitado pois há problemas com alguns conjuntos de dados
   spades.py -t ${threads} --cov-cutoff auto -m ${avail_mem} ${iontorrent} -s ${i} -o ${name}.spades.0.d > ${name}.spades.0.log.out.txt 2> ${name}.spades.0.log.err.txt

   echo "      Minia ..."
   # é necessário criar o diretório antes de rodar o minia (https://www.biostars.org/p/168676/)
   # kak.contigs.fa
   for k in {19..51..16}; do
	echo "         ${k} ..."
	mkdir -p ./${name}.minia.$k.d
   	minia -nb-cores ${threads} -kmer-size $k -max-memory ${avail_mem} -out-dir ${name}.minia.$k.d -out ${name}.minia.$k.d/kak -contig-max-len 100 -in ${i} > ${name}.minia.$k.log.out.txt 2> ${name}.minia.$k.log.err.txt
   done

   echo "      ABySS ..."
#   # kak-scaffolds.fa
   for k in {19..51..16}; do
	echo "         ${k} ..."
	mkdir -p ./${name}.abyss.$k.d
	rp=`readlink -f ${i}`
	cd ./${name}.abyss.$k.d

	# http://seqanswers.com/wiki/ABySS
	# j = n. processos
	# se = arquivo .fastq (single-end)
	# name = nome base
	# e = cobertura de k-mers mínimas para determinar a extensão dos contigs, cobertura abaixo deste limite o contig indica o fim do contig
	# c = cobertura de k-mers mínima para um unitig (montagem em que não há ambiguidade, pode compor vários contigs)
	# s = tamanho mínimo dos unitigs para construir contigs
	# n = número mínimo de pares requeridos para construir contigs
	# S = tamanho mínimo dos contigs para construir scaffolds (um determinado contig pode fazer parte de apenas um scaffold)
	abyss-pe k=${k} name=kak j=${threads} se="${rp}" e=2 c=2 s=100 N=1 n=1 S=100 scaffolds > ../${name}.abyss.$k.log.out.txt 2> ../${name}.abyss.$k.log.err.txt
	cd ../
   done

   echo "      Velvet ..."
   # contigs.fa
   for k in {19..51..16}; do
	echo "         ${k} ..."
	velveth ${name}.velvet.$k.d $k -fastq -long $i > ${name}.velveth.$k.log.out.txt 2> ${name}.velveth.$k.log.err.txt
	# exp_cov - cobertura esperada (automático)
	# cov_cutoff - limiar para cobertura (automático)
	# conserveLong - conservar sequências longas
	velvetg ${name}.velvet.$k.d -exp_cov auto -cov_cutoff auto -scaffolding yes -conserveLong yes -min_contig_lgth 100 > ${name}.velvetg.$k.log.out.txt 2> ${name}.velvetg.$k.log.err.txt
   done

   echo "      Ray ..."
   # Scaffolds.fasta
   for k in {19..51..32}; do
	echo "         ${k} ..."
	#  -merge-seeds - Fusiona sementes de alinhamento (seeds) para reduzir o tempo de execução
	mpiexec -n ${threads} -num-cores ${threads} /usr/local/bioinfo/Ray-2.3.1/ray-build/Ray -use-minimum-seed-coverage 1 -minimum-seed-length 100 -merge-seeds -k ${k} -s ${i} -o ${name}.ray.$k.d > ${name}.ray.$k.log.out.txt 2> ${name}.ray.$k.log.err.txt
   done

   groups=($(printf "%s\n" ${groups[@]} `basename ${i}  | cut -d '.' -f 1,2` | sort -u ))

done

# Gera estatísticas básicas sobre a montagem
# statTrimmedN, statSum, statMax, filename
# 
assemstats3.py ${cutoff} *.idba_ud.*.d/scaffold.fa *.newbler.*.d/454AllContigs.fna *.megahit.*.d/final.contigs.fa *.spades.*.d/scaffolds.fasta *.minia.*.d/kak.contigs.fa *.velvet.*.d/contigs.fa *.ray.*.d/Scaffolds.fasta *.abyss.*/kak-scaffolds.fa > assemstats3.out.txt 2> assemstats3.err.txt

for g in ${groups[@]}; do
	echo "   Evaluating assembly of group ${g} ..."
	calc-best-assembly.py -C ${cutoff} -q ${g}.{*.idba_ud.*.d/scaffold.fa,*.newbler.*.d/454AllContigs.fna,*.megahit.*.d/final.contigs.fa,*.spades.*.d/scaffolds.fasta,*.minia.*.d/kak.contigs.fa,*.velvet.*.d/contigs.fa,*.ray.*.d/Scaffolds.fasta,*.abyss.*/kak-scaffolds.fa} -o ${g}.best.fa
done > best-assemblies.txt

multi-rename.py ${asmname} *.best.fa > final-assembly.fa

assemblathon_stats.pl final-assembly.fa > assemblathon_stats.out.txt 2> assemblathon_stats.err.txt

cd ${curdir}
