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
#  Copyright (C) 2019  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
#  Faculdade de Ciências Agrárias e Veterinárias (FCAV)
#  Laboratório de Bioinformática (LB)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://www.fcav.unesp.br 
#

infile=$1

if [ ! ${infile} ]; then
	echo "[ERROR] Missing input fasta file" 1>&2
	exit
else
	if [ ! -e ${infile} ]; then
		echo "[ERROR] Wrong input fasta file (${infile})" 1>&2
		exit
	fi
fi

fwdprimer=$2

if [ ! ${fwdprimer} ]; then
	echo "[ERROR] Missing forward primer sequence" 1>&2
	exit
else
	if [[ ! ${fwdprimer} =~ ^[ACGTURYSWKMBDHVN]+$ ]]; then
		echo "[ERROR] Wrong forward primer sequence (${fwdprimer})" 1>&2
		exit
	fi
fi

revprimer=$3

if [ ! ${revprimer} ]; then
	echo "[ERROR] Missing reverse primer sequence" 1>&2
	exit
else
	if [[ ! ${revprimer} =~ ^[ACGTURYSWKMBDHVN]+$ ]]; then
		echo "[ERROR] Wrong reverse primer sequence (${revprimer})" 1>&2
		exit
	fi
fi

taboutfile=$4

if [ ! ${taboutfile} ]; then
	echo "[ERROR] Missing tabbed text output file" 1>&2
	exit
else
	if [ ! -e ${taboutfile} ]; then
		echo "[ERROR] Wrong tabbed text output file (${taboutfile})" 1>&2
		exit
	fi
fi

fastaoutfile=$5

cmdline=""

if [ ! ${fastaoutfile} ]; then
	echo "[WARNING] Missing fasta output file" 1>&2
else
	if [ ! -e ${fastaoutfile} ]; then
		echo "[ERROR] Wrong fasta output file (${fastaoutfile})" 1>&2
		exit
	else
		cmdline="-fastaout ${fastaoutfile}"
	fi
fi

totlen=$((${#fwdprimer}+${#revprimer}))

usearch11 -search_pcr2 ${infile} -fwdprimer ${fwdprimer} \
				 -revprimer ${revprimer} \
				 -strand both \
				 -maxdiffs 3 \
				 -minamp ${totlen} \
				 -tabbedout ${taboutfile} ${cmdline}

