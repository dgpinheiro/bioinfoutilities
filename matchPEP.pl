#!/usr/bin/perl
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

use strict;
use warnings;

my $indir=$ARGV[0];

die "Missing input directory" unless ($indir);

my $blastdir=$ARGV[1];

die "Missing blast directory" unless ($blastdir);

my $txtdir=$ARGV[2];

die "Missing txt directory" unless ($txtdir);

use File::Basename;
use Math::Combinatorics;

my $num_threads=14;

my @db;

print STDERR "Verifying fasta files on $indir ...\n";
foreach my $i (glob("$indir/*.fa")) {
	my $bn=basename(${i},'.fa');

	unless (-e "$indir/${bn}.nin") { 
	    print STDERR "Formatting $bn ...\n";
        `makeblastdb -in $i -dbtype prot -out $indir/$bn -title $bn`;
	}
	else {
		print STDERR "DB $bn already formatted !\n";
	}
	push(@db, $bn);
}

print STDERR "Generating reference blast scores ...\n";
foreach my $d (@db) {
    print STDERR "\t$d\n";
    my $blastcmd="blastp -db $indir/$d -query $indir/$d.fa -num_threads $num_threads -word_size 2 -seg no -matrix PAM30 -comp_based_stats F -out $blastdir/$d".'_x_'.$d.'.bls';
    #print $blastcmd,"\n";
    `$blastcmd`;
	my $extractREFcmd="./extractREF.pl -i $blastdir/$d".'_x_'.$d.".bls > $txtdir/$d".'_x_'.$d."_out.txt -l WARN 2> $txtdir/$d".'_x_'.$d.'_err.txt';
    #print $extractREFcmd,"\n";
    `$extractREFcmd`;
}

print STDERR "Making comparisons ...\n";
foreach my $ar_c (combine(2,@db)) {

    my $one = $ar_c->[0];
    my $two = $ar_c->[1];

	print STDERR "\t$two x $one ...\n";
	
	unless (-e $blastdir.'/'.$two.'_x_'.$one.'.bls' ) {
		print STDERR "\tBLASTP $two x $one\n";
		my $blastcmd="blastp -db $indir/$one -query $indir/$two.fa -num_threads $num_threads -word_size 2 -seg no -matrix PAM30 -comp_based_stats F -out $blastdir/$two".'_x_'.$one.'.bls';
        #print $blastcmd,"\n";
        `$blastcmd`;
    }
	
	unless (-e $txtdir.'/'.$two.'_x_'.$one.'.txt' ) {
            my $exactPEPcmd = "./extractEXACTpep.pl -i1 $indir/$two.fa -i2 $indir/$one.fa -b ../bed -l WARN > $txtdir/EXACT_$two".'_x_'.$one.'_out.txt 2> '.$txtdir."/EXACT_$two".'_x_'.$one.'_err.txt';
            #print $exactPEPcmd,"\n";
            `$exactPEPcmd`;

			my $extractPEPcmd="./extractBH.pl -i $blastdir/$two".'_x_'.$one.".bls -rq $txtdir/$two".'_x_'.$two."_out.txt -rs $txtdir/$one".'_x_'.$one."_out.txt -t 0.90 -e $txtdir/EXACT_$two".'_x_'.$one."_out.txt -b ../bed -l WARN > $txtdir/$two".'_x_'.$one."_out.txt 2> $txtdir/$two".'_x_'.$one.'_err.txt';
            #print $extractPEPcmd,"\n";
            `$extractPEPcmd`;
    }
	
    print STDERR "\t$one x $two ...\n";

	unless (-e $blastdir.'/'.$one.'_x_'.$two.'.bls' ) {
        print STDERR "\tBLASTP $one x $two\n";
		my $blastcmd="blastp -db $indir/$two -query $indir/$one.fa -num_threads $num_threads -word_size 2 -seg no -matrix PAM30 -comp_based_stats F -out $blastdir/$one".'_x_'.$two.'.bls';
        #print $blastcmd,"\n";
        `$blastcmd`;
	}

	unless (-e $txtdir.'/'.$one.'_x_'.$two.'.txt' ) {
			my $exactPEPcmd="./extractEXACTpep.pl -i1 $indir/$one.fa -i2 $indir/$two.fa  -b ../bed -l WARN > $txtdir/EXACT_$one".'_x_'.$two."_out.txt 2> $txtdir/EXACT_".$one.'_x_'.$two.'_err.txt';
            #print $exactPEPcmd,"\n";
            `$exactPEPcmd`;
            my $extractPEPcmd="./extractBH.pl -i $blastdir/$one".'_x_'.$two.".bls -rq $txtdir/$one".'_x_'.$one."_out.txt -rs $txtdir/$two".'_x_'.$two."_out.txt -t 0.90 -e $txtdir/EXACT_$one".'_x_'.$two."_out.txt -b ../bed -l WARN > $txtdir/$one".'_x_'.$two."_out.txt 2> $txtdir/${one}".'_x_'.$two.'_err.txt';
            #print $extractPEPcmd,"\n";
            `$extractPEPcmd`;
    }
}

