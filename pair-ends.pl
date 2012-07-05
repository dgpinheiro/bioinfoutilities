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
#  Copyright (C) 2010  Fundação Hemocentro de Ribeirão Preto
#
#  Laboratório de Genética Molecular e Bioinformática
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
# $Id$

=head1 NAME

    pair-ends.pl - Pair sequence "ends" from paired-end experiments after filtering steps.

=head1 SYNOPSIS
    
    # original files before filtering.: PEExp_1.fastq / PEExp_2.fastq
    # files after filtering...........: PEExp_1_cleaned.fastq / PEExp_2_cleaned.fastq
    # output paired-end files.........: PEExp_1_cleaned_paired.fastq / PEExp_2_cleaned_paired.fastq
    # output single-end files.........: PEExp_cleaned_single.fastq

    $ pair-ends.pl    -g1 PEExp_1.fastq -g2 PEExp_2.fastq \
    > -i1 PEExp_1_cleaned.fastq -i2 PEExp_2_cleaned.fastq \
    > -o1 PEExp_1_cleaned_paired.fastq -o2 PEExp_2_cleaned_paired.fastq -os PEExp_cleaned_single.fastq ;

=head1 ABSTRACT
    
=head1 DESCRIPTION
    
    Script to pair the sequences of two ordered fastq files from paired-end 
    experiments after filtering steps. The pairing is based on their original
    files before filtering procedures.

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2010 Regional Blood Center of Ribeirão Preto

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Getopt::Long;

my ($infile1, $infile2, $outfile1, $outfile2, $outfiles, $guidefile1, $guidefile2);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"  => sub { &Usage(); },
            "g1|guidefile1=s"=> \$guidefile1,
            "g2|guidefile2=s"=> \$guidefile2,
            "i1|inputfile1=s"=> \$infile1,
            "i2|inputfile2=s"=> \$infile2,
            "o1|outputfile1=s"=> \$outfile1,
            "o2|outputfile2=s"=> \$outfile2,
            "os|outputfiles=s"=> \$outfiles
) or &Usage();

die "Missing guide file 1 " unless ($guidefile1);
die "Wrong guide file 1" unless (-e $guidefile1);

die "Missing guide file 2" unless ($guidefile2);
die "Wrong guide file 2" unless (-e $guidefile2);

die "Missing input file 1" unless ($infile1);
die "Wrong input file 1" unless (-e $infile1);

die "Missing input file 2" unless ($infile2);
die "Wrong input file 2" unless (-e $infile2);

die "Missing output file 1" unless ($outfile1);
die "Missing output file 2" unless ($outfile2);
die "Missing output file s" unless ($outfiles);

open(ONEGUIDE, "<", $guidefile1) or die "Cannot open guide file 1: $!";
open(TWOGUIDE, "<", $guidefile2) or die "Cannot open guide file 2: $!";

open(ONEIN, "<", $infile1) or die "Cannot open input file 1: $!";
open(TWOIN, "<", $infile2) or die "Cannot open input file 2: $!";

open(ONEOUT, ">", $outfile1) or die "Cannot open output file 1: $!";
open(TWOOUT, ">", $outfile2) or die "Cannot open output file 2: $!";
open(SINGLEOUT, ">", $outfiles) or die "Cannot open output file s: $!";

my $discard_cnt = 0;
my $count = 0;
my $single_cnt = 0;
my $paired_cnt = 0;

my $inp1_seqid = <ONEIN>;
    
my $inp2_seqid = <TWOIN>;

$| = 1;
while(my $guidep1_seqid = <ONEGUIDE>) {
    my $guidep2_seqid = <TWOGUIDE>;

    # not used variables
    {
        my $guidep1_seq = <ONEGUIDE>;
        my $guidep1_qualid = <ONEGUIDE>;
        my $guidep1_qual = <ONEGUIDE>;
        
        my $guidep2_seq = <TWOGUIDE>;
        my $guidep2_qualid = <TWOGUIDE>;
        my $guidep2_qual = <TWOGUIDE>;
    }
     
    if (($inp1_seqid)&&($inp1_seqid eq $guidep1_seqid)) {
        my $inp1_seq = <ONEIN>;
        my $inp1_qualid = <ONEIN>;
        my $inp1_qual = <ONEIN>;
        
        if ($inp2_seqid eq $guidep2_seqid) {
            my $inp2_seq = <TWOIN>;
            my $inp2_qualid = <TWOIN>;
            my $inp2_qual = <TWOIN>;
            
            $paired_cnt++;
            
            print ONEOUT $inp1_seqid,$inp1_seq,$inp1_qualid,$inp1_qual;
            print TWOOUT $inp2_seqid,$inp2_seq,$inp2_qualid,$inp2_qual;
            
            $inp1_seqid = <ONEIN>;
            $inp2_seqid = <TWOIN>;

        }
        else {
            $single_cnt++;
            
            print SINGLEOUT $inp1_seqid,$inp1_seq,$inp1_qualid,$inp1_qual;
            $inp1_seqid = <ONEIN>;

        }
    }
    else {

        if (($inp2_seqid)&&($inp2_seqid eq $guidep2_seqid)) {
            my $inp2_seq = <TWOIN>;
            my $inp2_qualid = <TWOIN>;
            my $inp2_qual = <TWOIN>;
            
            $single_cnt++;
            
            print SINGLEOUT $inp2_seqid,$inp2_seq,$inp2_qualid,$inp2_qual;
            $inp2_seqid = <TWOIN>;

        }
        else {
            
            $discard_cnt++;
        }
    }

    $count++;
    if ($count % 10000 == 0) {
        print STDERR "Records: $count, Paired: $paired_cnt, Single: $single_cnt, Discarded: $discard_cnt                                                         \r";
    }

}

print STDERR "Records: $count, Paired: $paired_cnt, Single: $single_cnt, Discarded: $discard_cnt                                                         \n";

close(ONEGUIDE);
close(TWOGUIDE);

close(ONEIN);
close(TWOIN);

close(ONEOUT);
close(TWOOUT);
close(SINGLEOUT);

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2010 Regional Blood Center of Ribeirão Preto

Usage

        $0	[-h/--help] [-g1 PEExp_1.fastq -g2 PEExp_2.fastq 
                         -i1 PEExp_1_cleaned.fastq -i2 PEExp_2_cleaned.fastq 
                         -o1 PEExp_1_cleaned_paired.fastq -o2 PEExp_2_cleaned_paired.fastq -os PEExp_cleaned_single.fastq ]

Argument(s)

        -h      --help          Help
        
        -g1     --guidefile1    Guide file 1 (Original fastq p1 file - pre-filtering)
        -g2     --guidefile2    Guide file 2 (Original fastq p2 file - pre-filtering)

        -i1     --inputfile1    Input file 1 (Filtered fastq p1 file - post-filtering)
        -i2     --inputfile2    Input file 2 (Filtered fastq p2 file - post-filtering)

        -o1     --outputfile1   Output file 1
        -o2     --outputfile2   Output file 2
        -os     --outputfiles   Output file s

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

