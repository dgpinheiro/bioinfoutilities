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

    merge-ends-fastq.pl - Merge sequence "ends" from paired-end experiments creating another fastq file.

=head1 SYNOPSIS

    # Merge two fastq files from a paired-end experiment: PEExp_1.fastq and PEExp_2.fastq and creates
    # another fastq file. The sequence ids were merged using the specified separator "|".

    $ merge-ends.pl   -i1 PEExp_1.fastq -i2 PEExp_2.fastq \
    > -o PEExp_12.fastq -s "|"

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Merge two files from a paired-end experiment into another fastq file using
    an specific separator (<sep>) for ids from P1 and P2 sequences.

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

my ($infile1, $infile2, $outfile, $sep);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"  => sub { &Usage(); },
            "i1|inputfile1=s"=> \$infile1,
            "i2|inputfile2=s"=> \$infile2,
            "o|outputfile=s"=> \$outfile,
            "s|separator=s"=> \$sep
) or &Usage();

$sep||="\t";

die "Missing input file 1" unless ($infile1);
die "Wrong input file 1" unless (-e $infile1);

die "Missing input file 2" unless ($infile2);
die "Wrong input file 2" unless (-e $infile2);

die "Missing output file" unless ($outfile);

open(ONEIN, "<", $infile1) or die "Cannot open input file 1: $!";
open(TWOIN, "<", $infile2) or die "Cannot open input file 2: $!";

open(OUT, ">", $outfile) or die "Cannot open output file: $!";

my $count = 0;

$| = 1;
my $length;
while(my $one_seqid = <ONEIN>) {

    my $one_seq = <ONEIN>;
    my $one_qualid = <ONEIN>;
    my $one_qual = <ONEIN>;
    
    my $two_seqid = <TWOIN>;
    my $two_seq = <TWOIN>;
    my $two_qualid = <TWOIN>;
    my $two_qual = <TWOIN>;
    
    $length||=length($one_seq);
    
    if ((length($one_seq) != $length)||(length($two_seq) != $length)) {
        die "Error: a sequence in P1/P2 has different lengths";
    }
    
    chomp($one_seqid);
    chomp($one_seq);
    chomp($one_qualid);
    chomp($one_qual);
    chomp($two_seqid);
    chomp($two_seq);
    chomp($two_qualid);
    chomp($two_qual);

    print OUT $one_seqid.$sep.$two_seqid.
              "\n".
              $one_seq.$two_seq.
              "\n".
              $one_qualid.$sep.$two_qualid.
              "\n".
              $one_qual.$two_qual,"\n";
    

    $count++;
    if ($count % 10000 == 0) {
        print STDERR "Records: $count                                                         \r";
    }

}

print STDERR "Records: $count                                                         \n";

close(ONEIN);
close(TWOIN);

close(OUT);

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2010 Regional Blood Center of Ribeirão Preto

Usage

        $0	[-h/--help] [-i1 PEExp_1.fastq -i2 PEExp_2.fastq 
                         -o PEExp_12.fastq ] [-s "|"]

Argument(s)

        -h      --help          Help
        
        -i1     --inputfile1    Input file 1 (Ordered fastq p1 file - same order as p2)
        -i2     --inputfile2    Input file 2 (Ordered fastq p2 file - same order as p1)

        -o      --outputfile    Output file if fastq

        -s      --separator     Separator for ids (default: TAB)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

