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

    dedup-ends.pl - Remove duplicated sequences from structured file.

=head1 SYNOPSIS
    
    # remove duplicated sequences from structured file PEExp_12.txt with this separator "|" (-s)
    # (P1_seqid|P2_seqid|sequence|P1_qualid|P2_qualid|qual)
    # generating PEExp_12_unique.fastq.
    # Obs.: The lines of file (PEExp_12.txt) MUST BE sorted by sequence.

    $ dedup-ends.pl -i PEExp_12.txt \
    > -o PEExp_12_unique.fastq -s "|";

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Remove duplicated sequences from a structured file with an specific separator (<sep>)
    containing these data P1_seqid<sep>P2_seqid<sep>sequence<sep>P1_qualid<sep>P2_qualid<sep>qual
    The lines of input file MUST BE sorted by sequence.

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

my ($outfile, $infile, $sep);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"  => sub { &Usage(); },
            "o|outputfile=s"=> \$outfile,
            "i|inputfile=s"=> \$infile,
            "s|separator=s"=> \$sep
) or &Usage();

$sep||="\t";

die "Missing output file" unless ($outfile);

die "Missing input file" unless ($infile);
die "Wrong input file" unless (-e $infile);

open(OUT, ">", $outfile) or die "Cannot open output file: $!";

open(IN, "<", $infile) or die "Cannot open input file: $!";

my $count = 0;

$| = 1;
my $length;
my ($first_seqids,$first_qualids);
my ($last_seqs);
my (@last_quals);

{    
    my $refseq = <IN>;
    chomp($refseq);

    my ($one_seqid,$two_seqid,$seqs,$one_qualid,$two_qualid,$quals) = split(quotemeta($sep), $refseq);
    
    # FIRST LINE
    $first_seqids = join($sep, $one_seqid,$two_seqid);
    $first_qualids = join($sep, $one_qualid,$two_qualid);
    $last_seqs = $seqs;
    @last_quals = split(//, $quals);
}

while(my $refseq = <IN>) {
    chomp($refseq);

    my ($one_seqid,$two_seqid,$seqs,$one_qualid,$two_qualid,$quals) = split(quotemeta($sep), $refseq);
    $length||=length($seqs);
    
    if ((length($seqs) != $length)||(length($quals) != $length)) {
        die "Error: a sequence in file has different lengths";
    }
    
    if ($seqs eq $last_seqs) {
        my @q = split(//, $quals);
        for (my $i = 0; $i<length($quals); $i++) {
            if (ord($last_quals[$i]) < ord($q[$i])) {
                $last_quals[$i] = $q[$i];
            }
        }
    }
    else {
        print OUT join($sep, $first_seqids, $last_seqs, $first_qualids, join('',@last_quals)),"\n";

        $first_seqids = join($sep, $one_seqid,$two_seqid);
        $first_qualids = join($sep, $one_qualid,$two_qualid);
        $last_seqs = $seqs;
        @last_quals = split(//, $quals);
    }
    
    $count++;
    if ($count % 10000 == 0) {
        print STDERR "Records: $count                                                         \r";
    }
}

print OUT join($sep, $first_seqids, $last_seqs, $first_qualids, join('',@last_quals)),"\n";

print STDERR "Records: $count                                                         \n";

close(OUT);

close(IN);

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2010 Regional Blood Center of Ribeirão Preto

Usage

        $0	[-h/--help] [-i PEExp_12.txt 
                         -o PEExp_12_unique.txt ] [-s "|"]

Argument(s)

        -h      --help          Help
        
        -o      --outputfile    Output file 

        -i      --inputfile     Input file (from merge-ends.pl)

        -s      --separator     Separator (default: TAB)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

