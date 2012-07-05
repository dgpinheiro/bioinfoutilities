#!/usr/bin/perl 
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

    refGene2annot.pl - Script to convert the refGene.txt file to Goby annotation file.

=head1 SYNOPSIS

    # converts the refGene.txt file to Goby annotation file and prints to STDOUT, which was 
    # redirected to goby_gene_annotation.txt

    $ refGene2annot.pl -i refGene.txt > goby_gene_annotation.txt

=head1 ABSTRACT

=head1 DESCRIPTION

    Script to convert the refGene.txt file, downloaded from UCSC GoldenPath, to 
    annotation file used by Goby (http://campagnelab.org/software/goby).

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


my ($input);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "i|input=s"=>\$input
    ) or &Usage();

die "Missing refGene file" unless ($input);    
die "Wrong refGene file ($input)" unless (-e $input);    

open(IN, "<", $input) or die $!;

my @refGeneCols = ( 'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'id', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames' );

my @annotCols = ('chrom', 'strand', 'name2', 'exonName', 'exonStart', 'exonEnd');

my %refGene;

print join("\t", @annotCols),"\n";
while(<IN>) {
    chomp;
    @refGene{@refGeneCols} = split(/\t/, $_);
    
    $refGene{'strand'} = (($refGene{'strand'} eq '+' ) ? 1 : -1);

    chop($refGene{'exonStarts'});
    chop($refGene{'exonEnds'});
    
    my @exonStarts = split(',', $refGene{'exonStarts'});
    my @exonEnds = split(',', $refGene{'exonEnds'});
 
    my $exonStartsCount = scalar(@exonStarts);
    my $exonEndsCount = scalar(@exonEnds);

    die "Wrong exonStartsCount $exonStartsCount/$refGene{'exonCount'}" unless ($exonStartsCount == $refGene{'exonCount'});
    die "Wrong exonEndsCount $exonEndsCount/$refGene{'exonCount'}" unless ($exonEndsCount == $refGene{'exonCount'});
    
    for (my $i = 0; $i < $refGene{'exonCount'}; $i++) {
        my $exonCount = $i+1;
        print join("\t", @refGene{'chrom','strand','name2'}, $refGene{'name'}.'_e'.$exonCount, $exonStarts[$i], $exonEnds[$i] ),"\n";
    }
}

close(IN);

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2010 Regional Blood Center of Ribeirão Preto

Usage

        $0	[-h/--help] [-i refGene.txt]

Argument(s)

        -h      --help  Help
        -i      --input refGene.txt (e.g. http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

