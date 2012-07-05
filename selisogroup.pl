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

    selisogroup.pl  - Select sequences of one specified isogroup number from 454Isotigs.fna.

=head1 SYNOPSIS

    # Select only sequences of isogroup00005 and redirects to STDOUT

    $ selisogroup.pl   -f 454Isogroups.fna -g 5

=head1 ABSTRACT

=head1 DESCRIPTION

    Script to select sequences of one specified isogroup number from 454Isotigs.fna in fasta format.

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

my ($isogroupn,$file,$onlyisotigs,$onlycontigs);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "g|isogroup=s"=>\$isogroupn,
            "i|isotigs"=>\$onlyisotigs,
            "c|contigs"=>\$onlycontigs,
            "f|file=s"=>\$file
) or &Usage();

die "Missing FASTA file" unless ($file);
die "Wrong FASTA file ($file)" unless (-e $file);

die "Missing isogroup number" unless (defined $isogroupn);

my $isogroup = "isogroup".&LPad($isogroupn, 5, 0);

my $extra = '(?:';
if ($onlyisotigs) {
    $extra.= 'isotig';
    if ($onlycontigs) {
        $extra.='|contig';
    }        
}
else {
    if ($onlycontigs) {
        $extra.='contig';
    }
    else {
        $extra.='\S+';
    }
}
$extra.=')\d+';

#print '^>'.$extra.'\s+gene='.$isogroup."\n";
my $set;
open(IN, "<", $file) or die $!;
while(<IN>) {
    chomp;
    if ($_=~/^>/) { 
        if ($_=~/^>$extra\s+gene=$isogroup/) {
            $set = 1;
        }
        else {
            $set = undef;
        }
    }        
    if ($set) {
        print $_,"\n";
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

        $0	[-h/--help] [-f 454Isotigs.fna -g 5] [-i] [-c]

Argument(s)

        -h      --help      Help
        -g      --isogroup  Isogroup (e.g. 5)
        -i      --isotigs   Select isotigs
        -c      --contigs   Select contigs
        -f      --file      File 454Isotigs.fna
        
END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}



#---------------------------------------------------------------------
# LPad
#---------------------------------------------------------------------
# Pads a string on the left end to a specified length with a specified
# character and returns the result.  Default pad char is space.
#---------------------------------------------------------------------

sub LPad {

    my ($str, $len, $chr) = @_;

    $chr = " " unless (defined($chr));
        
    return substr(($chr x $len) . $str, -1 * $len, $len);

} # LPad

