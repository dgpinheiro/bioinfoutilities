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
# $Id$

=head1 NAME

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Arguments:

        -h/--help   Help
        -l/--level  Log level [Default: FATAL] 
            OFF
            FATAL
            ERROR
            WARN
            INFO
            DEBUG
            TRACE
            ALL

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2012 Universidade de São Paulo

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Readonly;
use Getopt::Long;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $gffile, $col, $regexp);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "g|gffile=s"=>\$gffile,
            "c|col=i"=>\$col,
            "r|regexp=s"=>\$regexp
    ) or &Usage();


if ($level) {
    my %LEVEL = (   
    'OFF'   =>$OFF,
    'FATAL' =>$FATAL,
    'ERROR' =>$ERROR,
    'WARN'  =>$WARN,
    'INFO'  =>$INFO,
    'DEBUG' =>$DEBUG,
    'TRACE' =>$TRACE,
    'ALL'   =>$ALL);
    $LOGGER->logdie("Wrong log level ($level). Choose one of: ".join(', ', keys %LEVEL)) unless (exists $LEVEL{$level});
    Log::Log4perl->easy_init($LEVEL{$level});
}

$col||=1;

my %gene;
my %parent;

$regexp||='^(\S+)';

# /data/tmp/files/rnas/LBDAv0.2.gff
open(GFF, "<", $gffile) or $LOGGER->logdie($!);
while(<GFF>) {
    next if ($_=~/^#/);
    chomp;

    my @d=split(/\t/, $_);
    if ($d[2] eq 'gene') {
        my ($ID) = $d[8] =~ /ID=([^;]+)/;
        $LOGGER->logdie("Missing ID ($d[8])") unless ($ID);
        $gene{$ID} = undef;
    } else {

        next if ($d[2] =~/exon|CDS|five_prime_UTR|three_prime_UTR/);

        my ($ID) = $d[8] =~ /ID=([^;]+)/;
        my ($Parent) = $d[8] =~ /Parent=([^;]+)/;

        $LOGGER->logdie("Missing ID ($d[8])") unless ($ID);
        $LOGGER->logdie("Missing Parent ($d[8])") unless ($Parent);

        if (exists $gene{$Parent}) {
            $gene{$Parent}->{$ID}=undef;
            $parent{$ID} = $Parent;
        } else {
            print "NOT FOUND\t$ID\t$Parent\n";
        }

    }
}
close(GFF);

open(IN, "<", $infile) or $LOGGER->logdie($!);
my @header;

while(<IN>) {
    chomp;
    if ($.==1) {
        @header=split(/\t/, $_);
        print join("\t", 'GENE', @header),"\n";
    } else {
        my %data;
        @data{ @header } = split(/\t/, $_);
        
        if ($data{ $header[ $col-1 ] } =~ /$regexp/) {
            my $ID=$1;
            my $gene = $parent{$ID};
            unless ($gene) {
               foreach my $g (keys %gene) {
                   foreach my $c (keys %{ $gene{$g} }) {
                        if ($c=~/$ID/) {
                            $gene=$g;
                            last;
                        }
                    }
                    last if ($gene);
                }
            }
            $LOGGER->logdie("Not found parent of $ID") unless ($gene);
            print join("\t", $gene, @data{ @header }),"\n";
        }
    }        
}

close(IN);

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file
        -c      --col       Column number [Default: 1]
        -g      --gffile    GFF
        -r      --regexp    Regular expression to capture from infile column selected

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

