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

my ($level, $lfc_threshold, $qvalue_threshold, $infile, $cmp);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "q|qvalue=f"=>\$qvalue_threshold,
            "f|lfc=f"=>\$lfc_threshold,
            "i|infile=s"=>\$infile,
            "c|cmp=s"=>\$cmp
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

$qvalue_threshold||=1;
$lfc_threshold||=0;

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

open(IN, "<", $infile) or $LOGGER->logdie($!);
my $header_string = <IN>;
chomp($header_string);
my @header=split(/\t/, $header_string);

my @analysis;
foreach my $h (@header) {
    if ($h =~/^(\S{3}x\S{3})\.fold\.change$/) {
        my $an=$1;
        if ($cmp) {
            if ($cmp eq $an) {
                push(@analysis, $an);
            }
        } else {
            push(@analysis, $an);
        }
    }
}

$LOGGER->logdie("Not found analysis") unless (scalar(@analysis));


print join("\t", @header),"\n";
while(<IN>) {
    chomp;
    my %data;
    @data{ @header } = split(/\t/, $_);
    my ($set_qvalue, $set_lfc);
    foreach my $an (@analysis) {
            if ($data{"$an.qvalue"} <= $qvalue_threshold) {
                $set_qvalue = 1;
            }
            if (abs($data{"$an.fold.change"}) >= $lfc_threshold) {
                $set_lfc  = 1;
            }
    }
    if (($set_qvalue)&&($set_lfc)) {
        print join("\t", @data{ @header }),"\n";
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
        -i      --infile    Input file (output from getDEMresult.pl)
        -q      --qvalue    q-value threshold [Default: 1]
        -f      --lfc       abs-log-fold-change threshold [Default: 0]
        
END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

