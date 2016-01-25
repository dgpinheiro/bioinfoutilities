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

my ($level, $simplefile, $tabfile, $colid, $colann, $newname, $hasheader, $colkey);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|simple=s"=>\$simplefile,
            "t|tab=s"=>\$tabfile,
            "c|colid=i"=>\$colid,
            "k|colkey=s"=>\$colkey,
            "a|colann=i"=>\$colann,
            "n|nanme=s"=>\$newname,
            "e|hasheader"=>\$hasheader
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


$LOGGER->logdie("Missing simple tsv file") unless ($simplefile);
$LOGGER->logdie("Wrong simple tsv file ($simplefile)") unless (-e $simplefile);

$LOGGER->logdie("Missing tabular file") unless ($tabfile);
$LOGGER->logdie("Wrong tabular file ($tabfile)") unless (-e $tabfile);

$colid||=1;
$colann||=2;

$newname||="Annotation";

open(TAB, "<", $tabfile) or $LOGGER->logdie($!);
my $header_line=<TAB>;
chomp($header_line);
my @header=split(/\t/, $header_line);

{
    my %check;
    @check{ @header } = (1..$#header);

    $LOGGER->logdie("The column name already exists in tabular file") if (exists $check{ $newname });

    $LOGGER->logdie("Not found column $colkey in tabular file") unless (exists $check{ $colkey });
}

my @newheader;
foreach my $h (@header) {
    push(@newheader, $h);
    if ( $h eq $colkey ) {
        push(@newheader, $newname);
    }
}

my %annot;

my @res;
while(<TAB>) {
    chomp;
    my %data;
    @data{ @header } = split(/\t/, $_);
    push(@res, \%data);
    $annot{ $data{ $colkey } } = undef;
}
close(TAB);


open(SIMPLE, "<", $simplefile) or $LOGGER->logdie($!);
if ($hasheader) {
    my $h = <SIMPLE>;
}

while(<SIMPLE>) {
    chomp;
    my (@simple) = split(/\t/, $_);
    if ( exists $annot{ $simple[$colid-1] } ) {
        $annot{ $simple[$colid-1] } = $simple[$colann-1];
        #print STDERR "(",$simple[$colid-1],")","\t",$simple[$colann-1],"\n";
    }
}
close(SIMPLE);

print join("\t", @newheader),"\n";
foreach my $hr_data ( @res ) {

    $hr_data->{ $newname } = $annot{ $hr_data->{ $colkey } }||'';

    print join("\t", @{$hr_data}{ @newheader }),"\n";
}

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help  Help
        -l      --level Log level [Default: FATAL]
        -t      --tab       Tabular file (count matrix)
        -k      --colkey    Column name in Tabular file to be used as key to merge
        -i      --simple    Simple annotation file
        -e      --hasheader Simple annotation file has header [Default: off]
        -c      --colid     Column number for ID (simple annotation file) [Default: 1]
        -a      --colann    Column number for Annotation (simple annotation file) [Default: 2]
        -n      --name      Annotation column name in output [Default: Annotation]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

