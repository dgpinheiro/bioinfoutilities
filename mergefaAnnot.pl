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

my ($level, $fafile, $tabfile, $colname, $newname);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "f|fa=s"=>\$fafile,
            "t|tab=s"=>\$tabfile,
            "c|column=s"=>\$colname,
            "n|nanme=s"=>\$newname
    ) or &Usage();

$newname||='Annotation';

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

$LOGGER->logdie("Missing fasta file") unless ($fafile);
$LOGGER->logdie("Wrong fasta file ($fafile)") unless (-e $fafile);

$LOGGER->logdie("Missing tabular file") unless ($tabfile);
$LOGGER->logdie("Wrong tabular file ($tabfile)") unless (-e $tabfile);

$LOGGER->logdie("Missing column name") unless ($colname);

my %annot;

open(FA, "<", $fafile) or $LOGGER->logdie($!);
while(<FA>) {
    chomp;
    if ($_=~/^>(\S+)\s+(.*)/) {
        $LOGGER->logdie("Duplicated identifier in fasta file ($1)") if (exists $annot{$1});
        $annot{$1} = $2;
    }
}
close(FA);

open(TAB, "<", $tabfile) or $LOGGER->logdie($!);
my $header_line=<TAB>;
chomp($header_line);
my @header=split(/\t/, $header_line);

{
    my %check;
    @check{ @header } = (1..$#header);

    $LOGGER->logdie("The column name already exists in tabular file") if (exists $check{ $newname });
    $LOGGER->logdie("Not found column $colname in tabular file") unless (exists $check{ $colname });
}

my @newheader;
foreach my $h (@header) {
    push(@newheader, $h);
    if ( $h eq $colname ) {
        push(@newheader, $newname);
    }
}

print join("\t", @newheader),"\n";

while(<TAB>) {
    chomp;
    my %data;
    @data{ @header } = split(/\t/, $_);
    $data{ $newname } = $annot{ $data{ $colname } }||'';
    print join("\t", @data{ @newheader }),"\n";
}
close(TAB);

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
        -f      --fa    Annotated fasta file
        -t      --tab   Tabular file (count matrix)
        -c      --col   ID column name
        -n      --name  Annotation column name [Default: Annotation]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

