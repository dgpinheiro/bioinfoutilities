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
#use utf8;
#use Encode;
#use Encode::Guess;

binmode(STDERR, ":utf8");
binmode(STDOUT, ":utf8");

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

use Spreadsheet::WriteExcel;

my ($level, $infile, $outfile, $hasheader);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "o|outfile=s"=>\$outfile,
            "h|header"=>\$hasheader
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

$outfile||='./output';

$outfile.='.xls';

$LOGGER->logdie("Missing intput file") unless ($infile);
$LOGGER->logdie("Wrong intput file ($infile)") unless (-e $infile);



# Create a new Excel workbook
my $workbook = Spreadsheet::WriteExcel->new($outfile);

$workbook->set_properties(utf8 => 1);

# Add a worksheet
my $worksheet = $workbook->add_worksheet();

#  Add and define a format
my $header_format = $workbook->add_format(); # Add a format
$header_format->set_bold();
$header_format->set_bg_color('cyan');
$header_format->set_align('center');
$header_format->set_bold();

# Write a formatted and unformatted string, row and column notation.
my $col = 0;
my $row = 0;

open(IN, '<:encoding(utf8)',$infile) or die $!;

if ($hasheader) {
    my $header = <IN>;
    chomp($header);
    my @head=split(/\t/, $header);

    $worksheet->set_column(0,$#head,30);

    foreach my $h (@head) {
        $worksheet->write($row, $col, $h, $header_format);
        $col++;
    }
	$row++;
}

while(<IN>) {
	chomp;
	my @data = split(/\t/, $_);
	
	$col=0;
	foreach my $d (@data) {
		if ($d =~/^[\+\-]?\d+(?:\.\d+)?$/) {
			$worksheet->write($row, $col, $d+0);
		} else {
			$worksheet->write($row, $col, $d);
		}
		$col++;
	}
	$row++;
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
        -i      --infile    Input file (Tab-separated-values format)
        -o      --outfile   Output file basename [Default: ./output.xls]
        -h      --header    Input file has header line

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

