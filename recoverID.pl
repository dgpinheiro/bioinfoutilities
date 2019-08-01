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

my ($level, $infile, $outfile, $hashfile, $fixedlen, $padchar, $sub, $replace);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "o|outfile=s"=>\$outfile,
            "h|hashfile=s"=>\$hashfile,
            "f|fixedlen=i"=>\$fixedlen,
            "c|char=s"=>\$padchar,
            "s|sub=s"=>\$sub,
            "r|replace=s"=>\$replace
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

use Storable;

$padchar||=" ";

$LOGGER->logdie("Missing hash file") unless ($hashfile);
$LOGGER->logdie("Wrong hash file ($hashfile)") unless (-e $hashfile);

my $hr_orig = retrieve($hashfile);

#foreach my $k (keys %{$hr_orig}) {
#    print $k,"\t",$hr_orig->{$k},"\n";
#}

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

$LOGGER->logdie("Missing output file") unless ($outfile);

open(IN, "<", $infile) or $LOGGER->logdie($!);

open(OUT, ">", $outfile) or $LOGGER->logdie($!);

while(my $l = <IN>) {
    chomp($l);
    foreach my $k ( keys %{ $hr_orig } ) {
        my $x = $hr_orig->{$k};
        if ($fixedlen) {
            if (length($x) < $fixedlen) {
                $x=&RPad($x, $fixedlen, $padchar);
            }
        }

        if (($sub)&&($replace)) {
            $x=~s/$sub/$replace/g;
        }
	$k=quotemeta($k);
        $l=~s/$k/$x/g;
    }
    print OUT $l,"\n";
}

close(IN);
close(OUT);

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
        -o      --outfile   Output file
        -h      --hashfile  Hash file
        -f      --fiexedlen Fixed length of strings (right pad character)
        -c      --char      Character for right padding [Default space]
        -s      --search    Search regular expression (used with -s)
        -r      --replace   Replace string (used with -s)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub RPad {

    my($str, $len, $chr) = @_;

    $chr = " " unless (defined($chr));
        
    return substr($str . ($chr x $len), 0, $len);

} # RPad

