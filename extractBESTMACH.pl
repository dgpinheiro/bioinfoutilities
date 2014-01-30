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

my ($level, $input, $reffiles, $threshold, $opt, $sbhfile, $auxfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|input=s"=>\$input,
            "r|ref=s"=>\$reffiles,
            "t|threshold=f"=>\$threshold,
            "o|opt=s"=>\$opt,
            "b|sbh=s"=>\$sbhfile,
            "a|aux=s"=>\$auxfile
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

$LOGGER->logdie("Missing reference file (subject)") unless ($reffiles);
$LOGGER->logdie("Wrong reference file (subject)") unless (-e $reffiles);

$LOGGER->logdie("Missing input file") unless ($input);
$LOGGER->logdie("Wrong input file") unless (-e $input);

my %conversor;
if ($auxfile) {
    open(AUX, "<", $auxfile) or $LOGGER->logdie($!);
    while(<AUX>) {
        chomp;
        my ($to, $from) = split(/\t/, $_);
        $conversor{$from} = $to;
    }
    close(AUX);
}

my %sbh;

if ($sbhfile) {
    open(SBH, "<", $sbhfile) or $LOGGER->logdie($!);

    while(<SBH>) {
        chomp;
        if ($_=~/Query=(\S+) Hit=(\S+) Score=(\S+) /) {
            my ($q, $h, $s) = ($1, $2, $3);
            $sbh{$q}->{ ((exists $conversor{$h}) ? $conversor{$h} : $h)   } = $s;
        }
    }

    close(SBH);
}


my %refscores;

open(REF, "<", $reffiles) or $LOGGER->logdie($!);

while(<REF>) {
    chomp;
    if ($_=~/Query=(\S+) Hit=(\S+) Score=(\S+) /) {
        my ($q, $h, $s) = ($1, $2, $3);
        if ($q ne $h) {
            $LOGGER->logdie("Wrong association found $_");
        }
        $refscores{$q} = $s;
    }
}

close(REF);

use File::Basename;

my ($b1, $b2) = split('_X_', basename($input, '.bls'));

#print $b1,"\t",$b2,"\n";

$threshold||=0.9;

use Bio::SearchIO; 

my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $input);

while ( my $result = $in->next_result ) {
    ## $result is a Bio::Search::Result::ResultI compliant object
    my $best_hit_score;
    my @allhits;

    while ( my $hit = $result->next_hit ) {
        $best_hit_score ||= $hit->bits();
        if ($hit->bits() > $best_hit_score) {
            $best_hit_score = $hit->bits();
        }
        push(@allhits, $hit);
    }
 
    if ($opt eq 'query') {   
        $LOGGER->logdie( "Not found refscoreq for " . $result->query_name() ) unless ( $refscores{ $result->query_name() } );
    }
    
    foreach my $hit (@allhits) {
    
        my $refval=0;
        if ($opt eq 'hit') {
            $LOGGER->logdie( "Not found refscoreq for " . $hit->name() ) unless ( $refscores{ $hit->name() } );
            $refval = $refscores{$hit->name()}
        }
        elsif ($opt eq 'query') {
            $refval = $refscores{$result->query_name()};
        }
        else {
            $refval = 'NA';
        }

        if ( abs( $hit->bits() / $best_hit_score ) >= $threshold ) {

              print
                           "Query=",   $result->query_name,
                           " Hit=",    $hit->name,
                           " Score=",  $hit->bits,
                           " Evalue=", $hit->significance,
                           " Threshold=", $threshold,
                           " Ratio=", sprintf( "%.2f", abs( $hit->bits() / $best_hit_score )),
                           " Reference(s) ($opt) =", $refval,
                           ((exists $sbh{$hit->name()}->{$result->query_name()}) ? "\t***RBH*** (Subject Best Hit Score=".$sbh{$hit->name()}->{$result->query_name()}.")" : ''),
                           "\n";

        }
        else {
            print STDERR "INSUFFICIENT", "\t",
              "Query=",      $result->query_name,
              " Hit=",       $hit->name,
              " Score=",     $hit->bits,
              " Evalue=",    $hit->significance,
              " Threshold=", $threshold,
              " Ratio=", sprintf( "%.2f", abs( $hit->bits() / $best_hit_score )),
              " Reference(s) ($opt) =", $refval,
              "\n";
        }

    }
}


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
        -i      --input     Input
        -r      --ref       Reference score file
        -t      --threshold Threshold   [Default: 0.9]
        -o      --opt       opt [Default: query]
        -b      --sbh       subject best hit
        -a      --aux       Aux

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

