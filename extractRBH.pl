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

my ($level, $infile, $minscore, $maxevalue, $rbh, $other);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=> \$infile,
            "s|minscore=i"=>\$minscore,
            "e|maxevalue=s"=>\$maxevalue,
            "r|rbh"=>\$rbh,
            "o|other"=>\$other
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

$minscore||= 0;
$maxevalue||= 100;

$minscore+=0;
$maxevalue+=0;

my %RBH_GB;
my %RBH_FB;
my %OTHER_GB;
my %OTHER_FB;

open(IN, "<", $infile) or $LOGGER->logdie($!);
while(<IN>) {
    chomp;
    my ($query);
    if ($_=~/Query=gnl\|Amel_4\.5\|([^- ]+)/) {
	    ($query)= $1;
    } elsif ($_=~/Query=([^- ]+)/) {
    	    ($query)= $1;
    }
    my ($hit)=$_=~/Hit=([^ ]+)/;
    my ($evalue)=$_=~/Evalue=([^ ]+)/;
    my ($score)=$_=~/Score=(\d+)/;

    if ($_=~/\*\*\*RBH\*\*\*/) {
        $RBH_GB{$query}=$hit;
        $RBH_FB{$hit}=$query;
    } else {
        if ( ($score >= $minscore) && ($evalue <= $maxevalue) ) {
            $OTHER_GB{$query}->{$hit} = {'score'=>$score, 'evalue'=>$evalue};
            $OTHER_FB{$hit}->{$query} = { 'score'=>$score, 'evalue'=>$evalue};
        }            
    }
}
close(IN);

my %already;

foreach my $gb (keys %RBH_GB) {
    print $gb,"\t",$RBH_GB{$gb},"\n" if ($rbh);
    $already{$gb}=undef;    
}

if ($other) {
    foreach my $gb (keys %OTHER_GB) {
        next if (exists $already{$gb});

        foreach my $hit (sort { $OTHER_GB{$gb}->{$b}->{'score'} <=> $OTHER_GB{$gb}->{$a}->{'score'} || $OTHER_GB{$gb}->{$a}->{'evalue'} <=> $OTHER_GB{$gb}->{$b}->{'evalue'} } keys %{ $OTHER_GB{$gb} }) {
            unless (exists $RBH_FB{ $hit }) {
                my $e;
                foreach my $tst (keys %{ $OTHER_FB{ $hit } }) {
                    next if ( $tst eq $gb );
                    $e=1;
                }

                unless ($e) {
                    print $gb,"\t",$hit,"\n";
                    last;
                }
            }
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

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -i      --infile        Input file
        -s      --minscore      Minimum score [Default: 0]
        -e      --maxevalue     Maximum evalue [Default: 100]
        -r      --rbh           RBH 
        -o      --other         Other

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

