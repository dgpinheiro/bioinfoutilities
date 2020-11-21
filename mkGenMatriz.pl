#!/usr/bin/env perl
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
#  Copyright (C) 2019  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
#  Faculdade de Ciências Agrárias e Veterinárias (FCAV)
#  Laboratório de Bioinformática (LB)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://www.fcav.unesp.br
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

Copyright (C) 2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

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

my ($level, %defines, $mapfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "define=s" => \%defines,
            "m|mapfile=s"=>\$mapfile
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
    
$LOGGER->logdie("Missing map file") unless ($mapfile);
$LOGGER->logdie("Wrong map file ($mapfile)") unless (-e $mapfile);


my %map;
open(MAP, "<", $mapfile) or $LOGGER->logdie($!);
while(<MAP>) {
    chomp;
    next if ($_=~/^#/);
    my ($id, $lane) = split(/\t/, $_);
    $map{$lane}=$id;
}
close(MAP);


$LOGGER->logdie("Missing define input file(s)") unless (scalar(keys %defines)>0);

my %data;

my %order_mw;

my @order_defines;

foreach my $k (keys %defines) {
    
    my $infile = $defines{$k};
    
    push(@order_defines, $k);

    $LOGGER->logdie("Wrong $k input file ($infile)") unless (-e $infile);


    open(IN, "<". $infile) or $LOGGER->logdie($!);
    while(<IN>) {
        chomp;
        next if (($_=~/^#/)||($_=~/^Lane #/));
        my ($lane, $band, $rf, $raw_vol, $cal_vol, $mw) = split(/\t/, $_);
        $lane=~s/\.$//;
        $lane+=0;
        $rf=~s/,/\./g;
        $rf+=0;
        $LOGGER->logdie("Not found lane ($lane) in map") unless (defined $map{$lane});
        
        next if ($map{$lane}=~/^MW|NA/);

        next if (($rf > 0.9)||
                ($rf < 0.1));
        
        $data{$k}->{$mw}->{$map{$lane}} = $raw_vol;
    }
    close(IN);

    @{$order_mw{$k}} = sort { $a <=> $b } keys %{$data{$k}}; 
}

my @order_lane = sort { $a <=> $b } keys %map; 

print 'ID';
foreach my $k (@order_defines) {
    foreach my $mw (@{$order_mw{$k}}) {
       print "\t", $k.'_'.$mw;
    }       
}
print "\n";

foreach my $lane (@order_lane) {
    next if ($map{$lane}=~/^MW|NA/);
    my @v;
    foreach my $k (@order_defines) {
        foreach my $mw (@{$order_mw{$k}}) {
            push(@v, $data{$k}->{$mw}->{$map{$lane}}||0);
        }
    }
    print join("\t", $map{$lane}, @v),"\n";
}

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -d      --define    Define Input file(s) (GelAnalyzer)
        -m      --mapfile   Map file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

