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
    Log::Log4perl->easy_init($WARN);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $indir, $minimum_free_energy_cutoff, $pvalue_cutoff, $nogu, $seed);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"  => sub { &Usage(); },
            "l|level=s" => \$level,
            "i|indir=s" => \$indir,
            "e|mfe=f"   =>\$minimum_free_energy_cutoff,
            "p|pvalue=f"=>\$pvalue_cutoff,
            "g|nogu"    =>\$nogu,
            "s|seed=s"  =>\$seed
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

$LOGGER->logdie("Missing input dir") unless ($indir);
$LOGGER->logdie("Wrong input dir ($indir)") unless (-d $indir);

my @seedcoord;
if ($seed) {
    $seed=~s/\s//g;
    if ($seed=~/^(\d+),(\d+)$/) {
        @seedcoord = ($1,$2);
        $LOGGER->logdie("Wrong seed $seed") if (($seedcoord[0] <= 0)||
                                                ($seedcoord[1] <= 0)||
                                                ($seedcoord[0] > $seedcoord[1]));
    }
    else {
        $LOGGER->logdie("Wrong seed coordinates");
    }
}

if ($nogu) {
    $LOGGER->logdie("Missing seed coordinate") unless ($seed);
}

$indir=~s/\/+$//;

$pvalue_cutoff||= 1;
$minimum_free_energy_cutoff||= -1;

while(my $f=glob("$indir/*.txt")) {
    unless (-d $f) {
        #print $f,"\n";
        open(IN, "<", $f) or $LOGGER->logdie($!);
        while(my $l=<IN>) {
            chomp($l);
            #target  name, query name, minimum free energy, position in target, alignment line 1, line 2, line 3, line 4
            my ($target, $target_length, $query_name, $query_length, $minimum_free_energy, $pvalue, $position_in_target,$alignment_line_1,$alignment_line_2,$alignment_line_3,$alignment_line_4, $extra) = split(/:/, $l);
            
#                    print $target, "\t", $query_name, "\t", $minimum_free_energy, "\t", $pvalue, "\n";
            if ($pvalue <= $pvalue_cutoff) {
                if ( $minimum_free_energy <= $minimum_free_energy_cutoff ) {
                    my $noguset;
                    if ($nogu) {
                        my $ral2 = reverse($alignment_line_2);
                        my $ral3 = reverse($alignment_line_3);
                        my @al2 = split(//, $ral2);
                        my @al3 = split(//, $ral3);
                        for (my $i = ($seedcoord[0]-1); $i<=($seedcoord[1]-1); $i++) {
                            if ( ($al2[$i] =~ /^\S+$/) && ($al3[$i] =~ /^\S+$/)) {
                                if ( (($al2[$i] eq 'G')&&($al3[$i] eq 'U')) ||
                                     (($al2[$i] eq 'U')&&($al3[$i] eq 'G')) ) {
                                    $noguset = 1;
                                    last;
                                }
                            }
                            else {
                                $LOGGER->logwarn("There isn't a perfect alignment in the seed position ($seed) in the alignment between $target and $query_name ($l)");
                                $noguset = 1;
                                last;
                            }                            
                        }
                    }

                    unless ($noguset) {

                        print $target, "\t", $query_name, "\t", $minimum_free_energy, "\t", $pvalue, "\t", $position_in_target, "\n";

                        #print $target,"\t",$target_length,"\n";
                        #print $query_name,"\t",$query_length,"\n";
                        #print "ENERGY: ",$minimum_free_energy,"\n";
                        #print "PVALUE: ",$pvalue,"\n";
                        #print "POS: ",$position_in_target,"\n";
                        #print $alignment_line_1,"\n";
                        #print $alignment_line_2,"\n";
                        #print $alignment_line_3,"\n";
                        #print $alignment_line_4,"\n";
                    } 
                } 
            }
        }
        close(IN);
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
        -i      --indir     Input directory (output of RNAhybrid)
        -e      --mfe       Minimum free energy cutoff (default: -1)
        -p      --pvaluef   p-value cutoff (default: 1)
        -g      --noGU      no GU alignment in the seed
        -s      --seed      Seed coordinates

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

