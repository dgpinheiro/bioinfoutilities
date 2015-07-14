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

my ($level, $infile, $promoter_size, $outbase, $chrsizefile, %chrsize);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "p|promotersize=i"=>\$promoter_size,
            "o|outbase=s"=>\$outbase,
            "s|chrsizefile=s"=>\$chrsizefile
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

$outbase||='./coords';

$promoter_size||=1000;

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

$LOGGER->logdie("Missing chr size file") unless ($chrsizefile);
$LOGGER->logdie("Wrong chr size file ($chrsizefile)") unless (-e $chrsizefile);

open(CHR, "<", $chrsizefile) or $LOGGER->logdie($!);
while(<CHR>) {
    chomp;
    my ($chr, $size) = split(/\t/, $_);
    $chrsize{ $chr} = $size;
}
close(CHR);

my %data;

open(IN, "<", $infile) or $LOGGER->logdie($!);
while(<IN>) {
    chomp;
    my ($chr, $src, $feat, $start, $end, $score, $strand, $frame, $info) = split(/\t/ ,$_);

    next if (($feat ne 'exon')&&($feat ne 'CDS'));

    my ($transcript_id) = $info=~/transcript_id \"([^\"]+)\"/;
    $LOGGER->logdie("Not found transcript id ($_)") unless ($transcript_id);
    my ($gene_id) = $info=~/gene_id \"([^\"]+)\"/;
    $LOGGER->logdie("Not found gene id ($_)") unless ($gene_id);
    
    $data{$transcript_id}->{'gene_id'} = $gene_id if (! exists $data{$transcript_id}->{'gene_id'});
    $data{$transcript_id}->{'chr'} = $chr if (! exists $data{$transcript_id}->{'chr'});
    if (! exists $data{$transcript_id}->{'strand'}) {
        $data{$transcript_id}->{'strand'} = $strand;
    } else {
        $LOGGER->logdie("Different strand for $transcript_id") if ($strand ne $data{$transcript_id}->{'strand'});
    }
    
    if ($feat eq 'exon') {
        $data{$transcript_id}->{'eleft'} = $start+0 if ((! exists $data{$transcript_id}->{'eleft'}) || ($data{$transcript_id}->{'eleft'}>$start));
        $data{$transcript_id}->{'eright'} = $end+0 if ((! exists $data{$transcript_id}->{'eright'}) || ($data{$transcript_id}->{'eright'}<$end));
        push(@{ $data{$transcript_id}->{'exons'} }, [$start+0, $end+0]);
    }
    if ($feat eq 'CDS') {
        $data{$transcript_id}->{'cleft'} = $start+0 if ((! exists $data{$transcript_id}->{'cleft'}) || ($data{$transcript_id}->{'cleft'}>$start));
        $data{$transcript_id}->{'cright'} = $end+0 if ((! exists $data{$transcript_id}->{'cright'}) || ($data{$transcript_id}->{'cright'}<$end));

        push(@{ $data{$transcript_id}->{'cds'} }, [$start+0, $end+0]);
    }
}
close(IN);


open(TSS, ">", $outbase.'_tss.bed') or $LOGGER->logdie($!);
open(EXON, ">", $outbase.'_exon.bed') or $LOGGER->logdie($!);
open(CDS, ">", $outbase.'_cds.bed') or $LOGGER->logdie($!);
open(FIVE, ">", $outbase.'_5utr.bed') or $LOGGER->logdie($!);
open(THREE, ">", $outbase.'_3utr.bed') or $LOGGER->logdie($!);
open(INTRON, ">", $outbase.'_intron.bed') or $LOGGER->logdie($!);
open(PROMOTER, ">", $outbase.'_promoter.bed') or $LOGGER->logdie($!);

foreach my $transcript_id (sort { &compn($data{$a}->{'chr'}, $data{$b}->{'chr'}) } keys %data) {

    my $gene_id = $data{$transcript_id}->{'gene_id'};

    my ($tss);
    my ($promoter_start, $promoter_end);
    
    my @exons;
    my @cds;
    my @fivep;
    my @threep;

    if ($data{$transcript_id}->{'strand'} eq '+') {

        $tss = $data{$transcript_id}->{'eleft'};

        $promoter_end = $tss;

        if ($tss < $promoter_size) {
            $promoter_start = 0;
        } else {
            $promoter_start = $promoter_end-$promoter_size;
        }            

        foreach my $ar_e (sort { $a->[0] <=> $b->[0] } @{ $data{$transcript_id}->{'exons'} }) {
            push(@exons, [ @{ $ar_e } ]);

            if (exists $data{$transcript_id}->{'cds'}) {
                my ($start_fivep, $end_fivep);
                my ($start_threep, $end_threep);
                if ($ar_e->[0] < $data{$transcript_id}->{'cleft'}) {
                    
                    $start_fivep = $ar_e->[0];
                    
                    if ($ar_e->[1] < $data{$transcript_id}->{'cleft'}) {
                        $end_fivep = $ar_e->[1];
                    } else {
                        $end_fivep = $data{$transcript_id}->{'cleft'};
                    }
                    
                    push(@fivep, [$start_fivep, $end_fivep]);
                    
                } 

                if ( ( $data{$transcript_id}->{'cright'} > $ar_e->[0] ) && 
                     ( $data{$transcript_id}->{'cright'} < $ar_e->[1] ) ) {
                    
                     $start_threep = $data{$transcript_id}->{'cright'};
                     $end_threep = $ar_e->[1];

                     push(@threep, [$start_threep, $end_threep]);

                } elsif ( $ar_e->[0] > $data{$transcript_id}->{'cright'} ) {
                     
                     $start_threep = $ar_e->[0];
                     $end_threep = $ar_e->[1];
                    
                     push(@threep, [$start_threep, $end_threep]);

                }
            }
        }

        if (exists $data{$transcript_id}->{'cds'}) {
            foreach my $ar_c (sort { $a->[0] <=> $b->[0] } @{ $data{$transcript_id}->{'cds'} }) {
                push(@cds, [ @{ $ar_c } ]);
            }
        }

    } elsif ($data{$transcript_id}->{'strand'} eq '-') {

        $tss = $data{$transcript_id}->{'eright'}-1;

        $promoter_start = $tss+1;
        if ( ($tss+$promoter_size) > $chrsize{ $data{$transcript_id}->{'chr'} }) {
            $promoter_end = $chrsize{ $data{$transcript_id}->{'chr'} };
        } else {
            $promoter_end = $promoter_start+$promoter_size;
        }            

        foreach my $ar_e (sort {$b->[1] <=> $a->[0]} @{ $data{$transcript_id}->{'exons'} }) {
            push(@exons, [ @{ $ar_e } ]);

            if (exists $data{$transcript_id}->{'cds'}) {

                my ($start_fivep, $end_fivep);
                my ($start_threep, $end_threep);

                if ($ar_e->[0] < $data{$transcript_id}->{'cleft'}) {
                    
                    $start_threep = $ar_e->[0];
                    
                    if ($ar_e->[1] < $data{$transcript_id}->{'cleft'}) {
                        $end_threep = $ar_e->[1];
                    } else {
                        $end_threep = $data{$transcript_id}->{'cleft'};
                    }
                    
                    push(@threep, [$start_threep, $end_threep]);
                    
                } 

                if ( ( $data{$transcript_id}->{'cright'} > $ar_e->[0] ) && 
                     ( $data{$transcript_id}->{'cright'} < $ar_e->[1] ) ) {
                    
                     $start_fivep = $data{$transcript_id}->{'cright'};
                     $end_fivep = $ar_e->[1];

                     push(@fivep, [$start_fivep, $end_fivep]);

                } elsif ( $ar_e->[0] > $data{$transcript_id}->{'cright'} ) {
                     
                     $start_fivep = $ar_e->[0];
                     $end_fivep = $ar_e->[1];
                    
                     push(@fivep, [$start_fivep, $end_fivep]);

                }
            }
        }
        
        if (exists $data{$transcript_id}->{'cds'}) {
            foreach my $ar_c (sort {$b->[1] <=> $a->[0]} @{ $data{$transcript_id}->{'cds'} }) {
                push(@cds, [ @{ $ar_c } ]);
            }
        }

    } else {
        $LOGGER->logdie("Error getting TSS for $transcript_id strand = $data{$transcript_id}->{'strand'}");
    }
    
    print TSS join("\t", $data{$transcript_id}->{'chr'}, $tss, $tss, join('/', 'tss', $gene_id, $transcript_id), 0, $data{$transcript_id}->{'strand'}),"\n";
    print PROMOTER join("\t", $data{$transcript_id}->{'chr'}, $promoter_start, $promoter_end, join('/', 'promoter', $gene_id, $transcript_id), 0, $data{$transcript_id}->{'strand'}),"\n";

    #print join("\t", $transcript_id, @{$data{$transcript_id}}{'chr', 'gene_id', 'strand', 'eleft', 'eright'}, $tss, $promoter_start, $promoter_end),"\n";
    my $ec=1;
    my @introns;
    foreach my $ar_e (@exons) {
        my ($exon_start, $exon_end) = @{$ar_e};
        
        if ( $#exons > 0 ) {
            if ($data{$transcript_id}->{'strand'} eq '+') {
                #print ">>>$exon_start - $exon_end \t $ec ($#exons)\n";
                # se não é o último exon
                if ( ($ec -1) < $#exons) {
                    #print "\tSET intron (".($ec).") START: $exon_end\n";
                    $introns[$ec-1]->[0] = $exon_end;
                }
                # se não é o primeiro exon
                if ( $ec > 1 ) {
                    #print "\tSET intron (".($ec-1).") END: $exon_start\n";
                    $introns[$ec-2]->[1] = $exon_start;
                }                    
            } else {
                #print ">>>$exon_start - $exon_end \t $ec ($#exons)\n";
                # se não é o último exon
                if ( ($ec -1) < $#exons) {
                    #print "\tSET intron (".($ec).") END: $exon_start\n";
                    $introns[$ec-1]->[1] = $exon_start;
                }
                # se não é o primeiro exon
                if ( $ec > 1 ) {
                    #print "\tSET intron (".($ec-1).") START: $exon_end\n";
                    $introns[$ec-2]->[0] = $exon_end;
                }                    
            }
        }

#            print "\t",join("\t", 'exon_'.$ec++, $exon_start, $exon_end),"\n";
            unless ($data{$transcript_id}->{'cds'}) {            
                print EXON join("\t", $data{$transcript_id}->{'chr'}, $exon_start, $exon_end, join('/', 'exon_'.$ec++, $gene_id, $transcript_id), 0, $data{$transcript_id}->{'strand'}),"\n";
            } else {
                $ec++;
            }
    }
    
    if ($data{$transcript_id}->{'cds'}) {
        my $cc=1;
        my $fc=1;
        my $tc=1;
        foreach my $ar_c (@cds) {
            my ($c_start, $c_end) = @{$ar_c};
            print CDS join("\t", $data{$transcript_id}->{'chr'}, $c_start, $c_end, join('/', 'cds_'.$cc++, $gene_id, $transcript_id), 0, $data{$transcript_id}->{'strand'}),"\n";
        }
        foreach my $ar_f (@fivep) {
            my ($f_start, $f_end) = @{$ar_f};
            print FIVE join("\t", $data{$transcript_id}->{'chr'}, $f_start, $f_end, join('/', '5utr_'.$fc++, $gene_id, $transcript_id), 0, $data{$transcript_id}->{'strand'}),"\n";
        }
        foreach my $ar_t (@threep) {
            my ($t_start, $t_end) = @{$ar_t};
            print THREE join("\t", $data{$transcript_id}->{'chr'}, $t_start, $t_end, join('/', '3utr_'.$tc++, $gene_id, $transcript_id), 0, $data{$transcript_id}->{'strand'}),"\n";
        }
    }

    if ($#introns>0) {
        my $ic = 1;
        foreach my $ar_i (@introns) {
            my ($intron_start, $intron_end) = @{$ar_i};
            
            $LOGGER->logdie('Problem ['.$intron_start.', '.$intron_end.') '.$transcript_id) unless ((defined $intron_start)&&(defined $intron_end));
#            print "\t",join("\t", 'intron_'.$ic++, $intron_start, $intron_end),"\n";
            print INTRON join("\t", $data{$transcript_id}->{'chr'}, $intron_start, $intron_end, join('/', 'intron_'.$ic++, $gene_id, $transcript_id), 0, $data{$transcript_id}->{'strand'}),"\n";
        }
    }
}

close(TSS);
close(PROMOTER);
close(INTRON);
close(EXON);
close(CDS);
close(FIVE);
close(THREE);

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
        -s      --chrsizefile   Chromosome size file
        -p      --promoter_size Promoter size in bp [Default: 1000]
        -o      --outbase       Output base [Default ./coords]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub compn {
    my ($x, $y) = @_;
    my ($xn,$yn);
    if ($x=~/^Group/) {
       ($xn)=$x=~/(\d+)$/;
       if ($x=~/^GroupUn/) {
            $xn+=1000;
       }
    } elsif ($x=~/^MT/) {
        $xn=0;
    } else {
        $LOGGER->logdie("Error: $x");
    }

    if ($y=~/Group/) {
       ($yn)=$y=~/(\d+)$/;
       if ($y=~/^GroupUn/) {
            $yn+=1000;
       }
    } elsif ($y=~/^MT/) {
        $yn=0;
    } else {
        $LOGGER->logdie("Error: $y");
    }
    return ($xn <=> $yn);
}

