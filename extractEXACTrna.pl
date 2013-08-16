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

my ($level, $input1, $input2, $beddir);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i1|input1=s"=>\$input1,
            "i2|input2=s"=>\$input2,
            "b|beddir=s"=>\$beddir
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

use File::Basename;

$LOGGER->logdie("Missing input1 file") unless ($input1);
$LOGGER->logdie("Wrong input1 file ($input1)") unless (-e $input1);

$LOGGER->logdie("Missing input2 file") unless ($input2);
$LOGGER->logdie("Wrong input2 file ($input2)") unless (-e $input2);

$LOGGER->logdie("Missing beddir") unless ($beddir);

my $b1 = basename($input1, '.fa');
my $b2 = basename($input2, '.fa');

my $bed1f = $beddir.'/'.$b1.'.bed';
my $bed2f = $beddir.'/'.$b2.'.bed';

$LOGGER->logdie("Not found input1 bed file ($bed1f)") unless (-e $bed1f);
$LOGGER->logdie("Not found input2 bed file ($bed2f)") unless (-e $bed2f);


use Bio::SeqIO;

my %data1;
{
    my $in1 = Bio::SeqIO->new(-file=>$input1, -format=>'FASTA');
    while(my $seq = $in1->next_seq() ) {
        $data1{ $seq->display_id() } = $seq->seq();
        $data1{ $seq->display_id() }=~s/\*//g;
    }
}

my %data2;
{
    my $in2 = Bio::SeqIO->new(-file=>$input2, -format=>'FASTA');
    while(my $seq = $in2->next_seq() ) {
        $data2{ $seq->display_id() } = $seq->seq();
        $data2{ $seq->display_id() }=~s/\*//g;
    }
}


my %rel;
foreach my $s1 (keys %data1) {
    foreach my $s2 (keys %data2) {
        if ($data1{$s1} eq $data2{$s2}) {
#            print $s1,"\t",$s2,"\n";
            $rel{$s1}->{$s2} = undef;
        }
    }
}

use File::Temp qw/tempdir/;

foreach my $s1 (keys %rel) {
    my $tmpdir = tempdir( CLEANUP => 1 );
    my @r = keys %{ $rel{$s1} };
    if (scalar(@r) > 1) {
        my ($b1) = $s1=~/\|?([^\|]+)$/;
        $LOGGER->logdie("Not found ACC for $s1") unless ($b1);
        $b1=~s/_Group.*//;
        $b1=~s/P/R/;
        
        `grep '$b1' $bed1f | cut -f 1,2,3,4,5,6 > $tmpdir/BED1.bed`;
        if (-z "$tmpdir/BED1.bed") {
            $LOGGER->logdie("Difficult to map $s1 to any of ".join(';', @r).' because '.$b1.' not found in '.$bed1f);
            next;
        }

        my @x;
        foreach my $p (@r) {
            my ($b2) = $p=~/\|?([^\|]+)$/;

            $LOGGER->logdie("Not found ACC for $p") unless ($b2);
            $b2=~s/_Group.*//;
            $b2=~s/P/R/;

            `grep '$b2' $bed2f | cut -f 1,2,3,4,5,6 > $tmpdir/BED2.bed`;
            if (-z "$tmpdir/BED2.bed") {
                $LOGGER->logdie("Difficult to map $s1 to any of ".join(';', @r).' because '.$b2.' not found in '.$bed2f);
                next;
            }
            `intersectBed -s -b $tmpdir/BED1.bed -a $tmpdir/BED2.bed -f 0.9 -r -wo > $tmpdir/INTER.txt`;
            open(INTER, "<", "$tmpdir/INTER.txt") or $LOGGER->logdie($!);
            while(<INTER>){
                chomp;
                my (@d) = split(/\t/, $_);
                $LOGGER->logdie("Not found overlap length") unless (defined $d[12]);
                push(@x, [ $p , $d[12] ]);
            }
            close(INTER);
        }
        
        if (scalar(@x) == 1) {
            print $s1,"\t",$x[0]->[0],"\n";
        }
        elsif (scalar(@x) > 1) {
            my @sel;
            my $best=0;
            # selecionar o de maior (ou igual) tamanho de intersecção
            foreach my $y (sort {$b->[1] <=> $a->[1]} @x) {
                if ($best <= $y->[1]) {
                    push(@sel, $y);
                    $best = $y->[1];
                }
                else {
                    last;
                }         
            }
            if (scalar(@sel)==1) {
                print $s1,"\t",$sel[0]->[0],"\n";
            }
            else {
                my %selhash;
                @selhash{ (map { $_->[0] } @sel) } = ( map { $_->[1] } @sel );

                $LOGGER->logwarn("Difficult to map $s1 to any of ".join(';', @r).' because '.join(';', map { $_."[".$selhash{$_}."]"  } keys %selhash).' have the same intersection length.');
                
            }
        }
        else {
                my %selhash;
                @selhash{@r} = @r;

                $LOGGER->logwarn("Not found intersection among $s1 and any of ".join(';', @r).'.');
        }
     }
    else {
        print $s1,"\t", $r[0],"\n";
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
        -i1     --input1    Input file 1
        -i2     --input2    Input file 2
        -b      --beddir    BED dir

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

