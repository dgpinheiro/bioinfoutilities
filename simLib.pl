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
#  Copyright (C) 2018  Universidade Estadual Paulista - UNESP
#
#  Universidade Estadual Paulista - UNESP
#  Laboratório de Bioinformática
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
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

Copyright (c) 2018 Universidade Estadual Paulista - UNESP

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

my ($level, $infile, $outdir, $abdfile, $prefix, $number, $insert_length, $fivepadapt, $threepadapt, $read_length, $runfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "o|outdir=s"=>\$outdir,
            "a|abdfile=s"=>\$abdfile,
            "p|prefix=s"=>\$prefix,
            "n|number=i"=>\$number,
            "s|inslen=i"=>\$insert_length,
            "fpa|fivepadapt=s"=>\$fivepadapt,
            "tpa|threepadapt=s"=>\$threepadapt,
            "rl|readlength=i"=>\$read_length,
            "rf|runfile=s"=>\$runfile
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

$number||=1000;
$prefix||='simLib';
$insert_length||=400;
$fivepadapt||='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG';
$threepadapt||='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT';
$read_length||=151;

$LOGGER->logdie("Missing simNGS runfile") unless ($runfile);
$LOGGER->logdie("Wrong simNGS runfile ($!)") unless (-e $runfile);

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($!)") unless (-e $infile);

my %data;
my @o;
open(IN, "<", $infile) or $LOGGER->logdie($!);
while(<IN>) {
    chomp;
    if ($_=~/^>(\S+)/) {
        $data{$1}=undef;
        push(@o, $1);
    }
}
close(IN);

use Bio::SeqIO;

my $seqin=Bio::SeqIO->new(-file=>$infile, -format=>'FASTA');
my %seq;
while (my $seq=$seqin->next_seq()) {
    $seq->description("");
    $seq{$seq->display_id()}=$seq;
}

$LOGGER->logdie("Missing abundance file") unless ($abdfile);
$LOGGER->logdie("Wrong abundance file ($!)") unless (-e $abdfile);

# see https://docstore.mik.ua/orelly/perl4/cook/ch02_10.htm
use Math::Random qw(random_normal);

$outdir||='.';

my $total=0;
open(ABD, "<", $abdfile) or $LOGGER->logdie($!);
while(<ABD>) {
    next if ($_=~/^\s*$/);
    chomp;
    my ($id, $m) = split(/\t/, $_);
    if ($id) {
        if (exists $data{$id}) {
            if ($data{$id}) {
                $LOGGER->logdie("Error in abundance file: abundance for sequence ($id) already specified: $data{$id}");
            }
            $total+=$m+0;
            $data{$id} = random_normal(1, $m+0, 0.01);
            $data{$id} = 0 if ($data{$id} < 0);
        }
    } else {
        $LOGGER->logdie("Wrong line in abundance file ($abdfile)! Not found id in ($_)");
    }
}
close(ABD);

my $rtotal=sprintf("%.1f", $total);
$LOGGER->logdie("Sum of proportions in abundance file ($abdfile) distinct than 1.0 ($rtotal)") if ($rtotal != 1.0);

use File::Temp qw/ tempfile tempdir /;


my $dir = tempdir( CLEANUP => 1 );
my ($fh_abd, $filename_abd) = tempfile( 'simLibABDXXXXX', DIR => $dir);
my ($fh_seq, $filename_seq) = tempfile( 'simLibFAXXXXX', DIR => $dir);

my $seqout=Bio::SeqIO->new(-fh=>\*$fh_seq, -format=>'FASTA');

foreach my $id (@o) {
    if (!defined $data{$id}) {
        $LOGGER->logdie("Sequence ($id) is not in abundance file ($abdfile)");
    }
    print {$fh_abd} $data{$id},"\n";
    $seqout->write_seq($seq{$id});
}

my $output_prefix="$outdir/$prefix";

my $cmd='cat '.$filename_seq.' | sed \'s/^>\(\S\+\).*/>\1/\' | simLibrary -s opposite -m '.$filename_abd.' -n '.$number.' -i '.$insert_length.' 2> '.$output_prefix.'.simLibrary.log.err | simNGS -a '.$fivepadapt.':'.$threepadapt.' -p paired '.$runfile.' -n '.$read_length.' 2> '.$output_prefix.'.simNGS.log.err | perl -slane \' if ($.%4==1) { if ($_=~/^@(\S+)\s+(\S+)\s+\(Strand ([+-]) Offset (\d+)--(\d+)\)/) { $_="\@$name:$1:$2:".(($3 eq "+") ? 1 : 0).":$4:$5"; }  else { die; } } print $_;\' -- -name="'.$prefix.'" | paste - - - - - - - - | bash -c \'tee >(cut -f 1-4 | tr "\t" "\n" > '.$output_prefix.'_R1.fastq'.') \' | cut -f 5-8 | tr "\t" "\n" > '.$output_prefix.'_R2.fastq';

#print "$cmd\n";
system($cmd);

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
    -o      --outdir        Output directory [Default: .]
    -a      --abdfile       Abundance file
    -p      --prefix        Prefix name
    -n      --number        Number of fragments [Default: 1000]
    -s      --inslen        Mean length of insert [Default: 400]
    -fpa    --fivepadapt    5' Adapter sequence [Default: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG']
    -tpa    --threepadapt   3' Adapter sequence [Default: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT']
    -rl     --readlen       Read length [Default: 151]
    -rf     --runfile       Run file 

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

