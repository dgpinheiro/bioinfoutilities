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

my ($level);

GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level
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


use Net::FTP;

my $ftp = Net::FTP->new( "ftp.ncbi.nih.gov", Debug => 0 )
  or die "Cannot connect to some.host.name: $@";

$ftp->login( "anonymous", '-anonymous@' )
  or die "Cannot login ", $ftp->message;

$ftp->cwd("/refseq/release/release-catalog/")
  or die "Cannot change working directory ", $ftp->message;

my @files = $ftp->ls()
  or die "get failed ", $ftp->message;

my $release;
foreach my $f (@files) {
    if ($f =~ /release(\d+)/) {
        $release = $1;
        last;
    }
}

$ftp->cwd("/genomes/Apis_mellifera/RNA/")
  or die "Cannot change working directory ", $ftp->message;

$ftp->binary();

$ftp->get("rna.gbk.gz")
      or die "get failed ", $ftp->message;

$ftp->quit;

$LOGGER->logdie("Not found release number") unless ($release);


print "RefSeq Release: $release\n";

my $cmd;

my $outdir = "./refseq.$release";
mkdir($outdir) unless (-e $outdir);

if (-e "./rna.gbk.gz") {
    `gunzip rna.gbk.gz`;
    `mv ./rna.gbk $outdir/`;
}


$cmd = "perl ../getNCBIseqs.pl -d nuccore -t 'PRJNA13343[BioProject]+AND+srcdb+refseq[Properties]+NOT+\"genomic+scaffold\"[Title]' -o $outdir/AmRefSeqDNA-$release.gbk -e gbwithparts -r 1";
`$cmd`;

$cmd = "perl ../getNCBIseqs.pl -d nuccore -t 'PRJNA13343[BioProject]+AND+srcdb+refseq[Properties]+AND+unplaced[Title]' -o $outdir/AmRefSeqDNA-unplaced-$release.gbk -e gbwithparts -r 500";
`$cmd`;

$cmd = "cat $outdir/AmRefSeqDNA-$release.gbk $outdir/AmRefSeqDNA-unplaced-$release.gbk > $outdir/AmRefSeqDNA-all-$release.gbk";
`$cmd`;

mkdir("$outdir/RNAs") unless (-e "$outdir/RNAs");

$cmd = "perl ../fullparseToGff.pl -i $outdir/AmRefSeqDNA-all-$release.gbk -o ./ORDER.txt -r $outdir/rna.gbk > $outdir/AmRefSeqDNA-all-$release.gff 2> $outdir/AmRefSeqDNA-all-$release.gff.stderr";
`$cmd`;

$cmd = 'gff3ToGenePred '.$outdir.'/AmRefSeqDNA-all-'.$release.'.gff '.$outdir.'/AmRefSeqDNA-all-'.$release.'.Gp';
`$cmd`;

$cmd = 'perl ../extractFromGenomeByGFF.pl -g apiMel4_5.fa -d ./apiMel4_5/ -f '.$outdir.'/AmRefSeqDNA-all-'.$release.'.gff -o ./refseq.'.$release.'/RNAs/ -p AmRefSeqRNA-all-'.$release.'_ -s refseq.'.$release.'/AmRefSeqRNA-all-'.$release.'.txt';
`$cmd`;

mkdir("$outdir/RNAs2000") unless (-e "$outdir/RNAs2000");

$cmd = 'perl ../extractFromGenomeByGFF.pl -m 2000 -g apiMel4_5.fa -d ./apiMel4_5/ -f '.$outdir.'/AmRefSeqDNA-all-'.$release.'.gff -o ./refseq.'.$release.'/RNAs2000/ -p AmRefSeqRNA-all-2000-'.$release.'_ -s refseq.'.$release.'/AmRefSeqRNA-all-2000-'.$release.'.txt';
`$cmd`;

mkdir("$outdir/RNAs1000") unless (-e "$outdir/RNAs1000");

$cmd = 'perl ../extractFromGenomeByGFF.pl -m 1000 -g apiMel4_5.fa -d ./apiMel4_5/ -f '.$outdir.'/AmRefSeqDNA-all-'.$release.'.gff -o ./refseq.'.$release.'/RNAs1000/ -p AmRefSeqRNA-all-1000-'.$release.'_ -s refseq.'.$release.'/AmRefSeqRNA-all-1000-'.$release.'.txt';
`$cmd`;

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help  Help
        -l      --level Log level [Default: FATAL]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

