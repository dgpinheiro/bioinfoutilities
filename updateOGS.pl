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

my ($level, $user, $password);

GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "u|user=s"=>\$user,
            "p|password=s"=>\$password
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

$LOGGER->logdie("Missing BeeBase user.") unless ($user);
$LOGGER->logdie("Missing BeeBase password.") unless ($password);

`wget --save-headers --timestamping --http-user=$user --http-password=$password http://www.hgsc.bcm.tmc.edu/collaborations/insects/bees/`;

my $release;
open(IN, "<", './index.html') or $LOGGER->logdie($!);
while(<IN>) {
    chomp;
    if ($_=~/Amel_4\.5_OGSv([0-9\.]+)/) {
        my $r=$1;
        if ($release) {
            if ($r > $release) {
                $release = $r;
            }
        }
        else {
            $release = $r;
        }
    }
        
}
close(IN);

$LOGGER->logdie("Not found release number") unless ($release);

print "OGS Release: $release\n";

my $outdir = "./OGSv$release";
mkdir($outdir) unless (-e $outdir);

my $cmd;

$cmd = "wget --directory-prefix=$outdir --timestamping --http-user=$user --http-password=$password http://www.hgsc.bcm.tmc.edu/collaborations/insects/bees/Amel_4.5_OGSv$release/amel_OGSv$release.gff3.gz";
`$cmd`;

$cmd = "gunzip $outdir/amel_OGSv$release.gff3.gz";
`$cmd`;



my @file = (
'CHR_LG1/ame_ref_Amel_4.5_chrLG1.fa',
'CHR_LG2/ame_ref_Amel_4.5_chrLG2.fa',
'CHR_LG3/ame_ref_Amel_4.5_chrLG3.fa',
'CHR_LG4/ame_ref_Amel_4.5_chrLG4.fa',
'CHR_LG5/ame_ref_Amel_4.5_chrLG5.fa',
'CHR_LG6/ame_ref_Amel_4.5_chrLG6.fa',  
'CHR_LG7/ame_ref_Amel_4.5_chrLG7.fa', 
'CHR_LG8/ame_ref_Amel_4.5_chrLG8.fa', 
'CHR_LG9/ame_ref_Amel_4.5_chrLG9.fa', 
'CHR_LG10/ame_ref_Amel_4.5_chrLG10.fa',
'CHR_LG11/ame_ref_Amel_4.5_chrLG11.fa',
'CHR_LG12/ame_ref_Amel_4.5_chrLG12.fa',
'CHR_LG13/ame_ref_Amel_4.5_chrLG13.fa',
'CHR_LG14/ame_ref_Amel_4.5_chrLG14.fa',
'CHR_LG15/ame_ref_Amel_4.5_chrLG15.fa',
'CHR_LG16/ame_ref_Amel_4.5_chrLG16.fa',
'CHR_MT/ame_ref_Amel_4.5_chrMT.fa',
'CHR_Un/ame_ref_Amel_4.5_chrUn.fa'
);

my $indir = "./Amel_4.5";
mkdir($indir) unless (-e $indir);
mkdir("$indir/fa") unless (-e "$indir/fa");
mkdir("$indir/tmp") unless (-e "$indir/tmp");

foreach my $df (@file) {
    $cmd = "wget --directory-prefix=$indir/ --timestamping ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/$df.gz";
    `$cmd`;
    my $f = $df;
    $f=~s/^CHR_[^\/]+\///;
    $cmd = "gunzip $indir/$f.gz -c > $indir/tmp/$f";
    `$cmd`;
}

$cmd = "./prepareNCBIchrfiles.pl -i Amel_4.5/tmp/ -o Amel_4.5/fa/";
`$cmd`;

$cmd="cat Amel_4.5/fa/*.fa > Amel_4.5/Amel_4.5.fa";
`$cmd`;

mkdir("$outdir/RNAs") unless (-e "$outdir/RNAs");

$cmd = 'perl ../extractFromGenomeByGFF.pl -d '.$indir.'/fa -f '.$outdir.'/amel_OGSv'.$release.'.gff3 -o '.$outdir.'/RNAs -p OGS-all-'.$release.'_ -s '.$outdir.'/OGS-all-'.$release.'.txt';
`$cmd`;

mkdir("$outdir/RNAs1000") unless (-e "$outdir/RNAs1000");

$cmd = 'perl ../extractFromGenomeByGFF.pl -m 1000 -d '.$indir.'/fa -f '.$outdir.'/amel_OGSv'.$release.'.gff3 -o '.$outdir.'/RNAs1000 -p OGS-all-1000-'.$release.'_ -s '.$outdir.'/OGS-all-1000-'.$release.'.txt';
`$cmd`;

mkdir("$outdir/RNAs2000") unless (-e "$outdir/RNAs2000");

$cmd = 'perl ../extractFromGenomeByGFF.pl -m 2000 -d '.$indir.'/fa -f '.$outdir.'/amel_OGSv'.$release.'.gff3 -o '.$outdir.'/RNAs2000 -p OGS-all-2000-'.$release.'_ -s '.$outdir.'/OGS-all-2000-'.$release.'.txt';
#`$cmd`;

my $bgoutdir = "./OGSv$release".'byGroup';
mkdir($bgoutdir) unless (-e $bgoutdir);


$cmd = 'perl scaff2superscaff.pl '.$outdir.'/amel_OGSv'.$release.'.gff3  > '.$bgoutdir.'/amel_OGSv'.$release.'.byGroup.gff3';
#`$cmd`;

mkdir("$bgoutdir/RNAs") unless (-e "$bgoutdir/RNAs");

my $ncbidir ='/srv/projects/bioinfoutilities/NCBI';

$cmd = 'perl ../extractFromGenomeByGFF.pl -d '.$ncbidir.'/apiMel4_5 -f '.$bgoutdir.'/amel_OGSv'.$release.'.byGroup.gff3 -o '.$bgoutdir.'/RNAs -p OGS-all-2000-'.$release.'_ -s '.$bgoutdir.'/OGS-all-'.$release.'.txt';
`$cmd`;

mkdir("$bgoutdir/RNAs2000") unless (-e "$bgoutdir/RNAs2000");

$cmd = 'perl ../extractFromGenomeByGFF.pl -m 2000 -d '.$ncbidir.'/apiMel4_5 -f '.$bgoutdir.'/amel_OGSv'.$release.'.byGroup.gff3 -o '.$bgoutdir.'/RNAs2000 -p OGS-all-2000-'.$release.'_ -s '.$bgoutdir.'/OGS-all-2000-'.$release.'.txt';
`$cmd`;


# Adjust gff
$cmd = 'echo "##gff-version 3" > '.$bgoutdir.'/OGSv'.$release.'.byGroup_adj.gff3';
`$cmd`;
$cmd = 'cat '.$bgoutdir.'/amel_OGSv'.$release.'.byGroup.gff3 >> '.$bgoutdir.'/OGSv'.$release.'.byGroup_adj.gff3';
`$cmd`;

$cmd = 'gff3ToGenePred '.$bgoutdir.'/OGSv'.$release.'.byGroup_adj.gff3 '.$bgoutdir.'/OGSv'.$release.'.byGroup_adj.Gp';
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

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -u      --user      BeeBase user
        -p      --password  BeeBase password

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

