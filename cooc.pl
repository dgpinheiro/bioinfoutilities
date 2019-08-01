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

use constant JAR_PATH => "/usr/local/bioinfo/CoNet3/lib/CoNet.jar";

use File::Spec;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level,$infile,$featfile,$outdir,$basename,$mem,%param,$setfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "f|featfile=s"=>\$featfile,
            "o|outdir=s"=>\$outdir,
            "b|basename=s"=>\$basename,
            "m|mem=i"=>\$mem,
            "s|setfile=s"=>\$setfile
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

$basename||='result';

$LOGGER->logdie("Missing input file - OTU table") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

$LOGGER->logdie("Missing feature file") unless ($featfile);
$LOGGER->logdie("Wrong feature file ($featfile)") unless (-e $featfile);

$LOGGER->logdie("Missing output directory") unless ($outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") unless (-d $outdir);

$LOGGER->logdie("Missing setting file") unless ($setfile);
$LOGGER->logdie("Wrong setting file ($setfile)") unless (-e $setfile);

$infile = File::Spec->rel2abs( $infile );
$featfile = File::Spec->rel2abs( $featfile );
$outdir = File::Spec->rel2abs( $outdir );

open(SET,"<",$setfile) or $LOGGER->logdie($!);
while(<SET>) {
    chomp;
    next if (($_=~/^\s*#/)||($_=~/^\s*$/));
    $LOGGER->logdie("Wrong setting file in line $. ($_)") unless ($_=~/^\s*([^=]+)\s*==\s*([^=]+)\s*$/);
    $param{$1} = $2;
    #print $1,"\t",$2,"\n";
}
close(SET);

$mem||='4000';

my $common_cmd = 'java -Xmx'.$mem.'m -cp '.JAR_PATH.' be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --features '.$featfile.'  --input '.$infile.' ';

my $cooc_cmd='';
foreach my $k ( 'transposefeatures','matchfeaturesamples','method','topbottom','nantreatment','nantreatmentparam','inputformat','matrixtype','stand','filter','ensemblemethods','multigraph','format','networkmergestrategy' ) {
    if (($param{$k} eq 'true')||($param{$k} eq 'false')) {
        $cooc_cmd.=' --'.$k.' ' if ($param{$k} eq 'true');
    } else {
        $cooc_cmd.=' --'.$k.' '.$param{$k}.' ';
    }
}

$LOGGER->logdie("Missing thresholdguessing parameter") unless ($param{'thresholdguessing'});
$LOGGER->logdie("Missing guessingparam parameter") unless ($param{'guessingparam'});

my $threshold_cmd='';
foreach my $k ( 'transposefeatures','matchfeaturesamples','method','topbottom','nantreatment','nantreatmentparam','inputformat','matrixtype','stand','filter','thresholdguessing','guessingparam','ensemblemethods' ) {
    if (($param{$k} eq 'true')||($param{$k} eq 'false')) {
        $threshold_cmd.=' --'.$k.' ' if ($param{$k} eq 'true');
    } else {
        $threshold_cmd.=' --'.$k.' '.$param{$k}.' ';
    }
}


print "Computing thresholds ...\n";

my $conffile=$outdir.'/CoNetConfig.txt';

open(OUT,">",$conffile) or $LOGGER->logdie($!);
print OUT "lineage_separator=".($param{'lineage_separator'}||'');
close(OUT);

my $cmd = $common_cmd.
' '.$threshold_cmd.
' --output '.$outdir.'/'.$basename.'_thresholds.txt'.
' --configfile '.$conffile.
'  > '.$outdir.'/'.$basename.'_thresholdguessing.out.txt '.
' 2> '.$outdir.'/'.$basename.'_thresholdguessing.err.txt';

`$cmd`;
#print $cmd;

print "Doing cooccurrence analysis ...\n";

$cmd = $common_cmd.
' '.$cooc_cmd.
' --ensembleparamfile '.$outdir.'/'.$basename.'_thresholds.txt '.
' --output '.$outdir.'/'.$basename.'_cooccurrenceNetwork.'.$param{'format'}.
'  > '.$outdir.'/'.$basename.'_cooccurrence_analysis.out.txt '.
' 2> '.$outdir.'/'.$basename.'_cooccurrence_analysis.err.txt';

`$cmd`;
#print $cmd,"\n";



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
        -i      --infile    Input file (OTU table)
        -f      --featfile  Feature file
        -o      --outdir    Output directory
        -m      --mem       Memory in MB [Default: 4000]
        -s      --setfile   Settings file
        -b      --basename  Basename [Default: result]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

