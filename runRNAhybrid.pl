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

use Bio::SeqIO;
use File::Temp qw/ tempfile tempdir /;
use Cwd qw(chdir getcwd abs_path);

my ($level, $qfile, $tfile, $pfile, $config, $outdir, $prefix);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "q|queries=s"=>\$qfile,
            "t|targets=s"=>\$tfile,
            "p|parameters=s"=>\$pfile,
            "c|config=s"=>\$config,
            "o|outdir=s"=>\$outdir,
            "r|prefix=s"=>\$prefix
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

$LOGGER->logdie("Missing query file") unless ($qfile);
$LOGGER->logdie("Wrong query file ($qfile)") unless (-e $qfile);

$LOGGER->logdie("Missing target file") unless ($qfile);
$LOGGER->logdie("Wrong target file ($tfile)") unless (-e $tfile);

$LOGGER->logdie("Missing parameter file") unless ($tfile);
$LOGGER->logdie("Wrong parameter file ($pfile)") unless (-e $pfile);

my %parameter;
open(PAR, "<", $pfile) or $LOGGER->logdie($!);
while(<PAR>){
    chomp;
    next if ($_=~/^target too long/);
    my ($mir, $ndatapoints, $locpar, $scalepar) = split(/\s/, $_);
    $parameter{$mir} =  {
                            'xi'=>$locpar,      # location parameter
                            'theta'=>$scalepar  # scale parameter
                        };
}
close(PAR);

$config||='-f 2,7 -c -e -20 -p 0.005 -g png -m 100000';

$prefix||='targets_';

if ($outdir) {
    $LOGGER->logdie("Wrong output directory ($outdir)") unless (-d $outdir);
}
else {
    $outdir = '.';
}

my $queriesin = Bio::SeqIO->new(-file=>$qfile, -format=>'FASTA');

{
    my $targetsin = Bio::SeqIO->new(-file=>$tfile, -format=>'FASTA');
}

my $tmpdir = tempdir( template=> 'rRhXXXX', CLEANUP => 1 );
print $tmpdir,"\n";

my $cwd = getcwd();
my $abs_path_tfile = abs_path($tfile);
while ( my $seq = $queriesin->next_seq() ) {
    my $mir = $seq->display_id();

    $LOGGER->logdie("Missing parameters for $mir") unless (exists $parameter{$mir});

    my $qfile = "$tmpdir/$mir.fa";
    {
        my $qout = Bio::SeqIO->new(-file=>">$qfile", -format=>'FASTA');
        $qout->write_seq($seq);
    }
    
    my $outfile = "$outdir/$prefix".$mir.'.txt';
    {
        open(OUT, ">", $outfile) or $LOGGER->logdie($!);
        close(OUT);
    }
    my $abs_path_outfile = abs_path($outfile);
    my $abs_path_qfile = abs_path($qfile);

    chdir($outdir);
    # RUN RNAhybrid
    `RNAhybrid -d $parameter{$mir}->{'xi'},$parameter{$mir}->{'theta'} $config -q $abs_path_qfile -t $abs_path_tfile > $abs_path_outfile`;
    chdir($cwd);
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
        -q      --queries       Queries (multi-fasta file with mature miRNAs)
        -t      --targets       Targets (multi-fasta file with target mRNAs)
        -o      --outdir        Output directory (default: .)
        -r      --prefix        Prefix of output files (default: targets_)
        -p      --parameters    Calibrated parameters (output from RNA calibrate for each miRNA)
        -c      --config        Configure RNAhybrid command line arguments except -d, -q and -t (default: '-f 2,7 -c -e -20 -p 0.005 -g png -m 100000'

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

