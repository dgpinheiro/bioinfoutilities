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
use File::Spec;
use File::chdir;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}


my ($level, $threads, $infile, $outfile, $params, $leavetmp, $tmpdirs);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "t|threads=i"=>\$threads,
            "i|infile=s"=>\$infile,
            "o|outfile=s"=>\$outfile,
            "p|params=s"=>\$params,
            "n|normtmp"=>\$leavetmp,
            "d|tmpdirs=s"=>\$tmpdirs
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

use POSIX ":sys_wait_h";


$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

$LOGGER->logdie("Missing output file") unless ($outfile);

if ($tmpdirs) {
    $LOGGER->logdie("Missing output file") unless (-d $tmpdirs);
} else {
    $tmpdirs='./';
}

$threads||=1;

$params||='';
$params.=' -parameter /usr/local/bioinfo/VersionDistribuableExecutable/DATA/DEFAULT.DAT ';
if ($params !~ /-L\s+\d+/) {
    $params.=' -L 150 ';
}
if ($params !~ /-E/) {
    $params.=' -E ';
}
if ($params !~ /-All/) {
    $params.=' -All ';
}

my $cleantmp;
if ($leavetmp) {
    $cleantmp = 0;
}
else {
    $cleantmp = 1;
}


use Bio::SeqIO;
use File::Temp qw/ tempfile tempdir /;

my $int = Bio::SeqIO->new(-file=>$infile, -format=>'FASTA');

my %SEQ;
my @batch;

{
    my $i = 0;
    while( my $seq = $int->next_seq() ) {
        if ($seq->length() >= 20) {
            $SEQ{ "SEQ$i" } = $seq->display_id();
            $batch[($i % $threads)]->{ "SEQ$i" } = $seq->seq();

            $i++;
        }            
    }
}

my $tmpdir = tempdir( CLEANUP => $cleantmp, DIR=>$tmpdirs );
$LOGGER->info("Working directory: $tmpdir");
    
my $seq_dir = File::Spec->rel2abs("$tmpdir/SEQS");
mkdir("$tmpdir/SEQS/") unless (-d "$tmpdir/SEQS");

for my $i (0..$#batch) {
    foreach my $seqid ( keys %{ $batch[$i] } ) {
        open( SEQ, ">", "$seq_dir/$seqid.fa" ) or $LOGGER->logdie($!);
        print SEQ ">", $seqid, "\n", $batch[$i]->{$seqid};
        close(SEQ);
    }
}

my $auxdir = tempdir( CLEANUP =>$cleantmp, DIR=>$tmpdirs );
$LOGGER->info("Auxiliary temp directory: $auxdir");

my $n = 0;

my $running = 0;
my @pids;
while ($n <= $#batch) {
    $running = 0;
    foreach my $t (glob("$auxdir/*.run")) {
        $running++;
    }
    while ((($running+1) <= $threads)&&($n<=$#batch)) {
        $LOGGER->info("Running $running (launch +1)");

        my $pid = fork();

        if (not defined $pid) {
            $LOGGER->logdie("Resources not available.");
        }
        elsif ($pid == 0) {
            $LOGGER->info("Running CHILD process...");

            # CHILD
            my $touch = File::Spec->rel2abs( $auxdir.'/'.$n.'.run' );
            open(TOUCH,">",$touch) or $LOGGER->logdie("$! ($touch)");
            print TOUCH $pid;
            close(TOUCH);
            sleep(1);

            # BEGIN
            open(STDOUT,">>","$auxdir/$n.log") or $LOGGER->logdie("$! ($auxdir/$n.log)");
            open(STDERR,">>","$auxdir/$n.log") or $LOGGER->logdie("$! ($auxdir/$n.log)");
                
            mkdir("$tmpdir/$n");
            open(ERR, ">", "$tmpdir/$n/ERR.txt") or $LOGGER->logdie($!);
            
            my $batch_dir = File::Spec->rel2abs("$tmpdir/$n");
            my $PWD = $CWD;
            $CWD = $batch_dir;

            mkdir("$batch_dir/ERR");
            mkdir("$batch_dir/OUT");
            
            foreach my $seqid ( keys %{ $batch[$n] } ) {
                    
                unless ( system("miRNAFold -s $seq_dir/$seqid.fa $params > ".$batch_dir.'/OUT/'.$seqid.'.txt') == 0 ) {

                    print ERR "Error with: $seqid\n";

                    open(MSG, ">", "$batch_dir/ERR/".$seqid.'.txt') or $LOGGER->logdie($!);
                    if ( $? == -1 ) {
                        print MSG "failed to execute: $!\n";
                    }
                    elsif ( $? & 127 ) {
                        printf MSG
                          "child died with signal %d, %s coredump\n",
                          ( $? & 127 ),
                          ( $? & 128 ) ? 'with' : 'without';
                    } else {
                         printf MSG "child exited with value %d\n", $? >> 8;
                    }
                    close(MSG);

                } else {
                    open(RES, "<", "$batch_dir/OUT/".$seqid.'.txt') or $LOGGER->logdie($!);
                    $|=1;
                    open(OUT, ">>", "$batch_dir/OUT.txt") or $LOGGER->logdie($!);
                    while(<RES>) {
                        chomp;
                        next if ($_=~/^#|^\./);
                        my (@all) = split(/\t/, $_);
                        $LOGGER->logdie("Not found $seqid") unless (exists $SEQ{$seqid});
                        print OUT join("\t", $SEQ{$seqid}, @all),"\n";
                    }
                    close(RES);
                    close(OUT);
                }
            }
            
            close(ERR);
            
            $CWD = $PWD;

            unlink($auxdir.'/'.$n.'.run');
            sleep(1);
            
            # END
            exit(0);
        } elsif ($pid == -1) {
            die;
        } else {
           $running++;
            $LOGGER->info("Increased running count $running (PID: $pid)");
            $n++;
            push(@pids, $pid); #parent stores children's pids
        }
    }
    sleep(5);
}

for my $pid (@pids) {
    waitpid $pid, 0;
}

$n = 0;
open(OUT, ">", $outfile) or $LOGGER->logdie($!);
while ($n <= $#batch) {
    open(IN, "<", "$tmpdir/$n/OUT.txt") or $LOGGER->logdie($!);
    while(<IN>) {
        print OUT $_;
    }
    close(IN);
    
    open(ERR, "<", "$tmpdir/$n/ERR.txt") or $LOGGER->logdie($!);
    while(<ERR>) {
        chomp;
        print STDERR $_;
    }
    close(ERR);
    print STDERR "\n";
    $n++;
}
close(OUT);

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
        -t      --threads   Threads [Default: 1]
        -i      --infile    Input file
        -o      --outfile   Output file
        -p      --params    Other miRNAFold parameters [Default: -L 150 -E -All]
        -n      --normtmp   Not remove tmp dir
        -d      --tmpdirs   Temporary directories [Default: ./]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

