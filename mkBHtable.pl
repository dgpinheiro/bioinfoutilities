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

use Bio::SeqIO;
use File::Spec;
use File::Basename;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level);

my ($rnatxtdir, $peptxtdir, $outdir, $mappeprnafile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "m|mappeprna=s"=>\$mappeprnafile,
            "r|rnatxtdir=s"=>\$rnatxtdir,
            "p|peptxtdir=s"=>\$peptxtdir,
            "o|outdir=s"=>\$outdir
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

$LOGGER->logdie("Missing pep2rna file") unless ($mappeprnafile);
$LOGGER->logdie("Wrong pep2rna file ($mappeprnafile)") unless (-e $mappeprnafile);

$LOGGER->logdie("Missing rna txt input directory") unless ($rnatxtdir);
$LOGGER->logdie("Wrong rna txt input directory ($rnatxtdir)") unless (-e $rnatxtdir);

$LOGGER->logdie("Missing pep txt input directory") unless ($peptxtdir);
$LOGGER->logdie("Wrong pep txt input directory ($peptxtdir)") unless (-e $peptxtdir);

if ($outdir) {
    $outdir=~s/\/+$//;
    $LOGGER->logdie("Wrong output directory ($outdir)") unless (-e $outdir);
}
else {
    $outdir='.';
}


$LOGGER->info("LOADING pep2rna...");

my %pep2rna;
my %rna2pep;
open(IN, "<", $mappeprnafile) or $LOGGER->logdie($!);

while(<IN>) {
    chomp;
    my ($db, $pep, $rna) = split(/\t/, $_);
    $pep2rna{$db}->{$pep} = $rna;
    $rna2pep{$db}->{$rna} = $pep;
}

close(IN);

$LOGGER->info("LOADING mappings:");


my %rnafromto;
$LOGGER->info("\tRNA ...");
&load_fromto($rnatxtdir, \%rnafromto, 1);

my @rnadb=sort keys %rnafromto;

my %pepfromto;
$LOGGER->info("\tPeptide ...");
&load_fromto($peptxtdir, \%pepfromto, 2);

my @pepdb=sort keys %pepfromto;

foreach my $qdb (@rnadb) {
    foreach my $sdb (@rnadb) {
        next if ($qdb eq $sdb);
        open(OUT,">", "$outdir/MAP_".$qdb.'_x_'.$sdb.'.txt') or $LOGGER->logdie($!);

        foreach my $qrna (keys %{ $rnafromto{$qdb}->{$sdb} }) {
            my $qpep=$rna2pep{$qdb}->{$qrna};

            print OUT $qrna,"\t",(($qpep) ? $qpep : '');

#            foreach my $srna (@rnadb) {
#                next if ($qdb eq $sdb);
                
                print OUT "\t", join(";", map { $_.'['.$rnafromto{$qdb}->{$sdb}->{$qrna}->{$_}.']' } keys %{ $rnafromto{$qdb}->{$sdb}->{$qrna} }),
                          "\t", ( ($qpep) ?  join(";", map { $_.'['.$pepfromto{$qdb}->{$sdb}->{$qpep}->{$_}.']' } keys %{ $pepfromto{$qdb}->{$sdb}->{$qpep} }) : '');
#            }
            print OUT "\n";
        }
        
        close(OUT);
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

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -r      --rnatxtdir     rna/txt input directory
        -p      --peptxtdir     pep/txt input directory
        -m      --mappeprnafile Mapping of pep to rna
        -o      --outdir        Output directory [Default: .]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub load_fromto {
    my ($dir,$hr_fromto,$type) = @_;
    foreach my $cmpref (glob("$dir/*_out.txt")) {
        my $bn=basename($cmpref);
        next if ($bn=~/^EXACT/);
        if ($bn=~/^(\S+)_x_(\S+)_out\.txt$/) {
            if ($1 ne $2) {
                my $qdb=$1;
                my $sdb=$2;
                #print STDERR "\t",$qdb," x ",$sdb,"\n";
                open(IN, "<", $cmpref) or $LOGGER->logdie($!);
                while(<IN>) {
                    chomp;
                    #query subject relationship
                    my ($qstr, $sstr, $r) = split(/\t/, $_);
                    my ($q) = &get_id($qstr, $type);
                    foreach my $esstr (split(';', $sstr)) {
                        my ($s) = &get_id($esstr, $type);
                        #print STDERR $qdb,"\t",$sdb,"\t",$q,"\t",$s,"\n";
                        $hr_fromto->{$qdb}->{$sdb}->{$q}->{$s} = $r;
                    }                        
                }
                close(IN);
            }
        }
    }
}

sub get_id {
    my ($rawid,$type) = @_;
    
    my $id;
    if ( $rawid =~ /(GB\d+\-[RP][A-Z](?:_\d+)?)/ ) {
        $id = $1;
        if ($type==1) {
            $id =~ s/-P/-R/;
        }
        elsif ($type==2) {
            $id =~ s/-R/-P/;
        }
        else {
            $LOGGER->logdie("Wrong type for get_id() function");
        } 
    }
    elsif ( $rawid =~ /([XN][MPR]_\d+)/ ) {
        if ($type==1) {
            ($id) = $_=~/([XN][MR]_\d+)/;
        }
        elsif ($type==2) {
            ($id) = $_=~/([XN]P_\d+)/;
        }
        else {
            $LOGGER->logdie("Wrong type for get_id() function");
        }
    }
    elsif ($rawid =~ /(LOC\d+#[^\|]+)/) {
        $id=$1;
    }
    elsif ( $rawid=~ /((?:ND\d+(?:\w\d*)?|CYTB|COX\d+|ATP\d+)#[^\|]+)/) {
        $id=$1;
    }
    else {
        $LOGGER->logdie("Pattern not considered: $rawid");
    }
    $LOGGER->logdie("Not found id ($rawid) type $type") unless ($id);
    return($id);
}

