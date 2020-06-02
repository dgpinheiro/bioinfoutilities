#!/usr/bin/env perl
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
#  Copyright (C) 2020  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho"
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

Copyright (c) 2020 Universidade Estadual Paulista "Júlio de Mesquita Filho"

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Getopt::Long;

use FileHandle;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}


my (	$level, 
	$infile, 
	$nbases,
	$outdir,
    $nthreads,
    %defines
	);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions(	"h|?|help" => sub { &Usage(); },
		"l|level=s"=> \$level,
		"i|infile=s"=>\$infile,
        "t|threads=s"=>\$nthreads,
        "d|define=s" => \%defines,
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

$nthreads||=2;

$LOGGER->logdie("Missing input file (kraken2/Kaiju output)") if (!defined $infile);
$LOGGER->logdie("Wrong input file ($infile)") if (! -e $infile);

$LOGGER->logdie("Missing output directory") if (! defined $outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") if (! -e $outdir);

my %valid;

my %fh;
foreach my $id (keys %defines) {
    $fh{$id} =  FileHandle->new(">$outdir/$defines{$id}");
    if (!defined $fh{$id}) {
        $fh{$id}->close();
    }
                
    my $taxonkit = `echo $id | taxonkit lineage -j $nthreads -t `;
    my $taxids = [ split("\t", $taxonkit) ]->[2];
    foreach my $k ( split(/;/, $taxids) ) {
        $valid{$k}->{$id} = undef;
    }
    my $downtaxa=`taxonkit list --show-rank --show-name --ids $id`;
    my @taxa = split("\n", $downtaxa);
    foreach my $t (@taxa) {
        my ($tid) = $t=~/^\s*(\d+)/;
        #print ">>>$tid\n";
        $valid{$tid}->{$id} = undef;
    }
    
}

open(IN, "<$infile") or $LOGGER->logdie($!);

my %lineage;

while(<IN>) {
    chomp;
    my ($status, $read, $taxinfo) = split(/\t/, $_);
    my ($taxid);
    if ($taxinfo=~/taxid (\d+)/) {
        $taxid=$1;
    } elsif ($taxinfo=~/^\d+$/) {
        $taxid=$taxinfo;
    } else {
        $LOGGER->logdie("TaxInfo column not recognized ($taxinfo)");
    }

    if ($status eq 'C') {
        
        if ($taxid) {
            
            if ($valid{$taxid}) {
                foreach my $k (keys %{ $valid{$taxid} }) {
                    
                    #print { $fh{$k} } join("\t", $status, $read, $taxid),"\n";
                    print { $fh{$k} } $read,"\n";
                    
                }
            }
            
        } else {
            $LOGGER->logdie("Not found taxid ($taxinfo)");
        }
        
    } else {
            foreach my $k (keys %defines) {
                #print { $fh{$k} } $status,"\t",$read,"\n";
                print { $fh{$k} } $read,"\n";
            }
    }        
}

close(IN);


# Subroutines

sub Usage {
    my ($msg) = @_;

    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2020 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

    -h  --help      Help
    -l  --level     Log level [Default: FATAL]
    -i  --infile    Input file
    -o  --outdir    Output directory
    -t  --threads   Number of threads [Default: 2]
    -d  --defines   Define base name of output file for each taxon id. Example -d 2="Bacteria.txt"

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

