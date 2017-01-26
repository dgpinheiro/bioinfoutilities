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

my ($level, $fafile, $rfafile, $cltfile, $outfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "f|fa=s"=>\$fafile,
            "r|rfa=s"=>\$rfafile,
            "c|clt=s"=>\$cltfile,
            "o|outfile=s"=>\$outfile
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

$LOGGER->logdie("Missing consensus fasta file") unless ($fafile);
$LOGGER->logdie("Wrong consensus fasta file ($fafile)") unless (-e $fafile);

$LOGGER->logdie("Missing representative fasta file") unless ($rfafile);
$LOGGER->logdie("Wrong representative fasta file ($rfafile)") unless (-e $rfafile);

$LOGGER->logdie("Missing cluster file") unless ($cltfile);
$LOGGER->logdie("Wrong cluster file ($cltfile)") unless (-e $cltfile);

$LOGGER->logdie("Missing output file") unless ($outfile);

my %cluster;
my $cltid;
open(CLT, "<", $cltfile) or $LOGGER->logdie($!);
while(<CLT>) {
    chomp;
    if ( $_ =~ /^>(.+)/ ) {
        $cltid = $1;
        $cltid =~ s/\s/_/g;
        $cltid =~ s/Cluster/cluster/;
    }
    else {
        if ($cltid) {
            if ( $_ =~ />([^\.\s]+)/ ) {
                my $member=$1;
                $cluster{$member} = $cltid;
            }
        }
    }
}
close(CLT);

my %annot;
open(REF, "<", $rfafile) or $LOGGER->logdie($!);
while(<REF>) {
    chomp;
    if ( $_ =~ /^>(\S+)\s*(.*)/ ) {
        my ($id, $desc) = ($1, $2);
        #print $id,"\t",$desc,"\n";
        if (exists $cluster{ $id }) {
            $annot{$cluster{$id}} = $id.'|'.($desc||'');
        } else {
            $LOGGER->logdie("Not found cluster id for reference fasta ($id)");
        }            
    }
}
close(REF);

open(OUT, ">", $outfile) or $LOGGER->logdie($!);

open(FA, "<", $fafile) or $LOGGER->logdie($!);
while(<FA>) {
    chomp;
    if ($_=~/^>(\S+)/) {
        my ($rawcltid) = $1;
        my ($cltid)=$rawcltid=~/(cluster_\d+)/;
        if ($cltid) {
            #print $cltid,"\n";
            if (exists $annot{$cltid}) {
                print OUT ">$rawcltid\t$annot{$cltid}\n";
            } else {
                $LOGGER->logdie("Not found annotation for cluster ($cltid)");
            }
        } else {
            $LOGGER->logdie("Not found proper cluster identifier ($rawcltid)");
        }
    } else {
        print OUT $_,"\n";
    }
}
close(FA);

close(OUT);

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -f      --fa        CD-HIT consensus fasta file
        -r      --rfa       CD-HIT cluster representative fasta
        -c      --clt       CD-HIT cluster file
        -o      --outfile   Output fasta file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

