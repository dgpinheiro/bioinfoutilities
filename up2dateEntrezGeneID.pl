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

my ($level, $infile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile
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

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

use LWP::Simple;
use File::Temp qw/tempdir/;
use File::Spec;

use XML::Simple qw(:strict);
#use Data::Dumper;

open(IN, "<", $infile) or $LOGGER->logdie($!); 

my %GeneID; 

my @order;

my %GeneIDstatus;

my $tempdir = tempdir ( DIR => File::Spec->tmpdir(), CLEANUP=>1 );
#print $tempdir,"\n";

while(<IN>){
    chomp;
    next if ($_=~/^#/);
    my $ID=$_;
    getstore("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&rettype=full&retmode=xml&id=$ID", $tempdir.'/'.$ID.'.xml');

    my $xml  = XMLin($tempdir.'/'.$ID.'.xml', ForceArray => ['Gene-commentary'], KeyAttr    => {});
    
#    my $xml  = XMLin("/tmp/cEiZEQ3bo_/100576238.xml", ForceArray => ['Gene-commentary'], KeyAttr    => {});
#    my $xml  = XMLin("/tmp/mXE8iPy5Hx/406069.xml", ForceArray => ['Gene-commentary'], KeyAttr    => {});
    my $status = $xml->{'Entrezgene'}->{'Entrezgene_track-info'}->{'Gene-track'}->{'Gene-track_status'}->{'value'};
#    foreach my $hr (@{ $xml->{'Entrezgene'}->{'Entrezgene_locus'}->{'Gene-commentary'}->[0]->{'Gene-commentary_products'}->{'Gene-commentary'} }) {
#        print STDERR $hr->{'Gene-commentary_accession'},"\n";
#    }

    $GeneIDstatus{$ID} = $status;

    if ($status eq 'live') {
        $GeneID{$ID} = $ID;
    } elsif ($status eq 'discontinued') {
        #print $status,"\n";
        $GeneID{$ID} = '';
    } elsif ($status eq 'secondary') {
        my $newID;
        foreach my $dbtag (@{  $xml->{'Entrezgene'}->{'Entrezgene_track-info'}->{'Gene-track'}->{'Gene-track_current-id'}->{'Dbtag'} }){
            if ($dbtag->{'Dbtag_db'} eq 'GeneID') {
                $newID = $dbtag->{'Dbtag_tag'}->{'Object-id'}->{'Object-id_id'}."\n";
            }
        }
        $LOGGER->logdie("Not found new GeneID for $ID") unless ($newID);
        chomp($newID);
        $GeneID{$ID} = $newID;
    } else {
        $LOGGER->logdie("Status unknown: $status");
    }
    
    push(@order,$ID);

    sleep(5);
} 

close(IN); 

foreach my $oldGeneID ( @order ) {
    print $oldGeneID, "\t", $GeneID{$oldGeneID}, "\t", $GeneIDstatus{$oldGeneID}, "\n";
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
        -i      --infile    Input file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

