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

    getNCBIseqs.pl 

=head1 SYNOPSIS

    perl getNCBIseqs.pl -d nuccore -t 'txid7460[Organism:exp]+AND+srcdb+refseq[Properties]+AND+biomol+RNA[Properties]' -o AmRefSeqRNA.fa

=head1 ABSTRACT

    Get sequences from NCBI using eutils (http://eutils.ncbi.nlm.nih.gov).    

=head1 DESCRIPTION
    
    Arguments:

        -h/--help           Help
        -d/--db             Database name (e.g. nuccore)
        -t/--term           Query term
        -o/--outfile        Output file path
        -e/--outfiletype    Output file type (e.g. fasta, genbank) DEFAULT: fasta
        -r/--retmax         Maximum sequence retrieval DEFAULT: 1

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2012 Universidade de São Paulo

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Getopt::Long;
use FileHandle;
use LWP::Simple;
use LWP::UserAgent;
use File::Temp qw/ tempfile tempdir /;
use File::Spec;
use Bio::SeqIO;

my ($query, $outfile, $db, $outftype, $retmax);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"          => sub { &Usage(); },
            "t|term=s"          => \$query,
            "o|outfile=s"       => \$outfile,
            "d|db=s"            => \$db,
            "e|outfiletype=s"   => \$outftype,
            "r|retmax=i"        => \$retmax
 ) or &Usage();

die "Missing term value for the query" unless ($query);
die "Missing db" unless ($db);
die "Missing output file" unless ($outfile);

$outftype||='fasta';

#assemble the esearch URL
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";


#post the esearch URL
print "ESEARCH URL ($url)...\n";

my $output = get($url);


#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

print "FOUND $count records...\n";

#open output file for writing
my $fhout = FileHandle->new(">$outfile");
unless ( defined $fhout ) {
    die "Can't open file ($outfile)!\n"
}
$fhout->autoflush(1);

#retrieve data in batches of 1 by default
$retmax||= 1;
#maximum tries 10
my $maxtries = 10;
for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=$outftype&retmode=text";
        my $try = 1;
        TRY:
        print " # Try $try: EFETCH URL $efetch_url\n";
        my $tempdir = tempdir ( 'getNCBIseqsXXXX', DIR => File::Spec->tmpdir());
        getstore($efetch_url,"$tempdir/$web");
        open(IN, "<", "$tempdir/$web") or die $!;      
        my $efetch_out='';
        while (<IN>) {
            $efetch_out.=$_;
        }
        close(IN);
        if ((!defined $efetch_out)||((defined $efetch_out)&&($efetch_out=~/<ERROR>([^>]+)<\/ERROR>/))) {
            print "Found an error: $1\t";
            if (!$efetch_out) {
                if ($try <= $maxtries) {
                    print "NEXT TRY ($try) ...\n";
                    $try++;
                    sleep(5);
                    goto TRY;
                }
                else {
                    die "Sorry, maximum tries ($maxtries) exceeded!";
                }
            }
        }
        else {
            my $seqio =  Bio::SeqIO->new(-file => "$tempdir/$web");
            my $c = 0;
            while ( my $seq = $seqio->next_seq() ) {
                $c++;
            }
            my $ec = ((($count-$retstart)>$retmax) ? $retmax : ($count-$retstart));
            if ($c!=$ec) {
                die "Some sequences not retrieved (Expected: $ec/Retrieved: $c)";
            }
        }

        print $fhout "$efetch_out";
        
}
$fhout->close;
        
my $seqio = Bio::SeqIO->new( -file => $outfile );
my $c = 0;
while ( my $seq = $seqio->next_seq() ) {
    $c++;
}
if ( $c != $count ) {
    die "Some sequences not retrieved (Expected: $count/Retrieved: $c)";
}



# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] -d <DB> -t <TERM> -o <OUTFILE> [-e <OUTFILETYPE>] [-r 10]

Argument(s)

        -h      --help          Help
        -d      --db            NCBI database (e.g. nuccore)
        -t      --term          Term (query)
        -o      --outfile       Output file
        -e      --outfiletype   output file type (e.g. fasta, genbank) default: fasta
        -r      --retmax        Maximum sequence retrieval default: 1

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

