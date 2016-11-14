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
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use FileHandle;
use POSIX 'isatty';

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($WARN);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level,$infile);

#Usage("Too few arguments") if $#ARGV < 0;
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

my $fhin;

if ($infile) {

    $LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

    $fhin = FileHandle->new;
    $fhin->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fhin = \*STDIN;
    } else {
        $LOGGER->logdie("Missing input file (-i/--infile) or STDIN data");
    }
}


my $tmpdir = tempdir(".gfftempXXXXX", DIR => ".", CLEANUP => 1 );

#print $tmpdir,"\n";

open(GENE, ">", $tmpdir.'/'.'gene.gff') or $LOGGER->logdie($!);
open(RNA, ">", $tmpdir.'/'.'rna.gff') or $LOGGER->logdie($!);
open(PARENT, ">", $tmpdir.'/'.'parent.gff') or $LOGGER->logdie($!);


my %gene;
my %rna;

my %gene_count;
my %rna_count;
my %parent_count;

my %parent;

while(<$fhin>) {
    chomp;
    if ($_=~/^#/) {
        print $_,"\n";
        print GENE $_,"\n";
        print PARENT $_,"\n";
        print RNA $_,"\n";
    } else {
        my @F = split(/\t/, $_);
        my ($ID) = $F[8]=~/ID=([^;]+)/;

        if ($F[2] eq 'gene') {
            print GENE $_,"\n";
            unless (exists $gene{ $ID }) {
                $gene_count{$ID} = 0;
            }
            $gene{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] } = { 'data'=>\@F, 'name'=>$ID.'_'.++$gene_count{$ID} };
        } elsif (($F[2] eq 'RNA')||($F[2] eq 'miRNA')) {
            print RNA $_,"\n";
            unless (exists $rna{ $ID }) {
                $rna_count{$ID} = 0;
            } 
            $rna{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] } = { 'data'=>\@F, 'name'=>$ID.'_'.++$rna_count{$ID} };
        } else {
            my ($Parent) = $F[8]=~/Parent=([^;]+)/;
            if (exists $gene{ $Parent }) {
                $LOGGER->logdie("Found a gene parent that is not an RNA/miRNA ($Parent) [$_]");
            }
            unless (exists $parent{ $ID }) {
                $parent_count{$ID} = 0;
            }
            $parent{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] } = { 'data'=>\@F, 'name'=>$ID.'_'.++$parent_count{$ID} };
            print PARENT join("\t", @{ $parent{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] }->{'data'} }),"\n";
        }
    }        
}

$fhin->close();

close(GENE);
close(PARENT);
close(RNA);

my $intergff = "$tmpdir/inter_gene_rna.txt";

unless (-e $intergff) {
    # sem -s (para pegar erros de strand entre coordenadas de genes e parents
    `intersectBed -b $tmpdir/gene.gff -a $tmpdir/rna.gff -wao > $intergff`;
}

$|=1;
foreach my $g (keys %gene) {
    foreach my $chr (keys %{ $gene{ $g } }) {
        foreach my $strand (keys %{ $gene{ $g }->{ $chr } }) {
            foreach my $s (keys %{ $gene{ $g }->{ $chr }->{ $strand } }) {
                foreach my $e (keys %{ $gene{ $g }->{ $chr }->{ $strand }->{ $s } }) {
                    my $name = $gene{ $g }->{ $chr }->{ $strand }->{ $s }->{ $e }->{'name'};
                    my $ar_data = $gene{ $g }->{ $chr }->{ $strand }->{ $s }->{ $e }->{'data'};

                    $ar_data->[8]=~s/ID=[^;]+/ID=$name/;
                    print join("\t", @{ $ar_data }),"\n";
                }
            }
        }
    }
}    

open(INTER, "<", $intergff) or $LOGGER->logdie($!);

while(<INTER>) {
    chomp;
    my @F=split(/\t/, $_);
    my ($ID) = $F[8]=~/ID=([^;]+)/;
    $LOGGER->logdie("Not found ID for an entry ($F[8]) [".$_."]") unless ($ID);
    my ($Parent) = $F[8]=~/Parent=([^;]+)/;
    $LOGGER->logdie("Not found Parent for an entry ($F[8]) [".$_."]") unless ($Parent);
    if ($F[17]) {
        my ($InterID) = $F[17]=~/ID=([^;]+)/;
        $LOGGER->logdie("Not found ID for gene entry ($F[17]) [".$_."]") unless ($InterID);

        if ( $Parent eq $InterID ) {
            my ($chr, $strand, $s, $e) = @F[9, 15, 12, 13];

            my $new_name = $gene{ $InterID }->{ $chr }->{ $strand }->{ $s }->{ $e }->{'name'};
            $LOGGER->logdie("Not found new name $InterID") unless ($new_name);
            $F[8]=~s/Parent=[^;]+/Parent=$new_name/;
            
            my $new_rna_name = $rna{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] }->{'name'};
            $F[8]=~s/ID=[^;]+/ID=$new_rna_name/;

            if ($strand ne $F[6]) {
                $LOGGER->logwarn("Distinct strand found for gene ($InterID) and RNA ($ID). Changed to gene strand!");
                $LOGGER->logdie("Not found data for RNA $ID") unless ($rna{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] }->{'data'});
                $rna{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] }->{'data'}->[6] = $strand;
                $F[6]=$strand;
            }

            print join("\t", @F[0..8]),"\n";
        }            
    } else {
        $LOGGER->logdie("Not found an intersection for RNA entry ($_)");
    }
}

close(INTER);

$intergff = "$tmpdir/inter_rna_parent.txt";

unless (-e $intergff) {
    # sem -s (para pegar erros de strand entre coordenadas de genes e parents
    `intersectBed -b $tmpdir/rna.gff -a $tmpdir/parent.gff -wao > $intergff`;
}

open(INTER, "<", $intergff) or $LOGGER->logdie($!);

while(<INTER>) {
    chomp;
    my @F=split(/\t/, $_);
    my ($ID) = $F[8]=~/ID=([^;]+)/;
    $LOGGER->logdie("Not found ID for an entry ($F[8]) [".$_."]") unless ($ID);
    my ($Parent) = $F[8]=~/Parent=([^;]+)/;
    $LOGGER->logdie("Not found Parent for an entry ($F[8]) [".$_."]") unless ($Parent);
    if ($F[17]) {
        my ($InterID) = $F[17]=~/ID=([^;]+)/;
        $LOGGER->logdie("Not found ID for gene entry ($F[17]) [".$_."]") unless ($InterID);

        if ( $Parent eq $InterID ) {
            my ($chr, $strand, $s, $e) = @F[9, 15, 12, 13];

            my $new_name = $rna{ $InterID }->{ $chr }->{ $strand }->{ $s }->{ $e }->{'name'};
            $LOGGER->logdie("Not found new name $InterID") unless ($new_name);
            $F[8]=~s/Parent=[^;]+/Parent=$new_name/;
            
            my $new_parent_name = $parent{ $ID }->{ $F[0] }->{ $F[6] }->{ $F[3] }->{ $F[4] }->{'name'};
            $F[8]=~s/ID=[^;]+/ID=$new_parent_name/;
            
            $LOGGER->logdie("Not found RNA info for $Parent") unless ($rna{ $Parent }->{ $F[9] }->{ $F[15] }->{ $F[12] }->{ $F[13] }->{'data'}->[6]);
            if ( $rna{ $Parent }->{ $F[9] }->{ $F[15] }->{ $F[12] }->{ $F[13] }->{'data'}->[6] ne $F[6]) {
                my ($gene) = $F[17]=~/Parent=([^;]+)/;
                $LOGGER->logwarn("Distinct strand found for gene ($gene) of rna ($InterID) and EXON/CDS ($ID). Changed to gene strand!");
                $F[6]=$rna{ $Parent }->{ $F[9] }->{ $F[15] }->{ $F[12] }->{ $F[13] }->{'data'}->[6];
            }

            print join("\t", @F[0..8]),"\n";
        }            
    } else {
        $LOGGER->logdie("Not found an intersection for EXON/CDS entry ($_)");
    }
}

close(INTER);

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
        -l      --level     Log level [Default: WARN]
        -i      --infile    Input file (GFF) - The GFF needs to contain gene/RNA/exon|CDS structure, 
                            where every exon or CDS entry must have an RNA entry and every rna entry 
                            must have a gene entry. We don't accept other entries than gene, RNA, miRNA, 
                            exon or CDS. [Default: STDIN]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

