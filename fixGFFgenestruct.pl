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

our %gene;
our %rna;

our %geneexon;

while(<$fhin>) {
    chomp;

    my @F = split(/\t/, $_);

	if ($_=~/^#/) { 
		print $_,"\n"; 
	} else { 
		if (($F[2] eq "gene")||($F[2] eq "pseudogene")) { 
			my ($ID)=$F[8]=~/ID=([^;]+)/;
            unless ($ID) {
                #$LOGGER->logdie("Missing ID: ".join("\t", @F)) unless ($ID);
                #($ID) = $F[8]=~/locus_tag=([^;]+(?:;part=[^;]+)?)/;
                my ($locus_tag) = $F[8]=~/locus_tag=([^;]+)/;
                $ID='gene-'.$locus_tag;
                #$ID=~s/;/_/g;
                #$ID=~s/\//_/g;
                $F[8]='ID='.$ID.';'.$F[8];
            }
            $F[2] = "gene"; # substitui pseudogene por gene
			push(@{ $gene{$ID} }, [@F]); 
		} elsif ( ($F[2]=~/(?:antisense_)?[tr]RNA/) && ($F[8]!~/Parent/) ) { 
			my $eight = $F[8]; 
			$F[8] =~ s/ID=([^;]+)/ID=gene_$1/; 
			my $id=$1; 
			$F[2] = "gene"; 
		
            print join("\t", @F),"\n";

			$F[2]="RNA"; 
			$F[8]=$eight; 
			$F[8].=";Parent=gene_$id"; 
		} elsif ($F[2]=~/\b((?:t?m|nc|sn|SRP_|lnc_|t|r|antisense_)RNA|RNASE_P_RNA|telomerase_RNA|RNase_MRP_RNA|guide_RNA|primary_transcript|transcript|\w?_gene_segment)\b/) { 
            my ($rna_track) = $1;
			$F[2]="RNA"; 
			my ($ID)=$F[8]=~/ID=([^;]+)/;
			my ($Parent)=$F[8]=~/Parent=([^;]+)/;
            unless ($Parent) { # RNA without gene Parent information
                my @G = @F;
                #$LOGGER->logdie("Missing Parent: ".join("\t",@F));
                $Parent='gene_'.$ID;
                $F[8].=';Parent='.$Parent;
                $G[2]='gene';
                $G[8]=~s/ID=[^;]+/ID=$Parent/;
                push(@{ $gene{$Parent} }, [@G]);
                print join("\t", @G),"\n";
            }
            
			if (exists $gene{$Parent}) { 
				for (my $i=0; $i<= $#{$gene{$Parent} }; $i++) {
                    if ($gene{$Parent}->[$i]) {
                        if ( ($F[3] == $gene{$Parent}->[$i]->[3]) && ($F[4] == $gene{$Parent}->[$i]->[4]) ) {
                            $gene{$Parent}->[$i]=undef;
                        }
                    }                        
				}
			}
			
			push(@{$rna{$ID}}, [@F]);

		} elsif ($F[2]=~/\b(CDS|exon)\b/) {  
			my ($Parent)=$F[8]=~/Parent=([^;]+)/; 
			next unless ($Parent);
            if (exists $gene{$Parent}) {
                my $RID='rna_'.$Parent;
				$F[8]=~s/Parent=([^;]+)/Parent=$RID/;
                if ($F[2] eq 'exon') {
                    $geneexon{$RID} = 1;
                }
			}
			if ($F[2] eq "exon") {
				delete($rna{$Parent});
			}
		} elsif ($F[2]=~/\b(miRNA)\b/) { 
			$F[2]=$1; 
		}
		print join("\t", @F),"\n"; 
	} 
}

$fhin->close();

foreach my $k (keys %gene) { 
    if (defined $gene{$k}) {
        foreach my $ar (@{ $gene{$k} }) {
            if ($ar) {
                $ar->[2]="RNA"; 
                $ar->[8]=~s/ID=([^;]+)/ID=rna_$1/; 
                $ar->[8].=";Parent=$k"; 
                print join("\t", @{$ar}),"\n"; 
                unless (exists $geneexon{'rna_'.$1}) {
                    $ar->[2]="exon"; 
                    $ar->[8]=~s/ID=([^;]+)/ID=exon_$1/; 
                    $ar->[8]=~s/Parent=([^;]+)/Parent=rna_$1/; 
                    print join("\t", @{$ar}),"\n";
                }
            }
        }
    } 
}

foreach my $k (keys %rna) { 
    foreach my $ar (@{ $rna{$k} }) {
        $ar->[2]="exon"; 
        $ar->[8]=~s/ID=([^;]+)/ID=exon_$1/; 
        $ar->[8]=~s/Parent=([^;]+)/Parent=rna_$1/; 
        print join("\t", @{ $ar }),"\n";
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

        -h      --help      Help
        -l      --level     Log level [Default: WARN]
        -i      --infile    Input file (GFF) - The GFF needs to contain only gene, RNA, exon or 
                            CDS features. [Default: STDIN]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}
