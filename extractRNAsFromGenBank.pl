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

    extractRNAsFromGenBank.pl

=head1 SYNOPSIS
        
    perl extractRNAsFromGenBank.pl -i AmRefSeqDNA.gbk -o ./output
    
=head1 ABSTRACT

    Extract sequences of 5' UTR, 3' UTR and CDS of mRNA sequences and sequences of non-coding
RNAs from GenBank chromosome files. The global statistics are displayed in STDERR output.

=head1 DESCRIPTION

    Arguments:

        -h/--help       Help
        -l/--level      Log level [Default: FATAL] 
            OFF
            FATAL
            ERROR
            WARN
            INFO
            DEBUG
            TRACE
            ALL
        -i/--infile     Input file path and name
        -o/--outbase    Output path and base name [Default: ./out]

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
use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}


my ($infile, $level, $outbase);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"      => sub { &Usage(); },
            "i|infile=s"    => \$infile,
            "o|outbase=s"   => \$outbase,
            "l|level=s"     => \$level
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
$LOGGER->logdie("Wrong input file") unless (-e $infile);

$outbase||='./emfgout';

my $in = Bio::SeqIO->new(-file=>$infile);
    
my %out;
$out{'3p'} =
  Bio::SeqIO->new( -file => '>' . $outbase . '3p.fa', -format => 'FASTA' );
$out{'5p'} =
  Bio::SeqIO->new( -file => '>' . $outbase . '5p.fa', -format => 'FASTA' );
$out{'CDS'} =
  Bio::SeqIO->new( -file => '>' . $outbase . 'CDS.fa', -format => 'FASTA' );
$out{'missedCDS'} =
  Bio::SeqIO->new( -file => '>' . $outbase . 'missedCDS.fa', -format => 'FASTA' );
$out{'chromosome'} =
  Bio::SeqIO->new( -file => '>' . $outbase . 'chromosome.fa', -format => 'FASTA' );
$out{'misc_RNA'} =
  Bio::SeqIO->new( -file => '>' . $outbase . 'misc_RNA.fa', -format => 'FASTA' );
$out{'ncRNA'} =
  Bio::SeqIO->new( -file => '>' . $outbase . 'ncRNA.fa', -format => 'FASTA' );
$out{'rRNA'} =
  Bio::SeqIO->new( -file => '>' . $outbase . 'rRNA.fa', -format => 'FASTA' );
    
my $tmpdir = tempdir( 'emfgXXXX', DIR=>'./', CLEANUP => 1 );

my %count = (
    'chromosome'  => 0,
    'gene'        => 0,
    'mRNA'        => 0,
    '5pCDS3pmRNA' => {'all'=>0, 'small'=>0},
    '5pCDSmRNA'   => {'all'=>0, 'small'=>0},
    'CDS3pmRNA'   => {'all'=>0, 'small'=>0},
    'CDSmRNA'     => {'all'=>0, 'small'=>0},
    'ncRNA'       => 0,
    'byGroup'     => {}
);

while ( my $seqobj = $in->next_seq()) {
    
    my ($group) = ($1||$2) if ($seqobj->description()=~/linkage group ([^, ]+)|(mitochondrion)/);
    
    $LOGGER->logdie("Not found linkage group")  unless ($group);

    @{ $count{'byGroup'}->{$group} }{ 'gene', 'mRNA', '5pCDS3pmRNA', '5pCDSmRNA', 'CDS3pmRNA', 'CDSmRNA', 'ncRNA' } = (0, 0, {'all'=>0, 'small'=>0}, {'all'=>0, 'small'=>0}, {'all'=>0, 'small'=>0}, {'all'=>0, 'small'=>0}, 0) unless (exists $count{'byGroup'}->{$group});
    
    $count{'chromosome'}++;
 
    $LOGGER->info("Extracting features from : ".$seqobj->display_id());
    my @features = $seqobj->get_SeqFeatures();
    
    my ($fh, $filename) = tempfile( 'emfgXXXX', DIR => $tmpdir);

    foreach my $feature ( @features ) {
            
            # Transcript
            if ( $feature->primary_tag() =~ /RNA/) {
                my ($transcript_id, @db_xref) = ('', 0);
                if ($feature->has_tag('transcript_id')) {
                    ($transcript_id) = $feature->each_tag_value('transcript_id');
                }
                if ($feature->has_tag('db_xref')) {
                    @db_xref = $feature->each_tag_value('db_xref');
                }
                my $geneid = '?';
                my $beebase = '?';
                foreach my $xref (@db_xref) {
                    if ($xref=~/GeneID:(\d+)/) {
                        $geneid = $1;
                    }
                    elsif ($xref=~/BEEBASE:(\S+)/) {
                        $beebase = $1;
                    }
                }
                my $rna_seq = $feature->spliced_seq();
                print $fh $geneid,"\t",$beebase,"\t",$feature->primary_tag(),"\t",$transcript_id,"\t",$rna_seq->seq(),"\t",$feature->strand(),"\t",join(';', @db_xref),"\n";
            }
            
            # CDS
            if ( $feature->primary_tag() eq 'CDS') {
                my ($protein_id, @db_xref) = ('', 0);
                if ($feature->has_tag('protein_id')) {
                    ($protein_id) = $feature->each_tag_value('protein_id');
                }
                if ($feature->has_tag('db_xref')) {
                    @db_xref = $feature->each_tag_value('db_xref');
                }
                my $geneid = '?';
                my $beebase = '?';
                foreach my $xref (@db_xref) {
                    if ($xref=~/GeneID:(\d+)/) {
                        $geneid = $1;
                    }
                    elsif ($xref=~/BEEBASE:(\S+)/) {
                        $beebase = $1;
                    }
                }
                my $cds_seq = $feature->spliced_seq();
                print $fh $geneid,"\t",$beebase,"\t",$feature->primary_tag(),"\t",$protein_id,"\t",$cds_seq->seq(),"\t",$feature->strand(),"\t",join(';', @db_xref),"\n";
            }
    }

    # Sorting features
    $LOGGER->info("\tSorting features in $filename to $filename.sort");
    `sort -k1,1g -k2,2 -k3,3r -n $filename > $filename.sort`;

    # Parsing features
    $LOGGER->info("\tParsing features in $filename.sort");
    open(FT, "<", "$filename.sort") or $LOGGER->logdie($!);

    my %RNA;
    my $lastgeneid = '';
    my @m;
    my %gene;
    my %mRNA;
    my %ncRNA;
    while (<FT>) {
        chomp;
        my ($geneid, $beebase, $type, $acc, $seq, $strand, $other) = split(/\t/, $_);

        $gene{$geneid} = 0 unless (exists $gene{$geneid});;
        $gene{$geneid}++;

        if ($lastgeneid ne $geneid) {
            @m = ();
            $lastgeneid = $geneid;
        }
        if ($type =~ /RNA/) {
            $RNA{$acc} = { 'seq'       =>  $seq,
                           'geneid'    =>  $geneid,
                           'beebase'   =>  $beebase,
                           'cds'       =>  undef,
                           'type'      =>  $type
                          };
            if ( $type eq 'mRNA' ) {
                push( @m, $acc );
                $mRNA{$acc} = [];
            }
            else {
                $ncRNA{$acc} = undef;
            }
        }
        else {
            my $found;
            foreach my $mRNA_acc (@m) {
                if ($RNA{$mRNA_acc}->{'seq'}=~/$seq/) {
                    $LOGGER->debug($mRNA_acc."\t".$acc."\n");
                    
                    $RNA{$mRNA_acc}->{'cds'}->{$acc} = {    'start' =>  $-[0]+1, 
                                                            'end'   =>  $+[0]   };

                    unless ($found) {
                        $found = 1;
                        push(@{$mRNA{$mRNA_acc}}, $acc);
                    }
                    else {
                        $LOGGER->logwarn("CDS ($acc) matches more than one mRNA");
                    }                    
                }
            }

            unless (defined $found) {
                $LOGGER->logwarn("Not found mRNA for this CDS ($acc). Trying to match the chromosome (".$seqobj->display_id().")");
                

                    $RNA{$seqobj->display_id().'.'.$acc} = {    'seq'       =>  (($strand==1) ? $seqobj->seq() : $seqobj->revcom()->seq()),
                                                                'geneid'    =>  $geneid||'-',
                                                                'beebase'   =>  $beebase||'-',
                                                                'cds'       =>  undef,
                                                                'type'      =>  'chromosome'
                                                    };
                
                    if ($RNA{$seqobj->display_id().'.'.$acc}->{'seq'}=~/$seq/) {
                        $RNA{$seqobj->display_id().'.'.$acc}->{'seq'}=$seq;
                        $RNA{$seqobj->display_id().'.'.$acc}->{'cds'}->{$acc} = {   'start' =>  0,
                                                                                    'end'   =>  length($seq)   };
                    }
                    else {
                        $LOGGER->logdie("Not found match for this CDS ($acc) in the chromosome (".$seqobj->display_id().")");
                    }
            }
        }
    }
    close(FT);
    
    $count{'gene'}+=scalar(keys %gene);
    $count{'byGroup'}->{$group}->{'gene'}+=scalar(keys %gene);
    $count{'mRNA'}+=scalar(keys %mRNA);       
    $count{'byGroup'}->{$group}->{'mRNA'}+=scalar(keys %mRNA);
    $count{'ncRNA'}+=scalar(keys %ncRNA);
    $count{'byGroup'}->{$group}->{'ncRNA'}+=scalar(keys %ncRNA);

    $LOGGER->info("\tWriting sequence objects");
    foreach my $RNA_acc (keys %RNA) {
        my $seqobj = Bio::Seq->new(   -seq => $RNA{$RNA_acc}->{'seq'},
                                      -display_id => $RNA_acc,
                                      -description=> "GeneID:$RNA{$RNA_acc}->{'geneid'}|BEEBASE:$RNA{$RNA_acc}->{'beebase'}"
                                            );
        
        
        my %valid;
        my %seqlength;

        if ($RNA{$RNA_acc}->{'type'} eq 'mRNA') {
            if (scalar( keys %{$RNA{$RNA_acc}->{'cds'}})) {
                foreach my $CDS_acc (keys %{$RNA{$RNA_acc}->{'cds'}}) {
                    
                    if ($RNA{$RNA_acc}->{'cds'}->{$CDS_acc}->{'start'} > 1) {
                        my $cdsobj = $seqobj->trunc(1, $RNA{$RNA_acc}->{'cds'}->{$CDS_acc}->{'start'}-1);
                        $cdsobj->display_id( $cdsobj->display_id().'|5p'.'|'.$CDS_acc );

                        $valid{'5p'}->{ $CDS_acc } = $cdsobj;
                    }
                    if ($RNA{$RNA_acc}->{'cds'}->{$CDS_acc}->{'end'} < $seqobj->length()) {
                        my $cdsobj = $seqobj->trunc($RNA{$RNA_acc}->{'cds'}->{$CDS_acc}->{'end'}+1, $seqobj->length());
                        $cdsobj->display_id( $cdsobj->display_id().'|3p'.'|'.$CDS_acc );
                        
                        $valid{'3p'}->{ $CDS_acc } = $cdsobj;
                    }
                     
                    my $cdsobj = $seqobj->trunc($RNA{$RNA_acc}->{'cds'}->{$CDS_acc}->{'start'}, $RNA{$RNA_acc}->{'cds'}->{$CDS_acc}->{'end'});
                    $cdsobj->display_id( $cdsobj->display_id().'|CDS'.'|'.$CDS_acc );
                    $valid{'CDS'}->{ $CDS_acc } = $cdsobj;
                    
                    if (!exists $seqlength{'CDS'}) {
                        $seqlength{'CDS'} = $cdsobj->length();
                    }
                    else {
                        $seqlength{'CDS'} = $cdsobj->length() if ($seqlength{'CDS'} > $cdsobj->length());
                    }
                }
                
            }
            else {
                    my $cdsobj = $seqobj;
                    $cdsobj->display_id( $cdsobj->display_id().'|missedCDS'.'|'.'-' );
                    $valid{'missedCDS'}->{ $RNA_acc } = $cdsobj;
            }
        }
        else {
            my $cdsobj = $seqobj;
            $cdsobj->display_id( $cdsobj->display_id().'|'.$RNA{$RNA_acc}->{'type'}.'|'.'-' );
            $valid{$RNA{$RNA_acc}->{'type'}}->{ $RNA_acc } = $cdsobj;
        }
        if (exists $valid{'CDS'}) {

            if ( exists $valid{'5p'} ) {
                if ( exists $valid{'3p'} ) {
                    $count{'5pCDS3pmRNA'}->{'all'}++;
                    
                    $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}->{'all'}++;
                    if ( $seqlength{'CDS'} < 60 ) {
                        $count{'5pCDS3pmRNA'}->{'small'}++;
                        $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}->{'small'}++;
                    }
                }
                else {
                    $count{'5pCDSmRNA'}->{'all'}++;
                    $count{'byGroup'}->{$group}->{'5pCDSmRNA'}->{'all'}++;
                    if ( $seqlength{'CDS'} < 60 ) {
                        $count{'5pCDSmRNA'}->{'small'}++;
                        $count{'byGroup'}->{$group}->{'5pCDSmRNA'}->{'small'}++;
                    }
                }
            }
            else {
                if ( exists $valid{'3p'} )
                {
                    $count{'CDS3pmRNA'}->{'all'}++;
                    $count{'byGroup'}->{$group}->{'CDS3pmRNA'}->{'all'}++;
                    if ( $seqlength{'CDS'} < 60 ) {
                        $count{'CDS3pmRNA'}->{'small'}++;
                        $count{'byGroup'}->{$group}->{'CDS3pmRNA'}->{'small'}++;
                    }
                }
                else {
                    $count{'CDSmRNA'}->{'all'}++;
                    $count{'byGroup'}->{$group}->{'CDSmRNA'}->{'all'}++;
                }
            }
        }
        

        foreach my $type (keys %valid) {
                foreach my $s (keys %{$valid{$type}}) {
                    print STDERR "$type" unless ($out{$type});
                    $out{$type}->write_seq( $valid{$type}->{$s} );
                } 
        }
    }
    
}


print STDERR "
*** GENERAL REPORT ***

Nº of groups scanned....: $count{'chromosome'}
Nº of genes found.......: $count{'gene'}
Nº of mRNAs found.......: $count{'mRNA'}
    5p CDS 3p...........: $count{'5pCDS3pmRNA'}->{'all'}
        Small...........: $count{'5pCDS3pmRNA'}->{'small'}
    5p CDS..............: $count{'5pCDSmRNA'}->{'all'}
        Small...........: $count{'5pCDSmRNA'}->{'small'}
    CDS 3p..............: $count{'CDS3pmRNA'}->{'all'}
        Small...........: $count{'CDS3pmRNA'}->{'small'}
    CDS.................: $count{'CDSmRNA'}->{'all'}
        Small...........: $count{'CDSmRNA'}->{'small'}
    ncRNA...............: $count{'ncRNA'}
";

print STDERR "
*** REPORT BY GROUP ***
";
foreach my $group (sort {$a cmp $b} keys %{$count{'byGroup'}} ) {
    print STDERR "
* $group

    Nº of genes found.......: $count{'byGroup'}->{$group}->{'gene'}
    Nº of mRNAs found.......: $count{'byGroup'}->{$group}->{'mRNA'}
        5p CDS 3p...........: $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}->{'all'}
            Small...........: $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}->{'small'}
        5p CDS..............: $count{'byGroup'}->{$group}->{'5pCDSmRNA'}->{'all'}
            Small...........: $count{'byGroup'}->{$group}->{'5pCDSmRNA'}->{'small'}
        CDS 3p..............: $count{'byGroup'}->{$group}->{'CDS3pmRNA'}->{'all'}
            Small...........: $count{'byGroup'}->{$group}->{'CDS3pmRNA'}->{'small'}
        CDS.................: $count{'byGroup'}->{$group}->{'CDSmRNA'}->{'all'}
            Small...........: $count{'byGroup'}->{$group}->{'CDSmRNA'}->{'small'}
        ncRNA...............: $count{'byGroup'}->{$group}->{'ncRNA'}
";

}



# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>] -i <INFILE> -o <OUTBASE>

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file path and name
        -o      --outbase   Output path and base name [Default: ./out]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}


