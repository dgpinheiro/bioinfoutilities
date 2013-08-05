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

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level);

my ($infile, $orderfile, $rnafile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "o|orderfile=s"=>\$orderfile,
            "i|infile=s"=>\$infile,
            "r|rnafile=s"=>\$rnafile
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

$LOGGER->logdie("Missing GenBank input file") unless ($infile);
$LOGGER->logdie("Wrong GenBank input file. File (\"$infile\") not found.") unless (-e $infile);

$LOGGER->logdie("Missing GenBank rna file") unless ($rnafile);
$LOGGER->logdie("Wrong GenBank rna file. File (\"$rnafile\") not found.") unless (-e $rnafile);


my $in = Bio::SeqIO->new(-file => $infile,-format=>'GenBank') or $LOGGER->logdie($!);

my %translate = (   'db_xref'=>'Dbxref',
                    'note'=>'Note'
                );

my %chr;
{
    open(ORDER, '<', $orderfile) or $LOGGER->logdie($!);

    while(<ORDER>) {
        chomp;
        my ($ncbi, $group) = split(/\t/, $_);
        $chr{$ncbi} = $group;
    }

    close(ORDER);
}

my %RNA;
{
    my $rnain = Bio::SeqIO->new(-file=>$rnafile) or $LOGGER->logdie ($!);
    my $gn_counter = 0;    
    my $pn_counter = 0;
    while(my $seq = $rnain->next_seq()) {
        my $transcript_id=$seq->display_id().'.'.$seq->version();

        my $productName_name = $seq->desc();
        my %productName;

        my $protein_id;

        my $geneName_name;
        my %geneName;

        my $locusLinkId;
        foreach my $feat ( $seq->get_SeqFeatures() ) {
            if ($feat->primary_tag() =~/gene/) {
                $LOGGER->logdie('Not found /gene value') unless ($feat->has_tag('gene'));
                ($geneName_name) = $feat->each_tag_value('gene');
                $LOGGER->logdie('Not found /db_xref value') unless ($feat->has_tag('db_xref'));

                foreach my $v ( $feat->each_tag_value('db_xref') ) {
                    if ($v =~/GeneID:(\d+)/){
                        $locusLinkId = $1;
                        last;
                    }
                }
                $LOGGER->logdie('Not found GeneID in /db_xref value(s)') unless ($locusLinkId);
            }
            elsif ($feat->primary_tag() =~/CDS/) {
                if ($feat->has_tag('protein_id')) {
                    ($protein_id) = $feat->each_tag_value('protein_id');
                    $RNA{$protein_id} = $transcript_id;

                    $LOGGER->logdie('Not found /protein_id for '.$transcript_id) unless ($protein_id);
                }
            }
        }
        unless (exists $geneName{$geneName_name}) {
            $geneName{$geneName_name} = $gn_counter;
        }
        unless (exists $productName{$productName_name}) {
            $productName{$productName_name} = $pn_counter;
        }
    }
}



my %data = ('count'=>{}, 'geneid'=>{}, 'locs'=>{});
my @full;
while(my $seq = $in->next_seq()) {

	foreach my $feat ( $seq->get_SeqFeatures() ) {
		if ($feat->primary_tag()=~/gene/) {
#			print ">>>>",$feat->{'_gsf_tag_hash'}->{'gene'}->[0],"\n";
			$data{'count'}->{ $feat->{'_gsf_tag_hash'}->{'gene'}->[0] }++;
			if ($feat->has_tag('db_xref')) {
				my $geneid;
				foreach my $db_xref ( $feat->each_tag_value('db_xref') ) {
					if ($db_xref=~/GeneID:(\d+)/) {
						$geneid = $1;
						last;
					}
				}
				
				if ($geneid) {
					$data{'geneid'}->{$feat->{'_gsf_tag_hash'}->{'gene'}->[0]}->{ $geneid } = undef;
                    $data{'locs'}->{ $geneid }->{ $seq->display_id().'.'.$seq->version() } = undef;
				}
				else {
					$LOGGER->logdie("Not found GeneID for gene ".$feat->{'_gsf_tag_hash'}->{'gene'}->[0]);
				}

			}
			else {
				$LOGGER->logdie("Not found db_xref for gene ".$feat->{'_gsf_tag_hash'}->{'gene'}->[0]);
			}
		}
	}
	push(@full, $seq);
}

#my %counter;

my %gene;
my %mRNA;
my %exon;
my %CDS;
my %tcounter;
my %ccounter;
foreach my $seq ( @full ) {

    my @comments = $seq->annotation->get_Annotations('comment');

    my $ncbi_id = $seq->display_id().'.'.$seq->version();

    foreach my $c (@comments) {
        if ($c->{'text'}=~/identical to (\S+)\./) {
            $ncbi_id = $1;
            last;
        }
    }
    
	foreach my $feat ( $seq->get_SeqFeatures() ) {
		if ($feat->primary_tag()=~/gene/) {

			my $gene = $feat->{'_gsf_tag_hash'}->{'gene'}->[0];

            $ccounter{$gene} = 1;
            $tcounter{$gene} = 1;

			if ($data{'count'}->{$gene} > 1) {
#				$counter{ $gene } = 0 unless (exists $counter{ $gene });
#				$counter{ $gene }++;
#				$feat->{'_gsf_tag_hash'}->{'gene'}->[0] = $feat->{'_gsf_tag_hash'}->{'gene'}->[0].'#'.$counter{ $feat->{'_gsf_tag_hash'}->{'gene'}->[0] };

				my $geneid;
				foreach my $db_xref ( $feat->each_tag_value('db_xref') ) {
					if ($db_xref=~/GeneID:(\d+)/) {
						$geneid = $1;
					}
				}

				unless (exists $data{'geneid'}->{$gene}->{$geneid}) {
					$LOGGER->logdie("Something wrong! The expected geneid $geneid and gene $gene don't match.");
				}
                
				$feat->{'_gsf_tag_hash'}->{'gene'}->[0] = $gene.'#'.$geneid;

                if ( scalar(keys %{$data{'locs'}->{$geneid}}) > 1 ) {
				    $feat->{'_gsf_tag_hash'}->{'gene'}->[0] = $gene.'#'.$geneid.'|'.$chr{$ncbi_id};
                }
	
				$data{'geneid'}->{$gene}->{ $geneid } = $feat->{'_gsf_tag_hash'}->{'gene'}->[0];

                $gene = $feat->{'_gsf_tag_hash'}->{'gene'}->[0];
			}
            
            my ($origgene) = $gene=~/([^#]+)(?:#\S+)?$/;
            my @complement = ('Name='.$origgene);
            foreach my $k (keys %{$feat->{'_gsf_tag_hash'}}) {
                next if ($k=~/gene/);

                my @v = map { if ($_) { s/;/,/g; } $_ } $feat->each_tag_value($k);

                push(@complement, ((exists $translate{$k}) ? $translate{$k} : $k).'='.((scalar(@v)) ? join(',', @v ) : '_no_value'));
            }

            $gene{$gene}->{'data'} = [$chr{$ncbi_id}, 'GenBank', 'gene', $feat->start(), $feat->end(), '.', ( ($feat->strand()) ?  ( ($feat->strand()>0) ? '+' : '-') : '.' ), '.', 'ID='.$feat->{'_gsf_tag_hash'}->{'gene'}->[0].';'.join(';', @complement) ];
            
#            print STDERR join("\t", @{$gene{$gene}->{'data'}}),"\n";

		}
		else {
			if ($feat->has_tag('gene')) {
				my $gene = $feat->{'_gsf_tag_hash'}->{'gene'}->[0];
				
				if ($data{'count'}->{$gene} > 1) {
				
					my $geneid;
					foreach my $db_xref ( $feat->each_tag_value('db_xref') ) {
						if ($db_xref=~/GeneID:(\d+)/) {
							$geneid = $1;
							last;
						}
					}

					unless (exists $data{'geneid'}->{$gene}->{$geneid}) {
						$LOGGER->logdie("Something wrong! The expected geneid $geneid and gene $gene don't match.");
					}

					$feat->{'_gsf_tag_hash'}->{'gene'}->[0] = $data{'geneid'}->{$gene}->{ $geneid };
                    $gene = $feat->{'_gsf_tag_hash'}->{'gene'}->[0];
				}

                if ($feat->primary_tag()=~/RNA|exon/) {
                    
                    my $transcript_id = $feat->{'_gsf_tag_hash'}->{'transcript_id'}->[0];
                    $transcript_id||=$gene.'#T'.++$tcounter{$gene};

                    my ($origgene) = $gene=~/([^#]+)(?:#\S+)?$/;

                    my @complement = ('Name='.$origgene);
                    foreach my $k (keys %{$feat->{'_gsf_tag_hash'}}) {
                        next if ($k=~/gene|transcript_id/);

                        my @v = map { if ($_) { s/;/,/g; } $_ } $feat->each_tag_value($k);

                        push(@complement, ((exists $translate{$k}) ? $translate{$k} : $k).'='.((scalar(@v)) ? join(',', @v ) : '_no_value'));
                    }
                    
                    my $origtranscript_id=$transcript_id;

                    if ( $gene=~/\|([^\|]+)$/) {
                        $LOGGER->logdie("The $1 and $chr{$ncbi_id} do not match for $transcript_id") unless ($1 eq $chr{$ncbi_id});
                        $transcript_id = $transcript_id.'|'.$chr{$ncbi_id};
                    }

                    $mRNA{$transcript_id}->{'data'} = [$chr{$ncbi_id}, 'GenBank', 'mRNA', $feat->start(), $feat->end(), '.', ( ($feat->strand()) ?  ( ($feat->strand()>0) ? '+' : '-') : '.' ), '.', 'ID='.$transcript_id.';'.'Parent='.$gene.';'.join(';', @complement) ];
                
#                    print STDERR join("\t", @{$mRNA{$transcript_id}->{'data'}}),"\n";
                    
                    push(@{$gene{$gene}->{'children'}}, $transcript_id);
                    
                    if ($feat->location->isa("Bio::Location::SplitLocationI")) {
                        my $exon_count = 1;
                        foreach my $l ($feat->location->each_Location()) {
            
                            my $exon_id = $transcript_id.'#E'.$exon_count;
                            $exon{$exon_id}->{'data'} = [$chr{$ncbi_id}, 'GenBank', 'exon', $l->start()+0, $l->end()+0, '.', ( ($l->strand()) ?  ( ($l->strand()>0) ? '+' : '-') : '.' ), '.', 'ID='.$exon_id.';'.'Parent='.$transcript_id ];

                            push(@{$mRNA{$transcript_id}->{'children-exon'}}, $exon_id);
                    
#                            print STDERR join("\t", @{$exon{$exon_id}->{'data'}}),"\n";
                            $exon_count++;                            
                        }
                    }


                }
                elsif ($feat->primary_tag()=~/CDS/) {
                    my $protein_id = $feat->{'_gsf_tag_hash'}->{'protein_id'}->[0];
                    unless ($protein_id) {
                        $protein_id = $gene.'#P'.++$ccounter{$gene};
                    }

                    my $transcript_id;
                    if (exists $RNA{$protein_id}) {
                        $transcript_id=$RNA{$protein_id};
                    }
                    else {
                        $transcript_id=$gene.'#T'.($ccounter{$gene}-1);
                        warn "Setting transcript $transcript_id to $protein_id\n";
                        $mRNA{$transcript_id}->{'data'} = [$chr{$ncbi_id}, 'GenBank', 'mRNA', $feat->start(), $feat->end(), '.', ( ($feat->strand()) ?  ( ($feat->strand()>0) ? '+' : '-') : '.' ), '.', 'ID='.$transcript_id.';'.'Parent='.$gene ];
                        push(@{$gene{$gene}->{'children'}}, $transcript_id);
                    }

                    my ($origgene) = $gene=~/([^#]+)(?:#\d+)?$/;

                    my @complement = ('Name='.$origgene);
                    foreach my $k (keys %{$feat->{'_gsf_tag_hash'}}) {
                        next if ($k=~/gene/);

                        my @v = map { if ($_) { s/;/,/g; } $_ } $feat->each_tag_value($k);
                        
                        if ($k ne 'translation') {
                            push(@complement, ((exists $translate{$k}) ? $translate{$k} : $k).'='.((scalar(@v)) ? join(',', @v ) : '_no_value'));
                        }                            
                    }
                    
                    my $origprotein_id = $protein_id;

                    if ( $gene=~/\|([^\|]+)$/) {
                        $LOGGER->logdie("The $1 and $chr{$ncbi_id} do not match for $protein_id") unless ($1 eq $chr{$ncbi_id});

                        $protein_id = $protein_id.'|'.$chr{$ncbi_id};
                        $transcript_id = $transcript_id.'|'.$chr{$ncbi_id};
                    }
                    
                    if ($feat->location->isa("Bio::Location::SplitLocationI")) {
                        my $cds_count = 1;
                        my $start_codon;

                        if ($feat->{'_gsf_tag_hash'}->{'start_codon'}->[0]) {
                            $start_codon = $feat->{'_gsf_tag_hash'}->{'start_codon'}->[0];
                        }
                        else {
                            $start_codon = 1;
                        }
                        foreach my $l ($feat->location->each_Location()) {
                            
                            my $frame = 0;
                            my $cds_length = ($l->end()-$l->start())+1;

                            my $cds_id = $protein_id.'#C'.$cds_count;
                            $CDS{$cds_id}->{'data'} = [$chr{$ncbi_id}, 'GenBank', 'CDS', $l->start()+0, $l->end()+0, '.', ( ($l->strand()) ?  ( ($l->strand()>0) ? '+' : '-') : '.' ), $frame, 'ID='.$cds_id.';'.'Parent='.$transcript_id.';'.join(';', @complement) ];
                            push(@{$mRNA{$transcript_id}->{'children-cds'}}, $cds_id);
#                            print STDERR join("\t", @{$cds{$cds_id}->{'data'}}),"\n";
                            $cds_count++;                            
                        }

                        my @frameorder;
                        # Plus (+)
                        if ($mRNA{$transcript_id}->{'data'}->[6] eq '+') {
                            @frameorder = @{$mRNA{$transcript_id}->{'children-cds'}};
                        }
                        else { # Minus (-)
                            @frameorder = ( sort { $CDS{$b}->{'data'}->[3] <=> $CDS{$a}->{'data'}->[3] } @{$mRNA{$transcript_id}->{'children-cds'}} );
                        }

                        foreach my $cds_id ( @frameorder ) {
                            my $frame = 0;
                            my ($s, $e) = ($CDS{$cds_id}->{'data'}->[3], $CDS{$cds_id}->{'data'}->[4]);
                            my $cds_length = ($e-$s)+1;

                                $frame = $start_codon-1;
                                my $aux = ($cds_length - $frame) % 3;
                                $start_codon  = 3-(($aux) ? $aux : 3)+1;
                            
                            $CDS{$cds_id}->{'data'}->[7] = $frame;
                        }

                    }
                    elsif (($feat->location->isa("Bio::Location::Simple"))||($feat->location->isa('Bio::Location::Fuzzy'))) {
                            my $l = $feat->location();
                            my $cds_count = 1;
                            my $frame = 0;
                            my $cds_length = ($l->end()-$l->start())+1;

                            my $cds_id = $protein_id.'#C'.$cds_count;
                            $CDS{$cds_id}->{'data'} = [$chr{$ncbi_id}, 'GenBank', 'CDS', $l->start()+0, $l->end()+0, '.', ( ($l->strand()) ?  ( ($l->strand()>0) ? '+' : '-') : '.' ), $frame, 'ID='.$cds_id.';'.'Parent='.$transcript_id.';'.join(';', @complement) ];
                            push(@{$mRNA{$transcript_id}->{'children-cds'}}, $cds_id);
                    }

                }

			}

		}
	}
#	$out->write_seq($seq);
}


print '##gff-version 3',"\n";
foreach my $g (keys %gene) {
    print join("\t", @{ $gene{ $g }->{'data'} }),"\n";
    foreach my $m (@{ $gene{ $g }->{'children'} }) {
        print join("\t", @{ $mRNA{ $m }->{'data'} }),"\n";
        if ((exists $mRNA{ $m })&&(exists $mRNA{ $m }->{'children-exon'})&&(scalar( @{ $mRNA{ $m }->{'children-exon'} } ) > 0)) {
            foreach my $e (@{ $mRNA{ $m }->{'children-exon'} }) {
                print join("\t", @{ $exon{ $e }->{'data'} }),"\n";
            }
        }
        else {
            my @exon = @{ $mRNA{ $m }->{'data'} };
            $exon[2] = 'exon';
            $exon[8]=~s/ID=([^;]+)/ID=$1#E1/;
            my $p = $1;
            $exon[8]=~s/Parent=[^;]+/Parent=$p/;
            print join("\t", @exon),"\n";
        }
        foreach my $c (@{ $mRNA{ $m }->{'children-cds'} }) {
            if ( ( exists $CDS{$c} ) && ( $CDS{$c}->{'data'} ) ) {
                print join( "\t", @{ $CDS{$c}->{'data'} } ), "\n";
                delete( $CDS{$c} );
            }
            else {
                warn 'We cannot find $CDS{"'.$c.'"} and $CDS{"'.$c.'"}->{"data"} from mRNA '.$m. ' and gene '.$g;
            }            
        }
    }
}

foreach my $c (keys %CDS) {
    print STDERR $c,"\n";
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
        -o      --orderfile Order file (ORDER.txt)
        -r      --rnafile   RNA file (rna.gbk)
        -i      --infil     Input file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

