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

my ($level, $infile, $gbsdir);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "g|gbsdir=s"=>\$gbsdir
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

use Bio::SeqIO;

use Storable;

my %rename;
unless (-e "chrnameinfo.dump") {
    foreach my $gbsfile (glob("$gbsdir/*.gbs")) {
        my $in = Bio::SeqIO->new(-file => $gbsfile,-format=>'GenBank') or $LOGGER->logdie($!);
        
        while(my $seq = $in->next_seq()) {
            my ($id, $desc) = ($seq->id().'.'.$seq->version(), $seq->desc());
            
            #print $id,"\t",$desc,"\n";

            my ($chr) = $desc =~/LG(\d+)/;

            unless ($chr) {
                ($chr) = $desc =~/(GroupUn\d+)/;

                unless ($chr) {
                    if ($desc=~/mitochondrion/) {
                        $rename{$id} = "MT";
                    } else {
                        $LOGGER->logdie("Not found chr ID");
                    }
                } else {
                    $rename{$id} = $chr; 
                }  
            } else {
                $rename{$id} = "Group$chr";
            }
        }
    }
    store \%rename, 'chrnameinfo.dump';
} else {
   %rename = %{ retrieve('chrnameinfo.dump') };
}


open(GFF, "<", $infile) or $LOGGER->logdie($!);

my %track;
my %trackID;

my %create_gene;
my %create_rna;
my %gene;
my %rnaparent;

while(<GFF>) {
    chomp;
    if ($_=~/^#/) {
	if ($_=~/##sequence-region (\S+)/) {
		my $nn= $rename{ $1 };
		$_=~s/##sequence-region \S+/##sequence-region $nn/;
	} 

        print $_,"\n";
	
    } else {
        my (@data) = split(/\t/, $_);
        
        next if (($data[2] eq 'region')||($data[2] =~ /match/));
        my ($id) = $data[8]=~/ID=([^;]+)/;
        my ($name) = $data[8]=~/Name=([^;]+)/;
        my ($parent) = $data[8]=~/Parent=([^;]+)/;
        my ($product) = $data[8]=~/product=([^;]+)/;
        my ($target) = $data[8]=~/Target=([^;]+)/;
        
        my $newparent;

        if ($parent) {
            my $lastparent;
            foreach my $s (keys %{ $track{$parent} }) {
                foreach my $e (keys %{ $track{$parent}->{$s} }) {
                    if (($data[3] >= $s)&&($data[4] <= $e)) {
                        $newparent=$track{$parent}->{$s}->{$e};
                    }
                    $lastparent=$track{$parent}->{$s}->{$e};
                }
            }
            unless ($newparent) {
#                $LOGGER->logwarn("Not found new parent for $parent");
#                $newparent=$lastparent.'_x';
#                my @tmp = @data;
#                $tmp[8]=~s/ID=[^;]+/ID=$newparent/;
#                $tmp[2]='gene';
#                print join("\t", $rename{ $tmp[0] }, @tmp[1..8]),"\n";
                $parent = undef;
                $data[8]=~s/Parent=[^;]+;?//;
                                
            } else {
                $data[8]=~s/Parent=[^;]+/Parent=$newparent/;
            }                
        }

        unless ($name) {
            if ($parent) {
                $name=$newparent.'_';
            } else {
                $name='';
            }

            $name.=$id;
        }
            my $newname;
                if (($data[2] =~/RNA|transcript/i)&&($name=~/^rna\d+/)&&($product)) {
                    $newname=$product.'_'.$name;
                } elsif (($data[2] =~/match/)&&($target)) {
                    if ($target) {
                        my ($tgname)=$target =~/^(\S+)/;
                        my $tg = $rename{ $data[0] }.'_'.$tgname;
                        $target=~s/\s+/\_/g;
                        unless (exists $create_gene{'genetarget_'.$tg}) {
                            $create_gene{'genetarget_'.$tg} = {'chr'=>$rename{ $data[0] }, src=>$data[1], 'type'=>'gene', 'start'=>$data[3], 'end'=>$data[4], 'score'=>0, 'strand'=>$data[6], 'frame'=>'.', 'opt'=>'ID=genetarget_'.$tg};
                        } else {
                            if ($data[3] < $create_gene{'genetarget_'.$tg}->{'start'}) {
                                $create_gene{'genetarget_'.$tg}->{'start'} = $data[3];
                            }
                            if ($data[4] > $create_gene{'genetarget_'.$tg}->{'end'}) {
                                $create_gene{'genetarget_'.$tg}->{'end'} = $data[4];
                            }

                        }
                        unless (exists $create_rna{'rnatarget_'.$tg}) {
                            $create_rna{'rnatarget_'.$tg} = {'chr'=>$rename{ $data[0] }, src=>$data[1], 'type'=>'RNA', 'start'=>$data[3], 'end'=>$data[4], 'score'=>0, 'strand'=>$data[6], 'frame'=>'.', 'opt'=>'ID=rnatarget_'.$tg.';Parent=genetarget_'.$tg};
                            $rnaparent{'rnatarget_'.$tg} = 'genetarget_'.$tg;
                        } else {
                            if ($data[3] < $create_rna{'rnatarget_'.$tg}->{'start'}) {
                                $create_rna{'rnatarget_'.$tg}->{'start'} = $data[3];
                            }
                            if ($data[4] > $create_rna{'rnatarget_'.$tg}->{'end'}) {
                                $create_rna{'rnatarget_'.$tg}->{'end'} = $data[4];
                            }

                        } 
                        $parent='genetarget_'.$tg;
                        $newname=$target;
                        $data[2] = 'exon';
                        $data[8].=';Parent=rnatarget_'.$tg;
                    } else {
                        $LOGGER->logdie("Not found target for $id");
                    }
                } else {
                    $newname=$name;
                }
            
            unless (exists $trackID{$newname}) {
                $trackID{$newname} = 0;
                $track{$id}->{$data[3]}->{$data[4]} = $newname;
#
#
#
#                $trackID{$name} = 0;
#                $track{$id}->{$data[3]}->{$data[4]} = $name;
            } else {
                $trackID{$newname}++;
                $newname.='-'.$trackID{$newname};
                $track{$id}->{$data[3]}->{$data[4]} = $newname;
#
#
#                
#                $trackID{$name}++;
#                $name.='-'.$trackID{$name};                
#                $track{$id}->{$data[3]}->{$data[4]} = $name;
            }
               
        $data[8]=~s/ID=[^;]+/ID=$newname/;
        
        if ($data[2] eq 'gene') {
            $gene{$newname}->{'gff'} = [ $rename{ $data[0] }, @data[1..8] ];
        } elsif ( $data[2] =~ /RNA|transcript/ ) {
            my ($parent)=$data[8]=~/Parent=([^;]+)/;
            unless ($parent) {
                unless ($newname=~/rnatarget/) {
                    unless (exists $create_gene{'gene_'.$newname}) {
                        $create_gene{'gene_'.$newname} = {'chr'=>$rename{ $data[0] }, src=>$data[1], 'type'=>'gene', 'start'=>$data[3], 'end'=>$data[4], 'score'=>0, 'strand'=>$data[6], 'frame'=>'.', 'opt'=>'ID=gene_'.$newname};
                        $data[8].=';Parent=gene_'.$newname;
                        $parent='gene_'.$newname;
                    }
                }                
            }
            push(@{ $gene{$parent}->{'rna'} }, [ $rename{ $data[0] }, @data[1..8] ]) if ($parent);
            $rnaparent{$newname} = $parent;
        } elsif ( $data[2] =~/exon/) {
            my ($parent)=$data[8]=~/Parent=([^;]+)/;
            #die join("\t", $rename{ $data[0] }, @data[1..8]) unless ($parent);
            unless (exists $rnaparent{ $parent }) {
                
                unless (exists $create_rna{'rna_'.$parent}) {
                        $create_rna{'rna_'.$parent} = {'chr'=>$rename{ $data[0] }, src=>$data[1], 'type'=>'RNA', 'start'=>$data[3], 'end'=>$data[4], 'score'=>0, 'strand'=>$data[6], 'frame'=>'.', 'opt'=>'ID=rna_'.$parent.';Parent='.$parent};
                        $rnaparent{ 'rna_'.$parent } = $parent;
                        $parent='rna_'.$parent;
                }
                else {
                        if ($data[3] < $create_rna{ 'rna_'.$parent }->{'start'}) {
                            $create_rna{'rna_'.$parent} ->{'start'} = $data[3];
                        }
                        if ($data[4] > $create_rna{'rna_'.$parent}->{'end'}) {
                            $create_rna{'rna_'.$parent}->{'end'} = $data[4];
                        }
                        $parent='rna_'.$parent;
                        $LOGGER->logdie("Error not expected") unless ($rnaparent{$parent});
                }
                $data[8]=~s/Parent=[^;]+/Parent=$parent/;
            } 
            push(@{ $gene{$rnaparent{$parent}}->{'exon'} }, [ $rename{ $data[0] }, @data[1..8] ]);

        } elsif ( $data[2] =~/cds/i) {
            my ($parent)=$data[8]=~/Parent=([^;]+)/;

            unless (exists $rnaparent{ $parent }) {
                
                unless (exists $create_rna{'rna_'.$parent}) {
                        $create_rna{'rna_'.$parent} = {'chr'=>$rename{ $data[0] }, src=>$data[1], 'type'=>'mRNA', 'start'=>$data[3], 'end'=>$data[4], 'score'=>0, 'strand'=>$data[6], 'frame'=>'.', 'opt'=>'ID=rna_'.$parent.';Parent='.$parent};
                        $rnaparent{ 'rna_'.$parent } = $parent;
                        $parent='rna_'.$parent;
                }
                else {
                        if ($data[3] < $create_rna{ 'rna_'.$parent }->{'start'}) {
                            $create_rna{'rna_'.$parent} ->{'start'} = $data[3];
                        }
                        if ($data[4] > $create_rna{'rna_'.$parent}->{'end'}) {
                            $create_rna{'rna_'.$parent}->{'end'} = $data[4];
                        }
                        $parent='rna_'.$parent;
                        $LOGGER->logdie("Error not expected") unless ($rnaparent{$parent});
                }
                $data[8]=~s/Parent=[^;]+/Parent=$parent/;
            } 
            push(@{ $gene{$rnaparent{$parent}}->{'cds'} }, [ $rename{ $data[0] }, @data[1..8] ]);

#            if (exists $rnaparent{ $parent }) {
#                push(@{ $gene{$rnaparent{$parent}}->{'cds'} }, [ $rename{ $data[0] }, @data[1..8] ]);
#            } else {
#                push(@{ $gene{$parent}->{'cds'} }, [ $rename{ $data[0] }, @data[1..8] ]);
#            }                
        }
        
        print join("\t", $rename{ $data[0] }, @data[1..8]),"\n";
    }
}

close(GFF);

foreach my $k (keys %create_gene) {
    print join("\t", @{ $create_gene{$k} }{'chr', 'src', 'type', 'start', 'end', 'score',  'strand', 'frame', 'opt'}),"\n";
    $gene{$k}->{'gff'} = [ @{ $create_gene{$k} }{'chr', 'src', 'type', 'start', 'end', 'score',  'strand', 'frame', 'opt'} ];
}

foreach my $k (keys %create_rna) {
    print join("\t", @{ $create_rna{$k} }{'chr', 'src', 'type', 'start', 'end', 'score',  'strand', 'frame', 'opt'}),"\n";
    my ($parent) = $create_rna{$k}->{'opt'}=~/Parent=([^;]+)/;
    push(@{ $gene{$parent}->{'rna'} }, [ @{ $create_rna{$k} }{'chr', 'src', 'type', 'start', 'end', 'score',  'strand', 'frame', 'opt'} ]);
}

foreach my $geneid (keys %gene) {
    my $parent;

    unless ($gene{$geneid}->{'rna'}) {
        die ">>>$geneid" unless ($gene{$geneid}->{'gff'});
        my @tmp = @{ $gene{$geneid}->{'gff'} };
        $tmp[2] = 'RNA';
        $tmp[8]=~s/ID=([^;]+)/ID=RNA_$geneid/;
        $tmp[8].=';Parent='.$geneid;
        print join("\t", @tmp),"\n";
        $parent='RNA_'.$geneid;
    }

    unless ($gene{$geneid}->{'exon'}) {
        unless ($parent) {
            ($parent) = $gene{$geneid}->{'rna'}->[0]->[8]=~/ID=([^;]+)/;
        }

        if ($gene{$geneid}->{'cds'}) {
            my $exoncount=1;
            foreach my $cdsar (@{ $gene{$geneid}->{'cds'} } ) {
                my @tmp = @{ $cdsar };
                $tmp[2] = 'exon';
                $tmp[8]=~s/ID=([^;]+)/ID=EXON_$1/;
                $tmp[8]=~s/Parent=([^;]+)/Parent=$parent/;;
                print join("\t", @tmp),"\n";
            }                    
        }
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
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file
        -g      --gbsdir    gbs dir

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

