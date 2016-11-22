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

	extractFromGenomeByGFF.pl

=head1 SYNOPSIS

	perl extractFromGenomeByGFF.pl 	-g /data/files/Amel_4.5/BeeBase/ -f amel_OGSv3.1_with_UTR.gff3 -o ./output \
					-p OGS_
	

=head1 ABSTRACT

	Extract sequences of  5' UTR, 3' UTR and CDS of mRNA sequences from chromosome files in FASTA format using
GFF coordinates. The global statistics are displayed in STDERR output.

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
         -g/--genome	Genome directory
         -f/--gff	GFF file
         -o/--outdir	Output directory
         -p/--outprefix	Output filenames' prefix

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
#use Bio::Factory::SeqAnalysisParserFactory;
use Bio::Tools::GFF;
use Bio::SeqIO;
use Bio::Seq;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $gfffile, $genomedir, $outdir, $outprefix, $genomefile, $statsfile, $threepminsize);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "d|genomedir=s"=>\$genomedir,
            "g|genomefile=s"=>\$genomefile,
            "s|statsfile=s"=>\$statsfile,
            "f|gff=s"=>\$gfffile,
            "o|outdir=s"=>\$outdir,
            "p|outprefix=s"=>\$outprefix,
            "m|threepminsize=i"=>\$threepminsize
    ) or &Usage();

$threepminsize||=0;

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

use File::Temp qw/ tempfile tempdir /;

$LOGGER->logdie('Missing genome directory or genome file') if ((!$genomedir)&&(!$genomefile));

if ($genomefile) {
    $LOGGER->logdie('Wrong genome file ('.$genomefile.')') unless (-e $genomefile);

    my $in = Bio::SeqIO->new(
        -file   => $genomefile,
        -format => 'fasta'
    );

    
    if ($genomedir) {
        unless (-e $genomedir) {
            mkdir($genomedir);
            $LOGGER->logdie('Cannot create genome directory ('.$genomedir.')') unless (-d $genomedir);
        }
    }
    else {
        $genomedir = tempdir( CLEANUP => 1 );
    }
    my @f = glob("$genomedir/*.fa");
    if ( scalar(@f) < 1) {
        while ( my $seq = $in->next_seq() ) {
            my $out = Bio::SeqIO->new(
                -file   => '>'.$genomedir.'/'.$seq->display_id().'.fa',
                -format => 'fasta'
            );
            $out->write_seq($seq);
        }
    }
}
else {
    if ($genomedir) {
        $LOGGER->logdie('Wrong genome directory ('.$genomedir.')') unless (-d $genomedir);
    }
}

$LOGGER->logdie('Missing gff file') unless ($gfffile);
$LOGGER->logdie('Wrong gff file ('.$gfffile.')') unless (-e $gfffile);

if ($outdir) {
    $LOGGER->logdie('Wrong output directory ('.$outdir.')') unless (-d $outdir);
}
else {
    $outdir='.';
}

$outprefix||='out_';

my %feature;
my %parent;

my $gffio = Bio::Tools::GFF->new(-file =>$gfffile, -gff_version => 3);
my $feat;
# loop over the input stream
my $c = 0;
my %chromosome;
my $idcounter = 0;
my %allgene;
while ($feat = $gffio->next_feature()) {
    $c++;
    die "Not found ".$feat->seq_id()." on $genomedir"  unless (-e "$genomedir/".$feat->seq_id().'.fa');

    my ($id,$parent);

    if ($feat->has_tag('ID')) {
        ($id) = $feat->get_tag_values('ID');
    }
    else {
        ($id) = 'UNK'.$idcounter++;
    }

        unless (exists $feature{$id}) {
            $feature{$id} = {   'primary_tag'   => $feat->primary_tag(),
                                'child'         => undef };
            push(@{$chromosome{ $feat->seq_id() }}, $id);
        }                            
    

    if ($feat->has_tag('Parent')) {
        ($parent) = $feat->get_tag_values('Parent');
        if ($id) {
            $parent{$id} = $parent unless (exists $parent{$id});
        }    

        push(@{$feature{$parent}->{'child'}->{ $feat->primary_tag() }}, $feat);
    }
    else {
        $allgene{ $feat->seq_id() }->{$id} = $feat;
    }
    
}

my %out;
my %count = (
    'scaffold'    => 0,
    'gene'        => 0,
    'mRNA'        => 0,
    '5pCDS3pmRNA' => 0,
    '5pCDSmRNA'   => 0,
    'CDS3pmRNA'   => 0,
    'CDSmRNA'     => 0,
    'ncRNA'       => 0,
    'byGroup'     => {}
);


foreach my $chr ( keys %allgene ) {


    my ($group) = $chr=~/^(Group(?:\d+|Un\d*)|MT)/i;

    @{ $count{'byGroup'}->{$group} }{ 'scaffold', 'gene', 'mRNA', '5pCDS3pmRNA', '5pCDSmRNA', 'CDS3pmRNA', 'CDSmRNA', 'ncRNA' } = (0, 0, 0, 0, 0, 0, 0, 0) unless (exists $count{'byGroup'}->{$group});
    
    $count{'scaffold'}++;
    $count{'byGroup'}->{$group}->{'scaffold'}++;

    my $in = Bio::SeqIO->new(
        -file   => "$genomedir/" . $chr . '.fa',
        -format => 'fasta'
    );
    
    my $seq = $in->next_seq();
    
    foreach my $id ( keys %{$allgene{ $chr } } ) {
        # gene
        if (!exists $parent{$id}) {
#            print $id,"\t",$feature{$id}->{'primary_tag'},"\n";
            $count{'gene'}++;
            $count{'byGroup'}->{$group}->{'gene'}++;
            
#            print STDERR "GENE: $id","\n";

            foreach my $mRNAfeat (@{ $feature{$id}->{'child'}->{ 'mRNA' } } ) {
                my ($mRNAid) = $mRNAfeat->get_tag_values('ID');
                
#                print STDERR "\t",$mRNAid,"\n";

                    
                my $strseq='';
                {
                    my @fe;
                    if ($mRNAfeat->strand() == -1) {
                        @fe = sort { $b->start() <=> $a->start() } @{ $feature{$mRNAid}->{'child'}->{'exon'} };
                    }
                    elsif ($mRNAfeat->strand() == 1) {
                        @fe = sort { $a->start() <=> $b->start() } @{ $feature{$mRNAid}->{'child'}->{'exon'} };
                    }
                    else {
                        $LOGGER->logdie("Cannot found strand for $mRNAid");
                    }
                    foreach my $f ( @fe ) {
                        #print "\t",$f->primary_tag(),"\n";
                        #print "\t",$seq->subseq($f->location()),"\n";
                        $strseq.=$seq->subseq($f->location());                    
                    }
                }
                
                my $extseq='';
                if ($threepminsize>0) {
                    my ($ts, $te) = ($mRNAfeat->start(), $mRNAfeat->end());
                    if ($mRNAfeat->strand() == -1) {
                        if ($ts > $threepminsize) {
                            $extseq = $seq->trunc($ts-1-($threepminsize-1), $ts-1)->revcom()->seq();
                        }
                        else {
                            if ($ts>1) {
                                $extseq = $seq->trunc(1, $ts-1)->revcom()->seq();
                            }
                        }
                    }
                    elsif ($mRNAfeat->strand() == 1) {
                        if ( ($te+1+$threepminsize) <= $seq->length() ) {
                            $extseq = $seq->trunc($te+1, $te+1+($threepminsize-1))->seq();
                        }
                        else {
                            if ( $te < $seq->length() ) {
                                $extseq = $seq->trunc($te+1, $seq->length())->seq();
                            }
                        }
                    }
                    else {
                        $LOGGER->logdie("Cannot found strand for $mRNAid");
                    }
                }
                
                my $smallest3plen = $threepminsize;
                if (exists $feature{$mRNAid}->{'child'}->{'CDS'}) {

                    my %eachCDS;
                    {
                        my @allcds = @{ $feature{$mRNAid}->{'child'}->{'CDS'} };
                        foreach my $cf (@allcds) {
                            my ($CDSid);
                            if ($cf->has_tag('ID')) {
                                ($CDSid) = $cf->get_tag_values('ID');
                                $CDSid=~s/-\d+$//;
                            }
                            else {
                                # OGS (each mRNA has only one CDS)
                                ($CDSid) =  $cf->get_tag_values('Parent');
                                $CDSid=~s/-R([ABC])/-P$1/;
                            }
                            my ($CDSgroup) = $CDSid=~/^([^#]+)/;
                            push(@{$eachCDS{ $CDSgroup }}, $cf);
                        }
                    }
                    foreach my $CDSgroup (keys %eachCDS) {
                        
                        my @fc;
                        if ($mRNAfeat->strand() == -1) {
                            @fc = sort { $b->start() <=> $a->start() } @{ $eachCDS{ $CDSgroup } };
                        }
                        elsif ($mRNAfeat->strand() == 1) {
                            @fc = sort { $a->start() <=> $b->start() } @{ $eachCDS{ $CDSgroup } };
                        }
                        else {
                            $LOGGER->logdie("Cannot found strand for $mRNAid");
                        }
                        my $cdsseq='';
                        foreach my $f ( @fc ) {
                            $cdsseq.=$seq->subseq($f->location());                    
                        }

                        {
                            my $seqobj = Bio::Seq->new(
                                -seq            => $cdsseq,
                                -display_id     => 'CDSmRNA'.'|'.$CDSgroup.'|'.$mRNAid.'|'.$chr,
                                -description    => ''
                            );

                            unless (exists $out{'CDSmRNA'}) {
                                $out{'CDSmRNA'} = Bio::SeqIO->new(
                                    -file   => ">$outdir/$outprefix".'CDSmRNA'.'.fa',
                                    -format => 'fasta'
                                );
                            }
                            $out{'CDSmRNA'}->write_seq($seqobj);
                        }

                        if ($strseq=~/(.*)$cdsseq(.*)/) {
                            my ($fivep,$threep) = ($1,$2);
                            if ($fivep) {

                                my $seqobj = Bio::Seq->new(
                                    -seq            => $fivep,
                                    -display_id     => 'UTR5p'.'|'.$CDSgroup.'|'.$mRNAid.'|'.$chr,
                                    -description    => ''
                                );
                                unless (exists $out{'UTR5p'}) {
                                    $out{'UTR5p'} = Bio::SeqIO->new(
                                        -file   => ">$outdir/$outprefix".'UTR5p'.'.fa',
                                        -format => 'fasta'
                                    );
                                }
                                $out{'UTR5p'}->write_seq($seqobj);

                            }

                            $smallest3plen = length($threep);

                            my $ext = 0;
                            if (length($extseq)>0) {
                                my $lenthreep = (($threep) ? length($threep) : 0);
                                if (length($extseq)>$lenthreep) {
                                    $ext = length($extseq)-$lenthreep;
                                } 
                                
                                $threep.=substr($extseq,0, $ext);
                            }

                            if ($threep) {

                                my $seqobj = Bio::Seq->new(
                                    -seq            => $threep,
                                    -display_id     => 'UTR3p'.(($ext) ? '-'.$ext : '').'|'.$CDSgroup.'|'.$mRNAid.'|'.$chr,
                                    -description    => ''
                                );
                                unless (exists $out{'UTR3p'}) {
                                    $out{'UTR3p'} = Bio::SeqIO->new(
                                        -file   => ">$outdir/$outprefix".'UTR3p'.'.fa',
                                        -format => 'fasta'
                                    );
                                }
                                $out{'UTR3p'}->write_seq($seqobj);

                            }
                            

                            if ($fivep) {
                                if ($threep) {
                                    $count{'5pCDS3pmRNA'}++;
                                    $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}++;
                                }
                                else {
                                    $count{'5pCDSmRNA'}++;
                                    $count{'byGroup'}->{$group}->{'5pCDSmRNA'}++;
                                }
                            }
                            elsif ($threep) {
                                if ($fivep) {
                                    $count{'5pCDS3pmRNA'}++;
                                    $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}++;
                                }
                                else {
                                    $count{'CDS3pmRNA'}++;
                                    $count{'byGroup'}->{$group}->{'CDS3pmRNA'}++;
                                }
                                    
                            }
                            else {
                                $count{'CDSmRNA'}++;
                                $count{'byGroup'}->{$group}->{'CDSmRNA'}++;
                            }
                            #print $CDSgroup."\t".$fivep."\n";
                            #print $CDSgroup."\t".$threep."\n";
                        }
                        else {
                            $LOGGER->logdie("Something wrong! Not match CDS ($CDSgroup) into mRNA ($mRNAid)");
                        }
                    }
                    
                    {
                         my $ext = 0;
                         my $threepcomp = '';
                         if (length($extseq)>0) {
                                my $lenthreep = $smallest3plen;
                                if (length($extseq)>$lenthreep) {
                                    $ext = length($extseq)-$lenthreep;
                                } 
                                
                                $threepcomp.=substr($extseq,0, $ext);
                        }
                        $count{'mRNA'}++;
                        $count{'byGroup'}->{$group}->{'mRNA'}++;
                        
                        my $seqobj = Bio::Seq->new(
                            -seq            => $strseq.$threepcomp,
                            -display_id     => 'mRNA'.(($ext) ? '-'.$ext : '').'|'.$mRNAid.'|'.$parent{$mRNAid}.'|'.$chr,
                            -description    => ''
                        );
                        unless (exists $out{'mRNA'}) {
                            $out{'mRNA'} = Bio::SeqIO->new(
                                    -file   => ">$outdir/$outprefix".'mRNA'.'.fa',
                                    -format => 'fasta'
                                );
                        }
                        $out{'mRNA'}->write_seq($seqobj);
                    }
                }
                else {
                    $count{'ncRNA'}++;
                    $count{'byGroup'}->{$group}->{'ncRNA'}++;
                    
                    my $seqobj = Bio::Seq->new(
                        -seq            => $strseq,
                        -display_id     => 'ncRNA'.'|'.$mRNAid.'|'.$parent{$mRNAid}.'|'.$chr,
                        -description    => ''
                    );
                    unless (exists $out{'ncRNA'}) {
                        $out{'ncRNA'} = Bio::SeqIO->new(
                                -file   => ">$outdir/$outprefix".'ncRNA'.'.fa',
                                -format => 'fasta'
                            );
                    }
                    $out{'ncRNA'}->write_seq($seqobj);
                }                    
            }

        }
    }
}

if ($statsfile) {
    open(STATS, ">", $statsfile) or $LOGGER->logdie("$!");

    print STATS "
*** GENERAL REPORT ***

Nº of scaffolds scanned.: $count{'scaffold'}
Nº of genes found.......: $count{'gene'}
Nº of mRNAs found.......: $count{'mRNA'}
    5p CDS 3p...........: $count{'5pCDS3pmRNA'}
    5p CDS..............: $count{'5pCDSmRNA'}
    CDS 3p..............: $count{'CDS3pmRNA'}
    CDS.................: $count{'CDSmRNA'}
Nº of ncRNAs found......: $count{'ncRNA'}
";

    print STATS "
*** REPORT BY GROUP ***
";

    foreach my $group (sort {$a cmp $b} keys %{$count{'byGroup'}} ) {
        print STATS "
* $group

    Nº of scaffolds scanned.: $count{'byGroup'}->{$group}->{'scaffold'}
    Nº of genes found.......: $count{'byGroup'}->{$group}->{'gene'}
    Nº of mRNAs found.......: $count{'byGroup'}->{$group}->{'mRNA'}
        5p CDS 3p...........: $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}
        5p CDS..............: $count{'byGroup'}->{$group}->{'5pCDSmRNA'}
        CDS 3p..............: $count{'byGroup'}->{$group}->{'CDS3pmRNA'}
        CDS.................: $count{'byGroup'}->{$group}->{'CDSmRNA'}
    Nº of ncRNAs found......: $count{'byGroup'}->{$group}->{'ncRNA'}
";
    }

    close(STATS);
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

        -h      --help  	Help
        -l      --level 	Log level [Default: FATAL]
        -g	--genomefile    Genome file
        -d	--genomedir	    Genome directory
        -f	--gff		    GFF file
        -o	--outdir	    Output directory
        -p	--outprefix	    Output filenames' prefix
        -s  --statsfile     Output Statistics file
        -m  --threepminsize 3' minimum size (extends if possible) [Default: 0]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

