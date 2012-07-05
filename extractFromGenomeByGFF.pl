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

my ($level, $gfffile, $genomedir, $outdir, $outprefix);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "g|genome=s"=>\$genomedir,
            "f|gff=s"=>\$gfffile,
            "o|outdir=s"=>\$outdir,
            "p|outprefix=s"=>\$outprefix
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

$LOGGER->logdie('Missing genome directory') unless ($genomedir);
$LOGGER->logdie('Wrong genome directory ('.$genomedir.')') unless (-d $genomedir);

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
while ($feat = $gffio->next_feature()) {
    $c++;
    die "Not found ".$feat->seq_id()." on $genomedir"  unless (-e "$genomedir/".$feat->seq_id().'.fa');

    my ($id,$parent);
    if ($feat->has_tag('ID')) {
        ($id) = $feat->get_tag_values('ID');
        unless (exists $feature{$id}) {
            $feature{$id} = {   'primary_tag'   => $feat->primary_tag(),
                                'child'         => undef };
            push(@{$chromosome{ $feat->seq_id() }}, $id);
        }                            
    }

    if ($feat->has_tag('Parent')) {
        ($parent) = $feat->get_tag_values('Parent');
        if ($id) {
            $parent{$id} = $parent unless (exists $parent{$id});
        }    

        push(@{$feature{$parent}->{'child'}->{ $feat->primary_tag() }}, $feat);
    }
    
    #last if ($c >= 26);
}

my %out;
my %count = (
    'scaffold'    => 0,
    'gene'        => 0,
    'mRNA'        => 0,
    '5pCDS3pmRNA' => {'all'=>0, 'small'=>0},
    '5pCDSmRNA'   => {'all'=>0, 'small'=>0}, 
    'CDS3pmRNA'   => {'all'=>0, 'small'=>0},
    'CDSmRNA'     => {'all'=>0, 'small'=>0},
    'ncRNA'       => 0,
    'byGroup'     => {}
);


foreach my $chr ( keys %chromosome ) {


    my ($group) = $chr=~/^(Group(?:\d+|Un))/;

    @{ $count{'byGroup'}->{$group} }{ 'scaffold', 'gene', 'mRNA', '5pCDS3pmRNA', '5pCDSmRNA', 'CDS3pmRNA', 'CDSmRNA', 'ncRNA' } = (0, 0, 0, {'all'=>0, 'small'=>0}, {'all'=>0, 'small'=>0}, {'all'=>0, 'small'=>0}, {'all'=>0, 'small'=>0}, 0) unless (exists $count{'byGroup'}->{$group});
    
    $count{'scaffold'}++;
    $count{'byGroup'}->{$group}->{'scaffold'}++;

    my $in = Bio::SeqIO->new(
        -file   => "$genomedir/" . $chr . '.fa',
        -format => 'fasta'
    );
    
    my $seq = $in->next_seq();

    foreach my $id ( @{ $chromosome{$chr} } ) {
        # gene
        if (!exists $parent{$id}) {
#            print $id,"\t",$feature{$id}->{'primary_tag'},"\n";
            $count{'gene'}++;
            $count{'byGroup'}->{$group}->{'gene'}++;
        }
        # mRNA
        if (exists $parent{$id}) {
            $count{'mRNA'}++;
            $count{'byGroup'}->{$group}->{'mRNA'}++;
            
            my %seqlength;
            foreach my $k (keys %{ $feature{$id}->{'child'} }) {
                my $strseq='';
                foreach my $f ( @{ $feature{$id}->{'child'}->{$k} } ) {
#                    print "\t",$f->primary_tag(),"\n";
#                    print "\t",$seq->subseq($f->location()),"\n";
                    $strseq.=$seq->subseq($f->location());                    
                }
                my $seqobj = Bio::Seq->new(
                    -seq            => $strseq,
                    -display_id     => $k.'|'.$id.'|'.$parent{$id}.'|'.$chr,
                    -description    => ''
                );
#                print '>',$seqobj->display_id(),"\n";
#                print $seqobj->seq(),"\n";
                unless (exists $out{$k}) {
                    $out{$k} = Bio::SeqIO->new(
                            -file   => ">$outdir/$outprefix".$k.'.fa',
                            -format => 'fasta'
                        );
                }
                $seqlength{$k} = length($strseq);
                $out{$k}->write_seq($seqobj);
            }



            if ( exists $feature{$id}->{'child'}->{'CDS'} ) {
                if ( exists $feature{$id}->{'child'}->{'five_prime_UTR'} ) {
                    if ( exists $feature{$id}->{'child'}->{'three_prime_UTR'} ) {
                        $count{'5pCDS3pmRNA'}->{'all'}++;
                        $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}->{'all'}++;
                        if ($seqlength{'CDS'} < 60) {
                            $count{'5pCDS3pmRNA'}->{'small'}++;
                            $count{'byGroup'}->{$group}->{'5pCDS3pmRNA'}->{'small'}++;
                        } 
                    }
                    else {
                        $count{'5pCDSmRNA'}->{'all'}++;
                        $count{'byGroup'}->{$group}->{'5pCDSmRNA'}->{'all'}++;
                        if ($seqlength{'CDS'} < 60) {
                            $count{'5pCDSmRNA'}->{'small'}++;
                            $count{'byGroup'}->{$group}->{'5pCDSmRNA'}->{'small'}++;
                        } 
                    }                        
                }
                else {
                    if ( exists $feature{$id}->{'child'}->{'three_prime_UTR'} ) {
                        $count{'CDS3pmRNA'}->{'all'}++;
                        $count{'byGroup'}->{$group}->{'CDS3pmRNA'}->{'all'}++;
                        if ($seqlength{'CDS'} < 60) {
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
            else {
                $count{'ncRNA'}++;
                $count{'byGroup'}->{$group}->{'ncRNA'}++;
            }
#            print $id,"\t",$parent{$id},"\t",$feature{$id}->{'primary_tag'},"\n";
            
        }
    }
}

print STDERR "
*** GENERAL REPORT ***

Nº of scaffolds scanned.: $count{'scaffold'}
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

    Nº of scaffolds scanned.: $count{'byGroup'}->{$group}->{'scaffold'}
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

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help  	Help
        -l      --level 	Log level [Default: FATAL]
        -g	--genome	Genome directory
        -f	--gff		GFF file
        -o	--outdir	Output directory
        -p	--outprefix	Output filenames' prefix

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

