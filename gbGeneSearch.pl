#!/usr/bin/env perl
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

my ($level, $gbfile, $idlist, $extra);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "g|gbfile=s"=>\$gbfile,
            "i|id=s"=>\$idlist,
            "e|extra=i"=>\$extra
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

$extra||=0;

$idlist=~s/ //g;
my @idlist = split(/,/, $idlist);

$LOGGER->logdie("Missing GeneID list") unless (scalar(@idlist)>0);

my %id;
@id{@idlist} = ();

use Bio::SeqIO;

my $seqin = Bio::SeqIO->new(-file => $gbfile, -format => "genbank" );
my $seq_object = $seqin->next_seq();

my $seqout = Bio::SeqIO->new(   -fh=>\*STDOUT,
                                -format=>'FASTA');

for my $feat_object ($seq_object->get_SeqFeatures) {
    if ($feat_object->primary_tag() eq 'gene') {


        if ($feat_object->has_tag("db_xref")) {
            for my $value ($feat_object->get_tag_values("db_xref")) {
                if ($value =~/GeneID:(\d+)/) {
                    my $GeneID=$1;
                    if (exists $id{$GeneID}) {
                        my $GeneName='';
                        if ( $feat_object->has_tag("gene") ) {
                            for my $value ( $feat_object->get_tag_values("gene") ) {
                                if ( $value =~ /(\S+)/ ) {
                                    $GeneName=$value;
                                }
                            }
                        } elsif ( $feat_object->has_tag("locus_tag") ) {
                            for my $value ( $feat_object->get_tag_values("locus_tag") ) {
                                if ( $value =~ /(\S+)/ ) {
                                    $GeneName=$value;
                                }
                            }
                        }
                        
                        #print $feat_object->location->start(),"\t",$feat_object->location->end(),"\t",$feat_object->location->strand(),"\n";
                        my $start = $feat_object->location->start();
                        my $end =   $feat_object->location->end();
                        my $chrlen = $feat_object->entire_seq->length();
                        $start = ( (($start-$extra)<1) ? 1 : $start-$extra );
                        $end = ( (($end+$extra)>$chrlen) ? $chrlen : $end+$extra );

                        my $seq = Bio::Seq->new(  -display_id => $seq_object->display_id().':'.(($GeneName ne '') ? $GeneName.':' : '').$GeneID,
                                                  -seq =>$feat_object->entire_seq->subseq($start, $end));
                        
                        $seqout->write_seq($seq);

                        delete($id{$GeneID});
                        if (scalar(keys %id) == 0) {                        
                            last;
                        }
                    }
                }                            
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
        -g      --gbfile    Input GenBank (gbwithparts) file 
        -i      --id        GeneID list
        -e      --extra     Extra bases

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

