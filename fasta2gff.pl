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

use POSIX 'isatty';
use FileHandle;
use Bio::SeqIO;
use Bio::FeatureIO;
use Bio::SeqFeature::Annotated;

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

my $fh;

if ($infile) {

    $LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

    $fh = FileHandle->new;
    $fh->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fh = \*STDIN;
    } else {
        $LOGGER->logdie("Missing input file (-i/--infile) or STDIN data");
    }
}

my $seqin = Bio::SeqIO->new(-fh=>$fh, -format=>"FASTA");

my $featureOut = Bio::FeatureIO->new(
    -format  => 'gff',
    -version => 3,
    -fh      => \*STDOUT,
    -validate_terms =>
      1,    #boolean. validate ontology terms online?  default 0 (false).
);

while(my $seq=$seqin->next_seq()) {
    my $feature_gene = Bio::SeqFeature::Annotated->new( -seq_id=>$seq->display_id(),
                                                        -name=>$seq->display_id(),
                                                        -source=>"DeNovoAssembly",
                                                        -type=>"gene",
                                                        -length=>$seq->length(),
                                                        -start=>1,
                                                        -end=>$seq->length(),
                                                        -strand=>'+'
                                                 );

    $feature_gene->add_tag_value('ID',$seq->display_id());

    $featureOut->write_feature( $feature_gene );

    my $feature_exon = Bio::SeqFeature::Annotated->new( -seq_id=>$seq->display_id(),
                                                        -name=>$seq->display_id(),
                                                        -source=>"DeNovoAssembly",
                                                        -type=>"exon",
                                                        -length=>$seq->length(),
                                                        -start=>1,
                                                        -end=>$seq->length(),
                                                        -strand=>'+'
                                                     );
    $feature_exon->add_tag_value('ID','exon_'.$seq->display_id());
    $feature_exon->add_tag_value('Parent',$seq->display_id());

    $featureOut->write_feature( $feature_exon );
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
        -i      --infile    Input fasta file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

