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

use Bio::SeqIO;

my $in = Bio::SeqIO->new(-file=>$infile);

while(my $seq=$in->next_seq()) {
    my $transcript_id = $seq->display_id() . '.' . $seq->version();

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
                    print ">pep|$protein_id|$geneName_name|".$seq->display_id()."\n";
                    my $seq='';
                    foreach my $s ( $feat->each_tag_value('translation') ) {
                          $seq.=$s;
                    }
                    $seq=~s/(.{50})/$1\n/g;
                    print $seq,"\n";

                    $LOGGER->logdie('Not found /protein_id for '.$transcript_id) unless ($protein_id);
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

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

