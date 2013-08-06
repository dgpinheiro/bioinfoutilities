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
    
    selSEQfromFASTA.pl - select sequences from multi fasta
    
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
        -f/--fasta  FASTA file            
        -i/--id     Id list/file
        
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
    Log::Log4perl->easy_init($WARN);
    $LOGGER = Log::Log4perl->get_logger($0);
}

use Bio::SeqIO;


my ($level, $infile, $idparam);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|id=s"=>\$idparam,
            "f|fasta=s"=>\$infile
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

my $seqin;

if ($infile) {
    $LOGGER->logdie("Wrong multi fasta file ($infile)") unless (-e $infile);

    $seqin = Bio::SeqIO->new(    -file=>$infile,
                                 -format=>'FASTA');
} else {
    $seqin = Bio::SeqIO->new(    -fh=>\*STDIN,
                                 -format=>'FASTA');
}

my %SELECTION;
foreach my $each_idparam (split(/\,/, $idparam)) {
    if (-e $each_idparam) {
        open(IN, "<", $each_idparam) or $LOGGER->logdie($!);
        while(<IN>) {
            chomp;
            $SELECTION{$_}=0;
        }
        close(IN);
    } else {
        $SELECTION{$each_idparam}=0;
    }
}



my $seqout = Bio::SeqIO->new(   -fh=>\*STDOUT,
                                -format=>'FASTA');

while(my $seq=$seqin->next_seq() ) {
    foreach my $ID (keys %SELECTION ) {
        my $qID = quotemeta($ID);
        if ($seq->display_id() =~ /$qID/) {
            $seqout->write_seq( $seq );
            $SELECTION{$ID}++;
        }
    }
}

foreach my $ID (keys %SELECTION) {
    unless ($SELECTION{$ID}) {
        $LOGGER->logwarn("Not found ID: $ID");
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
        -f      --fasta     Input fasta file or <STDIN>
        -i      --id        ID list/file 

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

