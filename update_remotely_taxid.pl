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
#  Copyright (C) 2019  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
#  Faculdade de Ciências Agrárias e Veterinárias (FCAV)
#  Laboratório de Bioinformática (LB)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://www.fcav.unesp.br
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

Copyright (C) 2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Readonly;
use Getopt::Long;

use constant MERGED_TAB_FILE_PATH=>"/usr/local/bioinfo/Krona/KronaTools/taxonomy/merged.tab";


use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile);

#Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile
    ) or &Usage();

use FileHandle;    
use POSIX 'isatty';
use Storable;

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
        $LOGGER->logdie("Missing input file or STDIN data");
    }
}

my %data;

$|=1;
while(<$fh>) {
    chomp;
    my $taxids=$_;
    $taxids=~s/^\s+//;
    $taxids=~s/^\s+$//;
    my @taxids=split(/\s+/, $taxids);
    foreach my $taxid (@taxids) {
        $taxid+=0;
        unless (exists $data{$taxid}) {
            my $cmd1="echo $taxid | GetTaxInfo.pl | cut -f 5 | tail -1";
            my $res1=`$cmd1`;
            chomp($res1);
            if (($res1=~/^\d+$/)&&($res1 == $taxid)) {
                my $cmd2='esearch -db taxonomy -query "'.$taxid.'[uid]" | efetch -db taxonomy -format docsum';
                sleep(1);
                my $res2=`$cmd2`;
                my $taxname=undef;
                my $newtaxid=undef;
                my $status;
                foreach my $l (split(/\n/, $res2)) {
                    unless ($status) {
                        if ($l=~/<Status>([^<]+)<\/Status>/) {
                            $status=$1;
                        }
                    } else {
                        if ($status eq 'active') {
                            if ($l=~/<TaxId>(\d+)<\/TaxId>/) {
                                $newtaxid=$1;
                                if ($newtaxid) {
        #                            $newtaxid=`echo "$taxname" | taxonkit name2taxid | cut -f 2`;
                                    chomp($newtaxid);
                                    last;
                                } else {
                                    $LOGGER->logdie("Not found esearch/efetch ($status) for $taxid.");
                                }
                            }
                        } elsif ($status eq 'merged') {                        
                            if ($l=~/<AkaTaxId>(\d+)<\/AkaTaxId>/) {
                                $newtaxid=$1;
                                if ($newtaxid) {
                                    chomp($newtaxid);
                                    last;
                                } else {
                                    $LOGGER->logdie("Not found esearch/efetch ($status) for $taxid.");
                                }
                            }
                        } elsif ($status eq 'deleted') {
                            if ($l=~/<ScientificName>([^<]+)<\/ScientificName>/) {
                                $taxname=$1;
                                if ($taxname) {
                                    $newtaxid=`echo "$taxname" | taxonkit name2taxid | cut -f 2`;
                                    chomp($newtaxid);
                                    last;
                                } else {
                                    $LOGGER->logdie("Not found esearch/efetch ($status) for $taxid.");
                                }
                            }

                        } else {
                            $LOGGER->logdie("TaxId $taxid Status $status");
                        }
                    }                        
                }
                chomp($newtaxid);
                $data{$taxid} = $newtaxid;
            } else {
                chomp($taxid);
                $data{$taxid} = $taxid;
            }
        }
        print $taxid,"\t",$data{$taxid},"\n";
    }
}

$fh->close();

# Subroutines

sub Usage {
    my ($msg) = @_;

	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

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

