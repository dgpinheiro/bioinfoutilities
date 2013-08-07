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

my ($level, $mirfile, $rnafoldfile, $check);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "m|mirfile=s"=>\$mirfile,
            "r|rnafoldfile=s"=>\$rnafoldfile,
            "c|check"=>\$check
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

$LOGGER->logdie("Missing mir file") unless ($mirfile);
$LOGGER->logdie("Wrong mir file ($mirfile)") unless (-e $mirfile);

$LOGGER->logdie("Missing RNAfold file") unless ($rnafoldfile);
$LOGGER->logdie("Wrong RNAfold file ($rnafoldfile)") unless (-e $rnafoldfile);

use Bio::SeqIO;

my %mature;
my $seqin = Bio::SeqIO->new(-file=>$mirfile, -format=>'FASTA');

while (my $seq = $seqin->next_seq() ) {
    $mature{$seq->display_id()} = $seq->seq();
}
                            
my %RESULT;

open(RFF, "<", $rnafoldfile) or $LOGGER->logdie($!);
my $last_id = '';
my $last_seq = '';
while(<RFF>) {
    chomp;
    if ( ($.%3)==1 ) {
        if ($_=~/^>(\S+)/) {
            $last_id=$1;
        }
        else {
            $LOGGER->logdie("Wrong pattern in RNAfold file. Not found the sequence id ($_).");
        }        
    }
    elsif ( ($.%3)==2 ) {
        if ($_=~/^[ACGU]+$/) {
            $last_seq = $_;
        }
        else {
            $LOGGER->logdie("Wrong pattern in RNAfold file. Not found the RNA sequence ($_).");
        }        
    }
    else {
        if ($_=~/^([\.\(\)]+)\s+/) {
            my $fold = $1;
            my $found = 0;
            my @f;
            foreach my $m (keys %mature) {
                my $pattern = $mature{$m};
                if ($last_seq=~/$pattern/) {
#                    print $last_id,"\t",$m,"\n";
                    $found++;
                    push(@f, $m);
                }
            }
            if ($found != 1) {
                $LOGGER->logwarn("Found $found for $last_id (".join(';', @f).")");
            }

            foreach my $m (@f) {
                my $pattern = $mature{$m};
                while( $last_seq=~m/$pattern/g ) {
                    my $end = pos($last_seq)-1;
                    my $start = $end-length($pattern)+1;
#                    print $m,"\t",$mature{$m},"\n",$last_seq,"\n",$start,'-',$end,"\n";
                    my @splittedseq = split(//, $last_seq);
                    my @splittedfold = split(//, $fold);
#                    print join(";", @splittedseq[$start..$end]),"\n";
#                    print join(";", @splittedfold[$start..$end]),"\n";
                    my $leftb = 0;
                    my $rightb = 0;
                    foreach my $b (@splittedfold[$start..$end]) {
                        if ($b eq '(') {
                            $leftb++;
                        }
                        elsif ($b eq ')') {
                            $rightb++;
                        }
                    }
                    my $result;
                    my $star_result;
                    my $star_sequence;
                    if ($leftb == $rightb) {
                        $LOGGER->logdie("Wrong result for $m and $mature{$m}");
                    }
                    else {
                        if ($leftb > $rightb) {
                            $result = '-5p';
                        }
                        else {
                            $result = '-3p';
                        }
                    }
#                    print "RESULT: $result\n";
#                    print $m,"\t",$mature{$m},"\n",$last_seq,"\n",$start,'-',$end,"\n",$result,"\n";
                    if ($check) {
                        if ($m=~/[35]p$/) {
                            $LOGGER->logdie("MIR: $m RESULT: $result") if ($m!~/$result/);
                        }                            
                    }
                    else {
                        if ($m !~ /(-[35]p)$/) {
                            $RESULT{$m.$result}->{$mature{$m}} = undef;
                        }
                        else {
                            my $lastres = $1;
                            $LOGGER->logdie("Wrong result: $lastres ne $result") unless ($lastres eq $result);
                            $RESULT{$m}->{$mature{$m}} = undef;
                        }
                    }                    
                }
            }
        }
        else {
            $LOGGER->logdie("Wrong pattern in RNAfold file. Not found the RNAfold result ($_).");
        }        
        
    }
}
close(RFF);

foreach my $id (sort {$a cmp $b} keys %RESULT) {
    my @r = keys %{ $RESULT{$id} };
    $LOGGER->logdie("This miRNA $id has more than one sequence: ".join(';', @r)) if (scalar(@r)>1);
    print '>'.$id."\n".$r[0]."\n";
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

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -m      --mirfile       mirbase mature miRNA
        -r      --rnafoldfile   RNAfold result with mirbase hairpin
        -c      --check         Checking mode

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

