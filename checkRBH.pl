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

use File::Basename;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $indir);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|indir=s"=> \$indir
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

$LOGGER->logdie("Missing input directory") unless ($indir);
$LOGGER->logdie("Wrong input directory ($indir)") unless (-e $indir);
    
my %group;

foreach my $f ( glob("$indir/MAP_*.txt") ) {
    next if ($f=~/_RBH\.txt/);

    my $bn=basename($f,'.txt');
    $bn=~s/^MAP_//;
    my ($g1, $g2) = split(/_x_/, $bn);
    #print "\t$bn\t$g1\t$g2\n";
    $group{$g1}->{$g2} = $f;
}

foreach my $k1 (keys %group) {
    foreach my $k2 (keys %{ $group{$k1} }) {
        print "Processing $k1 x $k2 ...\n";
        
        my %rna;
        my %pep;

        open(TWO, "<", $group{$k2}->{$k1}) or $LOGGER->logdie($!);
        while(<TWO>) {
            chomp;
            my ($refrna, $refpep, $relatedrna, $relatedpep) = split(/\t/, $_);
            foreach my $r ( split(/;/, $relatedrna) ) {
                    $r=~s/\[.+$//;
                    $rna{ $refrna }->{ $r } = undef;
            }
            foreach my $p ( split(/;/, $relatedpep) ) {
                    $p=~s/\[.+$//;
                    $pep{ $refpep }->{ $p } = undef;
            }
        }
        close(TWO);
        
        my $bn=basename($group{$k1}->{$k2}, '.txt');

        open(OUT, ">", $indir.'/'.$bn.'_RBH.txt') or $LOGGER->logdie($!);
        
        open(ONE, "<", $group{$k1}->{$k2}) or $LOGGER->logdie($!);
        while(<ONE>) {
            chomp;

            my ($refrna, $refpep, $relatedrna, $relatedpep) = split(/\t/, $_);
            my @rs;
            foreach my $r ( split(/;/, $relatedrna) ) {
                    my $t=$r;
                    $t=~s/\[.+$//;
                    if (exists $rna{ $t }->{ $refrna }) {
                        $r.='*';
                    }
                    push(@rs, $r);
            }
            $relatedrna=join(';', @rs);
            my @ps;
            foreach my $p ( split(/;/, $relatedpep) ) {
                    my $t=$p;
                    $t=~s/\[.+$//;
                    if (exists $pep{ $t }->{ $refpep }) {
                        $p.='*';
                    }
                    push(@ps, $p);
            }
            $relatedpep=join(';', @ps);

            print OUT join("\t", $refrna, $refpep, $relatedrna, $relatedpep),"\n";
        }
        
        close(ONE);
        
        close(OUT);
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

        -h      --help  Help
        -l      --level Log level [Default: FATAL]
        -i      --indir Input directory

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

