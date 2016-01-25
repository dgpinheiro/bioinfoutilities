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

my ($level, $infile, $colname, $normby, $filter);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "c|colname=s"=>\$colname,
            "n|norm=i"=>\$normby,
            "f|filter=i"=>\$filter
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


$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

$LOGGER->logdie("Missing column name with taxonomy in mpa format") unless ($colname);


open(IN, "<", $infile) or $LOGGER->logdie($!);
my $header_line=<IN>;
chomp($header_line);
my @header=split(/\t/, $header_line);

{
    my %check;
    @check{ @header } = (1..$#header);

    $LOGGER->logdie("Not found column $colname in input file") unless (exists $check{ $colname });
}

my %result;

my @hnumber;
my %total;

while(<IN>) {
    chomp;
    my %data;
    @data{ @header } = split(/\t/, $_);
    if (scalar(@hnumber) == 0) {
        foreach my $h (@header) {
            if ( $data{ $h } =~ /^\d+(?:\.\d+)?$/ ) {
                push(@hnumber, $h);
                $total{ $h } = 0;
            }
        }
    }
    if ($data{ $colname }) {
        my @taxlevel = (split(/\|/, $data{ $colname }));
        for (my $i=0; $i<= $#taxlevel; $i++) {
            my @v;
            for (my $j=0; $j<= $i; $j++) {
                push(@v, $taxlevel[$j]);
            }

            my $class=join('|', @v);
            
            $result{ $class }->{ 'level' } = scalar(@v) unless (exists $result{ $class }->{ 'level' });
            
            foreach my $h (@hnumber) {
                $result{ $class }->{ 'all' } = 0 unless (exists $result{ 'unclassified' }->{ 'all' });
                $result{ $class }->{ 'all' }+=$data{ $h };
                $result{ $class }->{ $h } = 0 unless (exists $result{ $class }->{ $h });
                $result{ $class }->{ $h }+= $data{ $h };
                $total{ $h }+=$data{ $h };
            }
        }    
        
    } else {
        $result{ 'unclassified' }->{ 'level' } = 0 unless (exists $result{ 'unclassified' }->{ 'level' });

        foreach my $h (@hnumber) {
            $result{ 'unclassified' }->{ 'all' } = 0 unless (exists $result{ 'unclassified' }->{ 'all' });
            $result{ 'unclassified' }->{ 'all' }+=$data{ $h };

            $result{ 'unclassified' }->{ $h } = 0 unless (exists $result{ 'unclassified' }->{ $h });
            $result{ 'unclassified' }->{ $h }+=$data{ $h };

            $total{ $h }+=$data{ $h };
        }
    }
}
close(IN);


my %normfactor;
if ($normby) {
    foreach my $h (@hnumber) {
        $normfactor{ $h } = $normby/$total{$h};
    }
} else {
    foreach my $h (@hnumber) {
        $normfactor{ $h } = 1;
    }
}

print join("\t", $colname, @hnumber),"\n";
foreach my $class (sort {$result{$b}->{'all'} <=> $result{$a}->{'all'}} keys %result) {
    if  ((!$filter)||($result{$class}->{'level'} == $filter)) { 
        my @v;
        foreach my $h (@hnumber) {
            push(@v, $result{$class}->{$h} * $normfactor{ $h });
        }
        print $class,"\t",join("\t", map { sprintf("%.2f",$_) } @v),"\n";
    }
}

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file (count matrix with taxonomy in mpa format)
        -c      --colname   Annotation column name of taxonomy in mpa format
        -n      --normby    Normalization by [Default: Off]
        -f      --filter    Filter by level [Default: Off]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

