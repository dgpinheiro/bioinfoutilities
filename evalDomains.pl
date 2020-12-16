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

use constant THRESHOLD_DEFAULT=>1;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, @field, $threshold);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "f|field=s"=>\@field,
            "t|threshold=i"=>\$threshold
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


$threshold||=THRESHOLD_DEFAULT;

$field[0]=~s/\s+//;

$LOGGER->logdie("Missing domain fields") unless ($field[0]);

@field=split(/,/, $field[0]);

$LOGGER->logdie("Threshold must be lower or equal than number of domain fields") unless ($threshold<=scalar(@field));

#print "\n\n",join(";", @field),"\n\n";

open(IN, "<", $infile) or $LOGGER->logdie($!);

while(<IN>) {
    chomp;
    next if ($.==1);

    my ($id, @allfield) = split(/\t/, $_);
    my %dom;

    for (my $if=$field[0]-2; $if <= $field[$#field]-2; $if++) {

        my $f = $allfield[$if];

        if ($f ne 'N') {
            $f=~s/\([^\)]+\)//g;
            $f=~s/\_[0-9A-Za-z]+//g;
            $f=~s/(?:\d+\.\d+)+(?:\.\d+)*//;
            $f=~s/^\+//;
            $f=~s/\+$//;
            if ($f) {
                my @eachf = split(/\+/, $f);
                foreach my $ef (@eachf) {
                    $dom{$ef}->{$if} = 0 unless (exists $dom{$ef});
                    $dom{$ef}->{$if}++;
                }                    
            }                
        }            
    }
    
    my @domains;
    foreach my $k (keys %dom) {
        push(@domains, $k) if (scalar(keys %{$dom{$k}})>=$threshold);
    }
    if (scalar(@domains)>=1) {
        print $id,"\t",join(";", @domains),"\t",join(";", map {$_.'('.scalar(keys %{$dom{$_}}).')'} @domains),"\n";
    }        
}    

close(IN);


# Subroutines

sub Usage {
    my ($msg) = @_;
    my $threshold=THRESHOLD_DEFAULT;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file
        -f      --field     Domain field numbers delimited by comma (at least one)
        -t      --threshold Minimum number of tools identifyiing a domain [Default: $threshold]        

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

