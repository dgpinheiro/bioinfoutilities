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

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $indir, $formula, $useddh);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|indir=s"=>\$indir,
            "f|formula=i"=>\$formula,
            "d|ddh"=>\$useddh
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

$formula||=2;

my $method = (($useddh) ? 'ddh' : 'distance');

my %availform = (1=>{'distance'=>4,  'ddh'=>2 },
                 2=>{'distance'=>8,  'ddh'=>6 },
                 3=>{'distance'=>12, 'ddh'=>10}
                );

$LOGGER->logdie("Wrong formula ($formula)") unless (exists $availform{$formula});
$LOGGER->logdie("Wrong formula/method ($formula/$method)") unless (exists $availform{$formula}->{$method});

my %distance;
foreach my $f (glob("$indir/ggdc_results_?*.csv")) {
    $LOGGER->info("Processing $f");
    open(IN, "<", $f) or $LOGGER->logdie($!);
    while(<IN>) {
        chomp;
        next if ($.<=2);
        my @F = split(/,/, $_);
        $LOGGER->info("Data: ".join("\t", $F[0],$F[1],$F[$availform{$formula}->{$method}]));
        $distance{$F[0]}->{$F[1]} = $F[$availform{$formula}->{$method}];
        $distance{$F[1]}->{$F[0]} = $F[$availform{$formula}->{$method}];
    }
    close(IN);
}

my @ord = sort {$a cmp $b} keys %distance;

print join("\t", @ord),"\n";
foreach my $g (@ord) {
    $distance{$g}->{$g} = 0;
    print join("\t", $g, @{$distance{$g}}{@ord}),"\n";
}

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
        -i      --indir     Input directory (GGDC results directory)
        -f      --formula   Formula (1, 2 or 3) [Default: 2]
        -d      --ddh       DDH instead of distance [Default: Off]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

