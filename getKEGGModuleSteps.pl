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
use Storable;
use File::Basename;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level);
my ($modfile,$dumpdir,$keggmod,$keggolist);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "f|moddef=s"=>\$modfile,
            "d|dumpdir=s"=>\$dumpdir,
            "m|mod=s"=>\$keggmod,
            "k|KO=s"=>\$keggolist
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

$dumpdir ||= '.';
my $moddump = $dumpdir.'/'.basename($0,'.pl').'.dump';

if ($modfile) {
    $LOGGER->logdie("Wrong module definition file ($modfile)") unless (-e $modfile);
    my %module;
    open(IN, "<", $modfile) or $LOGGER->logdie($!);
    while(<IN>) {
        chomp;
        my ($mod, $def) = split(/\t/, $_);
        $module{$mod} = $def;
    }
    close(IN);
    store \%module, $moddump;
    $LOGGER->info("Stored ".scalar(keys %module)." modules definitions ...");

}

$LOGGER->logdie("Missing module definition dump file") unless (-e $moddump);
my $hr_module = retrieve($moddump);

$LOGGER->info("Loaded ".scalar(keys %{$hr_module})." modules definitions ...");

my %K;
$LOGGER->logdie("Missing KEGG orthologs") unless ($keggolist);
$keggolist=~s/\s+//g;
@K{ split(/,/, $keggolist) } = ();

$LOGGER->logdie("Not found $keggmod") unless (exists $hr_module->{$keggmod});

my $definition = $hr_module->{$keggmod};

my %valid;
{
    my $definition_for_valid=$definition;
    $definition_for_valid=~s/[\(\)]//g;
    $definition_for_valid=~s/\+/ /g;
    $definition_for_valid=~s/\s+/,/g;
    @valid{ split(/,/, $definition_for_valid) } = ();
}    
$LOGGER->info("Found ".scalar(keys %valid)." KOs in the module $keggmod");
my $matched=0;
my @unmatched=();
foreach my $keggo (keys %K) {
    next if $keggo eq 'NA';
    if (exists $valid{$keggo}) {
        $matched++;
    } else {
        push(@unmatched, $keggo);
    }        
}
$LOGGER->info("There are ".$matched." matched KOs in the module $keggmod");
$LOGGER->logwarn( "WARNING: There is(are) ".scalar(@unmatched)." unmatched KO(s) in the module $keggmod: ".join(",", @unmatched)) if (scalar(@unmatched));


my @def=split(//, $definition);

my %stack;
my $lv=0;
my $n=0;
foreach my $d (@def) {

    if ($d eq '(') {
        $n++;
    }
        
    #    print $lv,"\t",$n,"\t",$d,"\n";
    for (my $w=$n; $w>=0; $w--) {
            #    print "\t",$lv,"\t",$w,"\t",$d,"\n";
            $stack{$lv}->[$w].=$d;
    }

    if ($d eq ')') {
        $n--;
     
    }        
        
    if (($d eq ' ')&&($n==0)) {
            $lv++;
    }
}

#foreach my $l (sort {$a <=> $b} keys %stack) {
#    for (my $n=0; $n<=$#{$stack{$l}}; $n++) {
#        print $l,"\t",$n,"\t",$stack{$l}->[$n],"\n";
#    }
#}

my $step_count=0;
my $mod_step_count=0;
foreach my $l (sort {$a <=> $b} keys %stack) {
    my $step = $stack{$l}->[0];
    $step=~s/^\s+//;
    $step=~s/\s+$//;

    my %element;

    my $defstr=$step;
    $defstr=~s/[\(\)]//g;
    $defstr=~s/[+ ]/,/g;
    my $defeval=$step;

    my %uniqueelement;
    @uniqueelement{ split(/,/, $defstr) } = ();
    foreach my $e ( keys %uniqueelement ) {
        $element{$e} = 1 if (exists $K{$e});
        $defeval=~s/$e/\$element{$e}/g;
    }
    $defeval=~s/\s+/&&/g;
    $defeval=~s/,/||/g;
    $defeval=~s/-/||/g;
    $defeval=~s/\+/&&/g;

    $LOGGER->info("Definition evaluation: $defeval");

    my $ret=((eval("($defeval)")) ? 1 : 0);
    $LOGGER->info("$step\t",$ret);
    $step_count++ if ($ret);
    $mod_step_count++;
}
print "$keggmod\t$step_count\t$mod_step_count","\n";



# Subroutines

sub Usage {
    my ($msg) = @_;
    
    my $dumpdir=basename($0,'.pl').'.dump';

	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -f      --moddef    File with modules definition (KEGG MODULE<TAB>MODULE DEFINITION)
        -d      --dumpdir   Directory path containing $dumpdir file [Default: current directory]
        -m      --mod       KEGG module
        -k      --KO        KEGG Orthologs separated by comma

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

