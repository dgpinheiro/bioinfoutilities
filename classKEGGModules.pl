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

my ($level,$infile,$modcol,$kcol,$restfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "m|modcol=i"=>\$modcol,
            "k|kcol=i"=>\$kcol,
            "r|kofile=s"=>\$restfile
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

my %restrict;


$modcol||=1;
$kcol||=2;

$LOGGER->logdie("Wrong modcol ($modcol)") if ($modcol < 1);
$LOGGER->logdie("Wrong kcol ($kcol)") if ($kcol < 1);

open(IN, "<", $infile) or $LOGGER->logdie($!);
my $header=<IN>;
chomp($header);
my @h=split(/\t/, $header);


if ($restfile) {
    $LOGGER->logdie("Wrong kofile ($restfile)") unless (-e $restfile);
    open(REST,"<",$restfile) or $LOGGER->logdie($!);
    while(<REST>) {
        chomp;
        $restrict{$_}=undef;
    }
    close(REST);    
    
    print join("\t", @h,'Steps_found','Steps_needed','Rest_KO','Rest_Steps_found'),"\n";

} else {
    print join("\t", @h,'Steps_found','Steps_needed'),"\n";
}

while(<IN>) {
    chomp;
    my @F = split(/\t/, $_);
    my $script_cmd='getKEGGModuleSteps.pl -m '.$F[$modcol-1].' -K '.$F[$kcol-1];
    my $script_res=`$script_cmd`;
    chomp($script_res);
    my ($mod_id, $steps_found, $steps_needed) = split(/\t/, $script_res);
    
    if ($restfile) {
        $F[$kcol-1]=~s/\s//g;
        my @ko=();
        foreach my $k ( split(/,/, $F[$kcol-1]) ) {
            push(@ko, $k) if (exists $restrict{$k});
        }
        my $k_param=join(',', @ko);
        $k_param||='NA';
        my $script_rest_cmd='getKEGGModuleSteps.pl -m '.$F[$modcol-1].' -K '.$k_param;
        my $script_rest_res=`$script_rest_cmd`;
        chomp($script_rest_res);
        my ($rest_mod_id, $rest_steps_found, $rest_steps_needed) = split(/\t/, $script_rest_res);

        print join("\t", @F, $steps_found, $steps_needed, join(',', @ko), $rest_steps_found),"\n" if ($rest_steps_found != $steps_found);

    } else {        
        print join("\t", @F, $steps_found, $steps_needed),"\n";
    }
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
        -i      --infile    Input file
        -m      --modcol    Column number for KEGG module ID [Default: 1]
        -k      --kcol      Column number for KEGG Orthologies [Default: 2]
        -r      --kofile    KOs file to restrict analysis

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

