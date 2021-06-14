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

use Getopt::Long;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $alternative);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "a|alternative"=>\$alternative
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

my @oflag = (
	   	'superkingdom',
		'kingdom',
		'phylum',
		'class',
		'order',
		'family',
		'genus',
		'species'
	);

my %flag = (	
	   	'superkingdom'=>'d',
		'kingdom'=>'k',
		'phylum'=>'p',
		'class'=>'c',
		'order'=>'o',
		'family'=>'f',
		'genus'=>'g',
		'species'=>'s',
		'subspecies/strain'=>'t',
		'subspecies'=>'S',
		'strain'=>'T'
	);



print join("\t", '#OTU ID', 'taxonomy'),"\n";

if ($alternative) {

	open(IN, "-|" , 'cut -f 2,3 '.$infile.' | taxonkit --data-dir '.$ENV{'KRAKEN2_DB_PATH'}.'/'.$ENV{'KRAKEN2_DEFAULT_DB'}.'/taxonomy lineage --taxid-field 2 -R 2>/dev/null') or $LOGGER->logdie($!);

	while(<IN>) {
		chomp;
		my @F=split(/\t/, $_);
		$F[0]=~s/\s//g;
		
		my @r;
		my (@lineage, @rank);
		if ($F[2]) {
			@lineage=split(/;/, $F[2]);
			@rank=split(/;/, $F[3]);
		} else {
			for (my $i=0; $i<=$#oflag; $i++) {
				$lineage[$i]='';
				$rank[$i]='';
			}
		}
		my %lr;
		for (my $i=0; $i<=$#lineage; $i++) {
			$lr{$rank[$i]} = $lineage[$i]||'';
		}
		
		foreach my $o (@oflag) {
			push(@r, $flag{$o}.'__'.($lr{$o}||'') );
		}
		print join("\t",$F[0], join('; ', @r)),"\n";
	}
	close(IN);
	
} else {
	open(IN, "-|",'cut -f 2,3 '.$infile.' | taxonkit --data-dir '.$ENV{'KRAKEN2_DB_PATH'}.'/'.$ENV{'KRAKEN2_DEFAULT_DB'}.'/taxonomy lineage --taxid-field 2 | taxonkit --data-dir '.$ENV{'KRAKEN2_DB_PATH'}.'/'.$ENV{'KRAKEN2_DEFAULT_DB'}.'/taxonomy reformat  --lineage-field 3  -f "k__{k}; p__{p}; c__{c}; o__{o}; f__{f}; g__{g}; s__{s}" 2>/dev/null') or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my @F=split(/\t/, $_);
		
		$F[0]=~s/\s//g;
		$F[3]||='k__; p__; c__; o__; f__; g__; s__';
                #$F[3]=~s/\|s__$//;
                #$F[3]=~s/\|g__$//;
                #$F[3]=~s/\|f__$//;
                #$F[3]=~s/\|o__$//;
                #$F[3]=~s/\|c__$//;
                #$F[3]=~s/\|p__$//;
                #$F[3]=~s/\|d__$//;

                print join("\t",@F[0,3]),"\n";
	}
	close(IN);
}

# Subroutines

sub Usage {
    my ($msg) = @_;
    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help  	Help
        -l      --level		Log level [Default: FATAL]
	-i	--infile	kraken2 output file
	-a	--alternative	Set alternative method

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

