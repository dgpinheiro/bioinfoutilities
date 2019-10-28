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
#  Copyright (C) 2018  Universidade Estadual Paulista - UNESP
#
#  Universidade Estadual Paulista - UNESP
#  Laboratório de Bioinformática
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
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

Copyright (c) 2018 Universidade Estadual Paulista - UNESP

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

my ($level,$gifile,$gffile,$taxid);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "g|gff=s"=>\$gffile,
            "i|gene_info=s"=>\$gifile,
            "t|tax_id=i"=>\$taxid
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

$LOGGER->logdie("Missing GFF") unless ($gffile);
$LOGGER->logdie("Wrong GFF ($gffile)") unless (-e $gffile);

$LOGGER->logdie("Missing gene_info file") unless ($gifile);
$LOGGER->logdie("Wrong gene_info file ($gifile)") unless (-e $gifile);

my %data;
open(GIF, "<", $gifile) or $LOGGER->logdie($!);
while(<GIF>) {
    next if ($.==1);
    chomp;
    my @F=split(/\t/, $_);
    if (($taxid) && ($taxid != $F[0])) {
        next;
    }
    $data{ $F[1] }->{'symbol'} = $F[2];
    $data{ $F[1] }->{'desc'} = $F[8];
}
close(GIF);

open(GFF, "<", $gffile) or $LOGGER->logdie($!);
while(<GFF>) {
    chomp;
    my @F = split(/\t/, $_);
    my @info = split(";", $F[8]) if ($F[8]);
    if ($info[0]) {
        next if ($info[0]!~/gene/); 
        $info[0]=~s/ID=//;  
        my ($gene_id)=$F[8]=~/GeneID:(\d+)/; 
        print join("\t",$info[0],(($gene_id) ? ($gene_id,$data{$gene_id}->{'symbol'},$data{$gene_id}->{'desc'})  : ('' x 3) )),"\n";
    }
}
close(GFF);

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
        -g      --gff       GFF
        -i      --gene_info gene_info
        -t      --tax_id    Taxonomy ID

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

