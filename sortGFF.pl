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
use FileHandle;
use POSIX 'isatty';

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $gffile);

#Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "g|gff=s"=>\$gffile
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

$infile=$gffile if ($gffile);

my $fhin;

if ($infile) {

    $LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

    $fhin = FileHandle->new;
    $fhin->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fhin = \*STDIN;
    } else {
        $LOGGER->logdie("Missing input file (-i/--infile) or STDIN data");
    }
}


my %data;

my %gff;
my %row;

my @root;


while(<$fhin>) {
    chomp;
    next if ($_=~/^#/);
    my @data=split(/\t/, $_);
    my ($id) = $data[8]=~/ID=([^;]+)/;
    my ($parent) = $data[8]=~/Parent=([^;]+)/;
    
    $gff{$id} = \@data;

    $row{$id} = undef;

    if ($parent) {
        push(@{ $data{$parent} }, $id);
    } else {
        push(@root, $id);
    }
}

$fhin->close();

print '##gff-version 3',"\n";
# gene
foreach my $id ( sort { $gff{$a}->[0] cmp $gff{$b}->[0] or $gff{$a}->[3] <=> $gff{$b}->[3] or $gff{$a}->[4] <=> $gff{$b}->[4] } @root ) {
    print join("\t", @{ $gff{ $id } }),"\n";
    delete($row{$id});
    if ($data{$id}) {
        # rna
        foreach my $subid ( sort { $gff{$a}->[0] cmp $gff{$b}->[0] or $gff{$a}->[3] <=> $gff{$b}->[3] or $gff{$a}->[4] <=> $gff{$b}->[4] } @{ $data{$id} }) {
            print join("\t", @{ $gff{ $subid } }),"\n";
            delete($row{$subid});
            # exon e cds
            if ($data{ $subid }) {
                foreach my $subsubid ( sort { $gff{$a}->[0] cmp $gff{$b}->[0] or $gff{$a}->[3] <=> $gff{$b}->[3] or $gff{$a}->[4] <=> $gff{$b}->[4] } @{ $data{$subid} }) {
                    print join("\t", @{ $gff{ $subsubid } }),"\n";
                    delete($row{$subsubid});
                }
            } else {
                $LOGGER->logdie("Not found children for $subid");
            }
        }
    } else {
        $LOGGER->logdie("Not found children for $id");
    }        
}

foreach my $id (keys %row) {
    $LOGGER->logwarn("Not found parent/child relation for $id");
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

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file (GFF) [Default: STDIN]
        -g      --gff       Input file (GFF) *** To be removed ***

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

