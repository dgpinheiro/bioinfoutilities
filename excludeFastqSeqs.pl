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
#  Copyright (C) 2010  Fundação Hemocentro de Ribeirão Preto
#
#  Laboratório de Genética Molecular e Bioinformática
#  Núcleo de Bioinformática
#  BiT -  Bioinformatics Team
#  Fundação Hemocentro de Ribeirão Preto
#  Rua Tenente Catão Roxo, 2501
#  Ribeirão Preto - São Paulo
#  Brasil
#  CEP 14051-140
#  Fone: 55 16 39639300 Ramal 9603
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://lgmb.fmrp.usp.br
#
# $Id$

=head1 NAME

    repair-ends.pl - Repair sequence "ends" from paired-end experiments after filtering steps.

=head1 SYNOPSIS
    
    # original files before filtering.: PEExp_1.fastq / PEExp_2.fastq
    # files after filtering...........: PEExp_1_cleaned.fastq / PEExp_2_cleaned.fastq
    # output paired-end files.........: PEExp_1_cleaned_paired.fastq / PEExp_2_cleaned_paired.fastq
    # output single-end files.........: PEExp_cleaned_single.fastq

    $ repair-ends.pl    -g1 PEExp_1.fastq -g2 PEExp_2.fastq \
    > -i1 PEExp_1_cleaned.fastq -i2 PEExp_2_cleaned.fastq \
    > -o1 PEExp_1_cleaned_paired.fastq -o2 PEExp_2_cleaned_paired.fastq -os PEExp_cleaned_single.fastq ;

=head1 ABSTRACT
    
=head1 DESCRIPTION
    
    Script to pair the sequences of two ordered fastq files from paired-end 
    experiments after filtering steps. The pairing is based on their original
    files before filtering procedures.

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2010 Regional Blood Center of Ribeirão Preto

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Getopt::Long;

my ($infile, $excfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"  => sub { &Usage(); },
            "i|infile=s"=> \$infile,
            "e|excfile=s"=> \$excfile
) or &Usage();

die "Missing input file" unless ($infile);
die "Wrong input file" unless (-e $infile);

die "Missing exclude file" unless ($excfile);
die "Wrong exclude file" unless (-e $excfile);

my %exclude;

open(EXC, "<", $excfile) or die $!;
while(<EXC>) {
    chomp;
    $exclude{$_}=undef;
}
close(EXC);

$| = 1;
open(IN, "<", $infile) or die $!;

my $set;
while(<IN>) {
    chomp;
    if ( ($. % 4) == 1) {
        my ($id) = $_=~/^@(.+)/;

        die "Wrong pattern ($_) not found an id" unless (defined $id);

        if (exists $exclude{$id}) {
            $set=undef;
        }
        else {
            $set=1;
        }
    }

    if ($set) {
        print $_,"\n";
    }
}

close(IN);

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2010 Regional Blood Center of Ribeirão Preto

Usage

        $0	[-h/--help] [-i InputFastq.fastq 
                         -e ExcludeIDs.fastq ]

Argument(s)

        -h      --help          Help
        
        -i      --infile        Input file (Original fastq file)
        -e      --excfile       File with list of IDs to exclude

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

