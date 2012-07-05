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

    fasta2fastq.pl - Convert sequences in fasta/qual files to fastq file format.

=head1 SYNOPSIS

    # create a fastq file from fasta and qual files
    
    $ fasta2fastq.pl -f infile.fasta -q infile.fasta.qual -fq outfile.fastq
    
    # create a fastq file from fasta file and quality default 25
    
    $ fasta2fastq.pl -f infile.fasta -qd 25 -fq outfile.fastq

    # The default ASCII offset is 33, the example below changes the offset to 64
    
    $ fasta2fastq.pl -f infile.fasta -qd 25 -o 64 -fq outfile.fastq

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Convert sequences in fasta/qual files to fastq file format using an
    specified quality file or a quality default value.

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

my ($fasta, $qual, $offset, $fastq, $qualdefault);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "f|fasta=s"=>\$fasta,
            "q|qual=s"=>\$qual,
            "qd|qualdefault=i"=>\$qualdefault,
            "o|offset=i"=>\$offset,
            "fq|fastq=s"=>\$fastq
    ) or &Usage();

$offset||=33;

die "Missing input fasta file" unless ($fasta);
die "Missing output fastq file" unless ($fastq);

die "Wrong input fasta file ($fasta)" unless (-e $fasta);

open(FASTA, "<", $fasta) or die $!;
open(FASTQ, ">", $fastq) or die $!;

$/ = ">";
{
    my $debris = <FASTA>;
}
my $count = 0;

if ($qual) {
    die "Wrong input qual file ($qual)" unless (-e $qual);
    open(QUAL, "<", $qual) or die $!;
    {
        my $debris = <QUAL>;
    }

    while (my $fastarec = <FASTA>) {
        $fastarec=~s/>$//;
        $fastarec=~s/\n+$//;

        my $qualrec = <QUAL>;
        die "Missing qual record" unless ($qualrec);
        $qualrec=~s/>$//;
        $qualrec=~s/\n+$//;

        my ($fastahead, $fastabody) = $fastarec=~/([^\n]+)\n(.*)/s;
        $fastabody=~s/\s+|\n+//sg;
        my ($fastaid) = $fastahead=~/^(\S+)/;
    #    print ">$fastaid\n$fastabody\n";
        
        my ($qualhead, $qualbody) = $qualrec=~/([^\n]+)\n(.*)/s;
        $qualbody=~s/\n+/ /sg;
        my ($qualid) = $qualhead=~/^(\S+)/;
    #    print ">$qualid\n$qualbody\n";
        $qualbody=~s/^\s+|\s+$//;
        my @qual = split(/\s+/, $qualbody);
    #    print join('-', @qual),"\n";

        die "FASTA and QUAL identifiers don't match each other" unless ($fastaid eq $qualid);
     
        print FASTQ '@',$fastaid,"\n",
                    $fastabody,"\n",
                    '+',$qualid,"\n",
                    join('', map { chr($_+$offset) } @qual),"\n";

        $count++;
        if ($count % 10000 == 0) {
            print STDERR "Converted records: $count                                                         \r";
        }
    }
    print STDERR "Converted records: $count                                                         \n";

    close(QUAL);
}
else {    
    unless (defined $qualdefault) {
        die "This script needs input qual file or quality default";
    }
    else {
        while (my $fastarec = <FASTA>) {
            $fastarec=~s/>$//;
            $fastarec=~s/\n+$//;

            my ($fastahead, $fastabody) = $fastarec=~/([^\n]+)\n(.*)/s;
            $fastabody=~s/\s+|\n+//sg;
            my ($fastaid) = $fastahead=~/^(\S+)/;
        #    print ">$fastaid\n$fastabody\n";
            
            my $qualid = $fastaid;
            my @qual = ( $qualdefault ) x length($fastabody);

            print FASTQ '@',$fastaid,"\n",
                        $fastabody,"\n",
                        '+',$qualid,"\n",
                        join('', map { chr($_+$offset) } @qual),"\n";

            $count++;
            if ($count % 10000 == 0) {
                print STDERR "Converted records with quality default ($qualdefault): $count                                                         \r";
            }
        }
        print STDERR "Converted records with quality default ($qualdefault): $count                                                         \n";

    }
}

close(FASTA);
close(FASTQ);

# Subroutines

sub Usage {
    my ($msg) = @_;
    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2010 Regional Blood Center of Ribeirão Preto

Usage

        $0	[-h/--help] [-f <input-fasta-file> -fq <output-fastq-file>] \
            [-q <input-qual-file> | -qd <quality-default>] [-o <ascii-offset>]

Argument(s)

        -h      --help          Help
        -f      --fasta         Input fasta file
        -q      --qual          Input qual file
        -qd     --qualdefault   Quality default
        -o      --offset        ASCII offset (default:33)
        -fq     --fastq         Output fastq file
    
END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

