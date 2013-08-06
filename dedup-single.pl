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

    dedup-single.pl - Remove duplicated sequences from structured file.

=head1 SYNOPSIS
    
    # remove duplicated sequences from structured file SE.txt with this separator "|" (-s)
    # SE1_seqid|SE1_seq|SE1_qualid|SE1_qual
    # Obs.: The lines of file (SE.txt) MUST BE sorted by sequence.

    $ dedup-single.pl -i SE.txt \
    > -o SE_cons.fastq -s "|";

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Remove duplicated sequences from a structured file with an specific separator (<sep>)
    containing these data P1_seqid<sep>sequence<sep>P1_qualid<sep>qual
    The lines of input file MUST BE sorted by sequence.

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

my ($outfile, $infile, $sep, $run_prefix);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"  => sub { &Usage(); },
            "o|outputfile=s"=> \$outfile,
            "i|inputfile=s"=> \$infile,
            "s|separator=s"=> \$sep,
	    "p|prefix=s"=>\$run_prefix
) or &Usage();

$sep||="\t";

$run_prefix||="SEQ";

die "Prefix must have 3 letters" unless ( length($run_prefix) == 3 );

our $to = (GenerateBase(36))[0];

die "Missing output file" unless ($outfile);

die "Missing input file" unless ($infile);
die "Wrong input file" unless (-e $infile);

open(OUT, ">", $outfile) or die "Cannot open output file: $!";

open(IN, "<", $infile) or die "Cannot open input file: $!";

my $count = 0;

$| = 1;


my $lastseq;
my @lastqual = ();

my $line = 0;

my $sequence_count = 0;

{
	my $fl = <IN>;
	chomp($fl);
	my ($seqid, $seq, $qualid, $qual) = split($sep, $fl);
	@lastqual = split(//, $qual);
	$lastseq = $seq;
    if ($seqid=~/_x(\d+)/) {
        $count+=$1;
    }
    else {
    	$count = 0;
    }        
}

while(<IN>) {
	chomp;
	my ($seqid, $seq, $qualid, $qual) = split($sep, $_);
	if ($lastseq ne $seq) {
		
		$sequence_count++;
		&print_seq(\*OUT, $sequence_count, $lastseq, \@lastqual, $count);

		@lastqual = split(//, $qual);
		$lastseq = $seq;
		$count = 0;
	}
	else {
		&update_qual(\@lastqual, $qual);
	}
    if ($seqid=~/_x(\d+)/) {
    	$count+=$1;
    }
    else {
        $count++;
    }
	$line++;
    	if ($line % 10000 == 0) {
        	print STDERR "Records: $line                                                         \r";
	}
}

$sequence_count++;
&print_seq(\*OUT, $sequence_count, $lastseq, \@lastqual, $count);

$line++;
print STDERR "Records: $line                                                         \n";

close(OUT);

close(IN);

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2010 Regional Blood Center of Ribeirão Preto

Usage

        $0	[-h/--help] [-i SE.txt 
                         -o SE_cons.fastq ] [-s "|"]

Argument(s)

        -h      --help          Help
        
        -o      --outputfile    Output file 

        -i      --inputfile     Input file (from merge-ends.pl)

        -s      --separator     Separator (default: TAB)
	
	-p	--prefix	Prefix of sequence ids (three letters)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub update_qual {
	my ($lqr, $qv) = @_;
	my @q = split(//, $qv);	
	for (my $i = 0; $i<length($qv); $i++) {
		if (ord($lqr->[$i]) < ord($q[$i])) {
			$lqr->[$i] = $q[$i];
		}
	}
}

sub print_seq {
	my ($fhr, $sequence_code, $s, $ar_qual, $c) = @_;
	#my $base36 = $to->($sequence_code);
	my $base36 = $sequence_code;
	my $sequence_name = (( '?' x (3-length($run_prefix))).$run_prefix).'_'.(( 0 x (6-length($base36))).$base36).'_x'.$c;
	print $fhr '@',$sequence_name,"\n",$s,"\n",'+',$sequence_name,"\n",join('', @{$ar_qual}),"\n";
}

# Change the base representation of a non-negative integer
sub GenerateBase
{
    my $base = shift;
    $base = 62 if $base > 62;
    my @nums = (0..9,'A'..'Z','a'..'z')[0..$base-1];
    my $index = 0;
    my %nums = map {$_,$index++} @nums;

    my $To = sub
    {
        my $number = shift;
        return $nums[0] if $number == 0;
        my $rep = ""; # this will be the end value.
        while( $number > 0 )
        {
            $rep = $nums[$number % $base] . $rep;
            $number = int( $number / $base );
        }
        return $rep;
    };

    my $From = sub
    {
        my $rep = shift;
        my $number = 0;
        for( split //, $rep )
        {
            $number *= $base;
            $number += $nums{$_};
        }
        return $number;
    };

    return ( $To, $From );
}
