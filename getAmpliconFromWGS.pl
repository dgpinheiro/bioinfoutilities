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

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

use Bio::SeqIO;

my ($level, $infile, $overwrite, $reffile, $outfile, $threads, $minlen);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$infile,
            "f|force"=>\$overwrite,
            "r|reffile=s"=>\$reffile,
            "o|outfile=s"=>\$outfile,
            "m|minlen=i"=>\$minlen,
            "t|threads=i"=>\$threads
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

$threads||=1;
$minlen||=10;

$LOGGER->logdie("Missing input SAM file") unless ($infile);
$LOGGER->logdie("Wrong input SAM file ($infile)") unless (-e $infile);

$LOGGER->logdie("Missing reference (fasta/fastq) file") unless ($reffile);
$LOGGER->logdie("Wrong reference (fasta/fastq) file ($reffile)") unless (-e $reffile);


use FileHandle;

my $fhout;

unless ($outfile) {
	$fhout = \*STDOUT;
} else {
	$fhout = FileHandle->new(">".$outfile);
        unless (defined $fhout) {
		$LOGGER->logdie("Not defined output file handle for $outfile: $!");
    	}
}

use Storable;

use File::Temp qw/tempdir tempfile/;

my $hr_data = {};

if ( (! defined $overwrite) && (-e "$infile.dump" )) {
	$hr_data = retrieve("$infile.dump");
} else {

	open(IN, "-|", "sort -u $infile") or $LOGGER->logdie($!);

	while(<IN>) {
		chomp;
		#print $.,"\t",$_,"\n";
		my (@field) = split(/\t/, $_);
		#E786R(V4)	256	ERR646488.41474243.1	40
		my $primer;
		if ($field[0]=~/\d+R/) {
			$primer = 'R';
		} elsif ($field[0]=~/\d+F/) {
			$primer = 'F';
		} else {
			$LOGGER->logdie("Couldn't find primer for $field[0]");
		}
		$hr_data->{ $field[2] }->{$primer}->{$field[0]} = [ $field[1], $field[3], length($field[9]) ];
	}

	close(IN);

	store $hr_data, "$infile.dump";
}


my  $tmpdir = tempdir( CLEANUP => 1 , DIR=>'./');
my  ($tmpfh, $tmpfilename) = tempfile( DIR => $tmpdir, CLEANUP=>1 );

foreach my $read (keys %{ $hr_data }) {
	print { $tmpfh } $read,"\n";
}

`pullseq -i $reffile -n $tmpfilename > $tmpdir/out.fastq`;

my $seqin = Bio::SeqIO->new(-file=>"$tmpdir/out.fastq",-format=>'FASTQ');
my $seqtmp = Bio::SeqIO->new(-file=>">"."$tmpdir/out_processed.fastq",-format=>'FASTQ');

my @selected;
while(my $seq = $seqin->next_seq() ) {
	my $seq_id = $seq->display_id();
	my $seq_str = $seq->seq();
	#print ">",$seq_id,"\n",$seq_str,"\n";
	my @p = keys %{ $hr_data->{$seq_id} };
	if (scalar(@p) > 1) {
		$LOGGER->logwarn("Found a read with more than one primer match:", $seq_id);
	} else {
		my @pname = keys %{ $hr_data->{$seq_id}->{$p[0]} };
		if (scalar(@pname) > 1) {
			$LOGGER->logwarn("Found a read with more than one primer $p[0] match:", $seq_id);
		} else {
			my ($bwf, $start, $len) = @{ $hr_data->{$seq_id}->{$p[0]}->{$pname[0]} };
			my $rev_bin = reverse(&dec2bin($bwf));
			my @bwf_info = split(//, $rev_bin);
			#print "\t", $pname[0], "\t", $bwf, ' [', $rev_bin,'] ', "\t", $start,"\n";
			my ($subseq);
			if ( $bwf_info[4] ) {
				# read reverse strand
				($subseq) = $seq->trunc(1, $start+$len-1);
			} else {
				# read forward strand
				($subseq) = $seq->trunc($start,$seq->length());
			}
			if ($subseq->length() > $minlen) {
				$subseq->display_id($seq_id.'.Amp');
				my $seq_desc = $subseq->description();
				$subseq->description($seq_desc.' Amplicon ('.$pname[0].')');
				
				$seqtmp->write_seq($subseq);

				push(@selected, $subseq);
			}
		}
	}
}

my $eachrefn=1000000;
open(PULLSEQ, "-|", "pullseq -i $reffile -n $tmpfilename -e -c") or $LOGGER->logdie($!);
my $seqinps = Bio::SeqIO->new(-fh=>\*PULLSEQ, -format=>'FASTA');
my $seqc=0;
my $seqoutpscount = 1;
my @outps = ("$tmpdir/out_exc_$seqoutpscount.fa");
my $seqoutps = Bio::SeqIO->new(-file=>">$tmpdir/out_exc_$seqoutpscount.fa", -format=>'FASTA');
while(my $seqps = $seqinps->next_seq()) {
	if ( (int($seqc/$eachrefn)+1) != $seqoutpscount) {
		$seqoutpscount=int($seqc/$eachrefn)+1;
		push(@outps, "$tmpdir/out_exc_$seqoutpscount.fa");
		$seqoutps = Bio::SeqIO->new(-file=>">$tmpdir/out_exc_$seqoutpscount.fa", -format=>'FASTA');
	}
	$seqoutps->write_seq($seqps);
	$seqc++;
}
close(PULLSEQ);

#`parasail_aligner -c 20 -d -l 30 -s 99 -q $tmpdir/out_processed.fasta -g $tmpdir/paraout.txt -f $tmpdir/out_exc.fasta -t $threads`;

#open(PA, "<", "$tmpdir/paraout.txt") or $LOGGER->logdie($!);
#while(<PA>) {
#	my ($index1, $index2, $length1, $length2, $score, $end_query, $end_ref, $matches, $similarities, $length) = split(/,/, $_);
#}
#close(PA);

$seqoutpscount = 1;
foreach my $seqoutfa (@outps) {
	`bbmap.sh -Xmx22g ref=$seqoutfa in=$tmpdir/out_processed.fastq outm=$tmpdir/BBmap_$seqoutpscount.sam nodisk cigar=false sam=1.3 noheader=true touppercase=true slow=true minid=0.10  local=true threads=$threads ambiguous=all overwrite=true secondary=true maxsites=1000 trimreaddescriptions=t idfilter=0.80 mappedonly=true`;
	$seqoutpscount++;
}

my $seqout = Bio::SeqIO->new(-fh=>$fhout,-format=>'FASTQ');

foreach my $subseq (@selected) {
	$seqout->write_seq($subseq);
}

`cat $tmpdir/BBmap_*.sam | cut -f 3 | nsort -u > $tmpdir/BBmap_list.txt`;

open(PULLSEQ, "-|", "pullseq -i $reffile -n $tmpdir/BBmap_list.txt") or $LOGGER->logdie($!);
my $seqinps2 = Bio::SeqIO->new(-fh=>\*PULLSEQ, -format=>'FASTQ');
while(my $seqps = $seqinps2->next_seq()) {
	$seqout->write_seq($seqps);
}
close(PULLSEQ);

$fhout->close();


# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help  	Help
        -l      --level 	Log level [Default: FATAL]
	-i	--infile	Input SAM file (without header)
	-o	--outfile	Output fastq file
	-m	--minlen	Minimum length [Default: 10]
	-t	--threads	Threads [Default: 1]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

