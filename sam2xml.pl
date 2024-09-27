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
use File::Basename;



use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile);

#Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=> \$infile
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

my $fh;

if ($infile) {

    $LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

    $fh = FileHandle->new;
    $fh->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fh = \*STDIN;
    } else {
        $LOGGER->logdie("Missing input file or STDIN data");
    }
}

my $diamond_version = 'v0.7.10.59';


my %query_length;
my %hit_length;

my %align;
my @qid;
my $hsp_num = 0;
my %hit_num;
my %subjct;
my %query;
while(<$fh>) {
    chomp;
 
    next if ($_ =~ /^@/);

    my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @extra) = split(/\t/, $_);   
    
    my ($bitscore, $rawscore, $reflength, $evalue, $identity, $qstart, $md, $frame);
    
    if ($extra[0] =~ /AS:i:(\d+)/) {
        $bitscore=$1;
    } else {
        $LOGGER->logdie('Not found AS flag');
    }

    if ($extra[2] =~ /ZL:i:(\d+)/) {
        $reflength=$1;
    } else {
        $LOGGER->logdie('Not found ZL flag');
    }

    if ($extra[3] =~ /ZR:i:(\d+)/) {
        $rawscore=$1;
    } else {
        $LOGGER->logdie('Not found ZR flag');
    }

    if ($extra[4]=~/ZE:f:(\S+)/) {
        $evalue=$1;
    } else {
        $LOGGER->logdie('Not found ZE flag');
    }

    if ($extra[5]=~/ZI:i:(\d+)/) {
        $identity=$1;
    } else {
        $LOGGER->logdie('Not found ZI flag');
    }

    if ($extra[6]=~/ZF:i:(-?\d+)/) {
        $frame=$1;
    } else {
        $LOGGER->logdie('Not found ZF flag');
    }

    if ($extra[7]=~/ZS:i:(\d+)/) {
        $qstart=$1;
    } else {
        $LOGGER->logdie('Not found ZS flag');
    }

    if ($extra[8]=~/MD:Z:(\S+)/) {
        $md=$1;
    } else {
        $LOGGER->logdie('Not found MD flag');
    }

    push(@qid, $qname) unless (exists $align{$qname});    

    if (exists $hit_num{$qname}) {
        unless (exists $align{$qname}->{$rname}) {
            $hit_num{$qname}++;
        }            
    } else {
        $hit_num{$qname} = 1;
    }

    unless (exists $align{$qname}->{$rname}->{'hsps'}) {
        $hsp_num = 0;
    }
    $hsp_num++;
    
    $align{$qname}->{$rname}->{'hit_num'} = $hit_num{$qname};
    
#    print $qname,"\t",$rname,"\n";   

    #print $seq,"\n";
    my @qseq = split(//, $seq);
    my @rseqtmp = &rmins(\@qseq, $cigar);

    my $qseqstr = '';
    my $rseqstr = '';
    
    my @ed;
    my $edc = 0;
    foreach my $e (split(//, $md)) {

        if ($e=~/\d/) {
            if (defined $ed[$edc]) {
                if (($ed[$edc]=~/^[A-Z]+$/)||($ed[$edc]=~/^\^/)) {
                    $edc++;
                }
            }                
            $ed[$edc].=$e;
        } elsif ($e=~/[A-Z]/i) {
            if (defined $ed[$edc]) {
                if ($ed[$edc]=~/^\d+$/) {
                    $edc++;
                }
            }                
            $ed[$edc].=$e;
        } elsif ($e=~/\^/) {
            if (defined $ed[$edc]) {
                if (($ed[$edc]=~/^\d+$/)||($ed[$edc]=~/^[A-Z]+$/)) {
                    $edc++;
                }
            }                
            $ed[$edc]=$e;
        }
    }
    
    my $rseqt = '';
    my $w = 0;
    foreach my $e (@ed) {
        #print STDERR "$qname / $rname >>>",$e,"\n";
        if ($e =~ /^(\d+)$/) {
            my $n = $1;
            for (my $i = 0; $i<$n; $i++) {
                $rseqt.=$rseqtmp[$w];
                $w++;
            }
        } elsif ($e =~/^([A-Z]+)$/ ) {
            $rseqt.=$1;
            my $el = length($1);
            $w+=$el;
        } elsif ($e=~/^\^([A-Z]+)/) {
            $rseqt.=$1;
            my $el = length($1);
            #$w+=$el;
        }
    }
    
    #$LOGGER->logdie("Error, not produced original rseq") unless ($rseqt);

    my @rseq = split(//, $rseqt);

    #$LOGGER->logdie("Error, not produced splitted original rseq") unless (scalar(@rseq));
    
    my ($pq, $pr) = (0, 0);
    my $alnlen = 0;
    my $ngaps = 0;
    while ($cigar=~/(\d+)([MDI])/g) {
        my ($n, $o) = ($1, $2);
        if ($o =~/M/) {
            for (my $i = 0; $i<$n; $i++) {
                $qseqstr.=$qseq[$pq];
                #$LOGGER->logdie("$qname $rname PR: $pr ($rseqt)") unless ($rseq[$pr]);
                $rseqstr.=$rseq[$pr];

                $pq++;
                $pr++;
            }
        } elsif ($o =~/D/) {
            for (my $i =0; $i<$n; $i++) {
                $qseqstr.='-';
                $rseqstr.=$rseq[$pr];
                $pr++;
                $ngaps++;
            }
        } elsif ($o =~/I/) {
            for (my $i =0; $i<$n; $i++) {
                $qseqstr.=$qseq[$pq];
                $rseqstr.='-';
                $pq++;
                $ngaps++;
            }
            
        }
        $alnlen++;
    }

#    print $qseqstr,"\n\n",$rseqstr,"\n\n\n";
    
    $hit_length{ $rname } = $reflength;
    $query_length{ $qname } = length($seq);

    push(@{ $align{$qname}->{$rname}->{'hsps'} }, {'hsp_num'=>$hsp_num, 'bit_score'=>$bitscore, 'raw_score'=>$rawscore, 'evalue'=>$evalue, 'query-from'=>$qstart, 'query-to'=>length($seq)+$qstart-1, 'hit-from'=>$pos, 'hit-to'=>$pos+length($rseqt)-1, 'hsp_id'=>$identity, 'hsp_gaps'=>$ngaps, 'align-len'=>$alnlen, 'frame'=>$frame, 'query-seq'=>$qseqstr, 'hit-seq'=>$rseqstr } );

    $subjct{ $rname } = undef;        
    $query{ $qname } = undef;

}
$fh->close();
$fh = undef;

printf("<?xml version=\"1.0\"?>\n");
printf("<BlastOutput>\n");
printf("  <BlastOutput_program>blastp</BlastOutput_program>\n");
printf("  <BlastOutput_version>diamond $diamond_version</BlastOutput_version>\n");
printf("  <BlastOutput_reference>Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nat Methods. 2015 Jan;12(1):59-60. doi: 10.1038/nmeth.3176. Epub 2014 Nov 17. PubMed PMID: 25402007.</BlastOutput_reference>\n");
printf("  <BlastOutput_db>n/a</BlastOutput_db>\n");
printf("  <BlastOutput_query-ID>n/a</BlastOutput_query-ID>\n");
printf("<BlastOutput_iterations>\n");

my $iter_num = 1;
foreach my $query_id (@qid) {

    printf("<Iteration>\n");
    printf("  <Iteration_iter-num>%d</Iteration_iter-num>\n", $iter_num++);
    printf("  <Iteration_query-ID>%s</Iteration_query-ID>\n", $query_id);
    printf("  <Iteration_query-def>%s</Iteration_query-def>\n", $query_id);
    $LOGGER->logdie("Not found length for query sequence $query_id") unless (exists $query{$query_id});
    printf("  <Iteration_query-len>%s</Iteration_query-len>\n", $query_length{ $query_id });
    printf("<Iteration_hits>\n");
    
    foreach my $hit_id (sort { $align{ $query_id }->{$a}->{'hit_num'} <=> $align{ $query_id }->{$b}->{'hit_num'} } keys %{ $align{ $query_id } }) {
            
        printf("<Hit>\n");
        printf("  <Hit_num>%d</Hit_num>\n", $align{ $query_id }->{ $hit_id }->{'hit_num'});
        printf("  <Hit_id>%s</Hit_id>\n", $hit_id );
        printf("  <Hit_def>%s</Hit_def>\n", $hit_id );
        $LOGGER->logdie("Not found length for database sequence $hit_id") unless (exists $subjct{$hit_id});
        printf("  <Hit_len>%d</Hit_len>\n", $hit_length{ $hit_id });
        printf("  <Hit_hsps>\n");
        
        foreach my $hr_hsp ( @{ $align{ $query_id }->{ $hit_id }->{'hsps'} }) {

            printf("    <Hsp>\n");
            printf("      <Hsp_num>%d</Hsp_num>\n", $hr_hsp->{'hsp_num'});
            printf("      <Hsp_bit-score>%f</Hsp_bit-score>\n", $hr_hsp->{'bit_score'});
            printf("      <Hsp_score>%f</Hsp_score>\n", $hr_hsp->{'raw_score'});
            printf("      <Hsp_evalue>%f</Hsp_evalue>\n", $hr_hsp->{'evalue'});
            printf("      <Hsp_query-from>%d</Hsp_query-from>\n", $hr_hsp->{'query-from'});
            printf("      <Hsp_query-to>%d</Hsp_query-to>\n", $hr_hsp->{'query-to'});
            printf("      <Hsp_hit-from>%d</Hsp_hit-from>\n", $hr_hsp->{'hit-from'});
            printf("      <Hsp_hit-to>%d</Hsp_hit-to>\n", $hr_hsp->{'hit-to'});
            printf("      <Hsp_query-frame>%d</Hsp_query-frame>\n",$hr_hsp->{'frame'});
            printf("      <Hsp_hit-frame>1</Hsp_hit-frame>\n");
            printf("      <Hsp_identity>%f</Hsp_identity>\n", $hr_hsp->{'hsp_id'});
            printf("      <Hsp_positive>?</Hsp_positive>\n");
            printf("      <Hsp_gaps>%d</Hsp_gaps>\n", $hr_hsp->{'hsp_gaps'});
            printf("      <Hsp_align-len>%d</Hsp_align-len>\n", $hr_hsp->{'align-len'});
            printf("      <Hsp_qseq>%s</Hsp_qseq>\n", $hr_hsp->{'query-seq'});
            printf("      <Hsp_hseq>%s</Hsp_hseq>\n", $hr_hsp->{'hit-seq'});
            printf("      <Hsp_midline>?</Hsp_midline>\n");
            printf("    </Hsp>\n");
        }

        printf("  </Hit_hsps>\n");
        printf("</Hit>\n");
            
    } 
    printf("</Iteration_hits>\n");
    printf("</Iteration>\n");
}

printf("</BlastOutput_iterations>\n");
printf("</BlastOutput>\n");

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
        -i      --infile    SAM file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub rmins {
    my @r;
    my ($ar_qseq, $cigarstr) = @_;
    my $p = 0;
    while ($cigarstr=~/(\d+)([MDI])/g) {
        my ($n, $o) = ($1, $2);
        if ($o eq 'M') {
            for (my $i=0; $i<$n; $i++) {
                push(@r, $ar_qseq->[$p++]);
            }
        } elsif ($o eq 'I') {
            for (my $i=0; $i<$n; $i++) {
                $p++;
            }
        }
    }
    return(@r);
}

