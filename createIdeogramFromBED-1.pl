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

use File::Basename;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $gfile, $infile, $outbase, $title);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "g|genome=s"=>\$gfile,
            "i|infile=s"=>\$infile,
            "o|outbase=s"=>\$outbase,
            "t|title=s"=>\$title
    ) or &Usage();

use File::Basename;

$title||=basename($infile);

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

$outbase||="./ideogram";

use Bio::SeqIO;

my $in = Bio::SeqIO->new(-file=>$gfile, -format=>'FASTA');

use R;
use Cwd;

&R::initR("--silent");
&R::library("RSPerl");

&R::eval('suppressMessages(library(package="GenomicRanges"))');
&R::eval('suppressMessages(library(package="ggbio"))');

unless (-e "$gfile.RData") {
    my %glength;
    while (my $seq = $in->next_seq()) {
        $glength{$seq->display_id()} = $seq->length();
    }

#    foreach my $s (keys %glength) {
#            print $s,"\t",$glength{$s},"\n";
#    }

    my @kglength = sort { &compn($a, $b) } keys %glength;
    #print STDERR join(";",@kglength),"\n";
    &R::eval(
    'gr <<- GRanges(seqnames = Rle(c('.join(', ', map { '"'.$_.'"' } @kglength ).')),
                            ranges = IRanges(1, end = c('.join(', ', @glength{@kglength}).'),names=c('.join(', ', map { '"'.$_.'"' } @kglength ).'))              
                            )'
    );
    &R::eval('seqlengths(gr)[c('.join(', ', map { '"'.$_.'"' } @kglength ).')] <<- c('.join(', ', @glength{@kglength}).')');
    &R::eval('save(list = ls(all=TRUE,envir=globalenv()), file = "'."$gfile.RData".'")');
}
else {
    &R::eval('load(file="'."$gfile.RData".'",envir = globalenv())');
}

#&R::eval('print(head(gr))');

my %input;
my %ref;

open(IN, "<", $infile) or $LOGGER->logdie($!);
while(<IN>) {
    chomp;
    my (@f) = split(/\t/, $_);
    my ($chr, $start, $end, $peakname, $score, $strand) = split(/\t/, $_);
    $input{ $chr }->{basename($peakname)} = {'start' => $start, 'end'=> $end, 'score'=>$score, 'strand'=>(($strand eq ".") ? "*" : $strand)};
    $ref{$chr} = undef;
}
close(IN);
    

#foreach my $k (sort { &compn($a, $b) } keys %ref) {
#    print STDERR $k,"\n";
#}

my @refs = sort { &compn($a, $b) } keys %ref;

&R::eval(
        'grref <<- GRanges(seqnames = Rle(as.character(seqnames(gr[c('.join(', ', map { '"'.$_.'"' }  @refs).')]))),
                               ranges = ranges( gr[c('.join(', ', map { '"'.$_.'"' }  @refs).')] ),
                                                      '.((scalar(@refs)) ? 'chrnames = as.factor(c('.join(', ', map { '"'.$_.'"' }  @refs).'))' : 'chrnames = as.factor(as.character(seqnames(gr[c('.join(', ', map { '"'.$_.'"' }  @refs).')])))' ).'
                                                                                  )'
                                                                                      );


&R::eval('seqlengths(grref)[names(seqlengths(grref))] <<- seqlengths(gr)[names(seqlengths(grref))]');

my %result;

foreach my $chr (@refs) {
    foreach my $peakname (keys %{ $input{$chr} }) {
        push(@{ $result{'seq_id'} }, $chr );
        push(@{ $result{'start'} }, $input{$chr}->{$peakname}->{'start'} );
        push(@{ $result{'end'} }, $input{$chr}->{$peakname}->{'end'} );
        push(@{ $result{'strand'} }, $input{$chr}->{$peakname}->{'strand'} );
        push(@{ $result{'score'} }, $input{$chr}->{$peakname}->{'score'} );
        push(@{ $result{'ID'} }, $peakname);
#        print ">>>>".$input{$chr}->{$peakname}->{'score'}."\n";
#        print ">>>>".$input{$chr}->{$peakname}->{'strand'}."\n";
    }
}

close(IN);


&R::eval(
    'grmap <<- GRanges(seqnames = Rle(c('.join(', ', map { '"'.$_.'"' } @{ $result{'seq_id'} } ).')),
                            ranges = IRanges(c('.join(', ', @{ $result{'start'} }).'), end = c('.join(', ', @{ $result{'end'} }).')),
                            strand =Rle(strand(c('.join(', ', map { '"'.$_.'"' } @{ $result{'strand'} } ).'))),
                            score  =Rle(c('.join(', ', @{ $result{'score'} } ).')),
                            ID=Rle(c('.join(', ', map { '"'.$_.'"' } @{ $result{'ID'} } ).'))
                            )'
        );

my $pwd = getcwd();
my $pngfile1 = $outbase.'_karyogram.png';
my $pngfile2 = $outbase.'_circle1.png';
my $pngfile3 = $outbase.'_circle2.png';

&R::eval('seqlengths(grmap)[names(seqlengths(grmap))] <<- seqlengths(gr)[names(seqlengths(grmap))]');
&R::eval('setwd("'.$pwd.'")');
&R::eval('seqlevels(grref, force=TRUE) <<- seqlevels(grref)[order(as.numeric(gsub("Group|MT","",seqlevels(grref))))]');
&R::eval('save(list = ls(all=TRUE,envir=globalenv()), file = "'."test.RData".'")');

&R::eval('png(filename="'.$pngfile1.'", bg="white", res=300, width=3000, height=3000)
p<<-ggplot() + layout_karyogram(grref, fill="gray56") + layout_karyogram(grmap,geom="rect", color="gold", fill="red")+ggtitle("'.$title.'")
plot(p)
dev.off()');

&R::eval("grmapplus <<- grmap[which(as.character(strand(grmap))=='+')]");
&R::eval("grmapminus <<- grmap[which(as.character(strand(grmap))=='-')]");
&R::eval("grmapnone <<- grmap[which(as.character(strand(grmap))=='*')]");
&R::eval('
png(filename="'.$pngfile2.'", bg="white", res=300, width=3000, height=3000)
p <<- ggplot() + layout_circle(grref, geom = "text", aes(label = seqnames), vjust = 0, radius = 40, trackWidth = 2) + layout_circle(grref, geom = "ideo", fill = "azure2", radius = 30,trackWidth = 4)
p <<- p + layout_circle(grref, geom = "scale", size = 2, radius = 35, trackWidth = 2)
if ( dim(elementMetadata(grmapplus))[1] > 0 ) {
    p <<- p + layout_circle(grmapplus, geom = "rect", color = "steelblue", radius = 26, trackWidth = 3)
}
if ( dim(elementMetadata(grmapminus))[1] > 0 ) {
    p <<- p + layout_circle(grmapminus, geom = "rect", color = "seagreen", radius = 23, trackWidth = 3)
}
if ( dim(elementMetadata(grmapnone))[1] > 0 ) {
    p <<- p + layout_circle(grmapnone, geom = "rect", color = "gold", radius = 23, trackWidth = 3)
}

p <<- p+ggtitle("'.$title.'")
plot(p)
dev.off()');

&R::eval('
png(filename="'.$pngfile3.'", bg="white", res=300, width=3000, height=3000)
p <<- ggplot() + layout_circle(grref, geom = "text", aes(label = seqnames), vjust = 0, radius = 40, trackWidth = 2) + layout_circle(grref, geom = "ideo", fill = "azure2", radius = 30,trackWidth = 4)

p <<- p + layout_circle(grref, geom = "scale", size = 2, radius = 35, trackWidth = 2)

if ( dim(elementMetadata(grmapplus))[1] > 0 ) {
    p <<- p + layout_circle(grmapplus, geom = "bar", aes(color=score, y=score), color = "steelblue", radius = -26, trackWidth = 3)
}
if ( dim(elementMetadata(grmapminus))[1] > 0 ) {
    p <<- p + layout_circle(grmapminus, geom = "bar", aes(color=score, y=score), color = "seagreen", radius = -23, trackWidth = 3)
}
if ( dim(elementMetadata(grmapnone))[1] > 0 ) {
    p <<- p + layout_circle(grmapnone, geom = "bar", aes(color=score, y=score), color = "gold", radius = -29, trackWidth = 6)
}

p <<- p+ggtitle("'.$title.'")
plot(p)
dev.off()');

#

# Subroutines

sub compn {
    my ($x, $y) = @_;
    my ($xn,$yn);
    if ($x=~/^Group/) {
       ($xn)=$x=~/(\d+)$/;
       if ($x=~/^GroupUn/) {
            $xn+=1000;
       }
    } elsif ($x=~/^MT/) {
        $xn=0;
    } else {
        $LOGGER->logdie("Error: $x");
    }

    if ($y=~/Group/) {
       ($yn)=$y=~/(\d+)$/;
       if ($y=~/^GroupUn/) {
            $yn+=1000;
       }
    } elsif ($y=~/^MT/) {
        $yn=0;
    } else {
        $LOGGER->logdie("Error: $y");
    }
    return ($xn <=> $yn);
}

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
        -g      --genome    Genome multi-fasta file
        -i      --infile    Input file (bed coordinates)
        -o      --outbase   Output base filenames
        -t      --title     Plot title [Default: <Input file basename>]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}


