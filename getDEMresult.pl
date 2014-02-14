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

use Math::Combinatorics;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $samples);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=> \$infile,
            "s|samples=s"=>\$samples
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

my %data;

$LOGGER->logdie("Missing samples") unless ($samples);

my @sample = split(/\,/, $samples);
for (my $i = 0; $i <= $#sample; $i++) {
	$sample[$i]=~s/^\s+|\s+$//;
}

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

open(IN, "<", $infile) or $LOGGER->logdie($!);

my @header;
my %line;

my %total;

while(<IN>) {
	chomp;
        my @v = split(/\t/, $_);
	if ($.==1) {
		@header = @v;
		my %test;
		@test{@header} = ();
		foreach my $s (@sample) {
			$LOGGER->logdie("Wrong sample ($s)") unless (exists $test{$s});
		}
	}
	else {
		@line{@header} = @v;
		foreach my $s (@sample) {
            $line{ $s }||=0;
            $line{ $s }+=1;
			$data{ $line{'#miRNA'} }->{ $s }+= $line{ $s };
            $total{ $s }+=$line{ $s };
		}
	}
}

close(IN);

#foreach my $s (keys %total) {
#    print "TOTAL $s: $total{$s}\n";
#}

#print join("\n", map { join " ", @$_ } combine(2,@sample)),"\n";

my @DE = combine(2,@sample);
my @DE_fold_header = map { $_->[0].'x'.$_->[1].'.fold.change' } @DE;
my @DE_pvalue = map { $_->[0].'x'.$_->[1].'.pvalue' } @DE;
my @DE_qvalue = map { $_->[0].'x'.$_->[1].'.qvalue' } @DE;


my @tmp;
my %mirnap;
foreach my $mirna (keys %data) {
    my @foldchange;
    my @pvalue;
    foreach my $cmp_ar ( @DE ) {
#        print $cmp_ar->[0],"\t",$cmp_ar->[1],"\n";

        my $v1 = $data{$mirna}->{ $cmp_ar->[0] };
        my $v2 = $data{$mirna}->{ $cmp_ar->[1] };
        
        my $t1 = $total{ $cmp_ar->[0] };
        my $t2 = $total{ $cmp_ar->[1] };
        
        my $vn1 = sprintf("%.2f", ((1000000/$t1)*$v1));
        my $vn2 = sprintf("%.2f", ((1000000/$t2)*$v2));
        
        unless ($data{$mirna}->{ $cmp_ar->[0].'.norm' }) {
            $data{$mirna}->{ $cmp_ar->[0].'.norm' } = $vn1;
        } 
        unless ($data{$mirna}->{ $cmp_ar->[1].'.norm' }) {
            $data{$mirna}->{ $cmp_ar->[1].'.norm' } = $vn2;
        } 

        push(@foldchange, sprintf("%.2f", &log2( ($vn1/$vn2) )));
        
        my @winflat = `winflat -xvalue $v1 -yvalue $v2 -diff $t1 $t2`;
        my $pvalue = 0;
        if ( $v2 < $v1 ) {
            ($pvalue) = $winflat[0] =~ /\) = (\S+)/;
        }
        else {
            ($pvalue) = $winflat[1] =~ /\) = (\S+)/;
        }

        push(@pvalue, sprintf("%.3f", $pvalue) );
        
        $mirnap{ $mirna }->{$cmp_ar->[0].'x'.$cmp_ar->[1].'.pvalue'} = $pvalue;
    }
    push(@tmp, [$mirna, @{ $data{$mirna} }{@sample}, @{ $data{$mirna} }{ map { $_.'.norm' } @sample} , @foldchange, @pvalue]);
}

#print join("\t", '#miRNA', @sample, @DE_fold_header, @DE_pvalue),"\n";

#foreach my $art (@tmp) {
#    my $mirna = $art->[0];
#    print join("\t", @{$art}),"\n";
#}


use File::Temp qw/ tempfile tempdir /;


use R;

use File::Basename;

&R::initR("--silent","--no-save");
&R::library("RSPerl");

$|=1;

my $tmpdir = tempdir( CLEANUP => 1 );
my ($tmpfh, $tmpfile) = tempfile( DIR => $tmpdir );

autoflush { $tmpfh } 1;

print { $tmpfh } join("\t", 'mirna', @DE_pvalue),"\n";
foreach my $mirna (keys %mirnap) {
    print { $tmpfh } join("\t", $mirna, @{$mirnap{$mirna}}{@DE_pvalue}),"\n";
}

#print $tmpdir,"\t",$tmpfile,"\n";
&R::eval('pvalue.df <<- read.table(file="'.$tmpfile.'", header=TRUE, stringsAsFactors = FALSE)');
&R::eval('analysis <<- setdiff(colnames(pvalue.df), "mirna")');
&R::eval('analysis.final <<- gsub("pvalue", "qvalue", analysis)');
&R::eval('analysis.final.tmp <<- paste(analysis.final,"tmp", sep=".")');
&R::eval('
for (a in 1:length(analysis.final.tmp)) {
   pvalue.df[, analysis.final.tmp[a] ] <<- p.adjust(pvalue.df[, analysis[a] ], method="fdr")
}        
');
&R::eval('
    if ( length(analysis.final.tmp) > 1) {
        pvalue.df[, analysis.final] <<- t(apply(pvalue.df[, analysis.final.tmp], 1, function(x){ return(as.numeric(p.adjust(x, method="fdr"))) } ))
    } else {
        pvalue.df[, analysis.final] <<- pvalue.df[, analysis.final.tmp]
    }    
    ');
#&R::eval('pvalue.df[, setdiff(colnames(pvalue.df), c(analysis.final.tmp, "mirna"))] <- round(pvalue.df[, setdiff(colnames(pvalue.df), c(analysis.final.tmp, "mirna"))],3)');

&R::eval('write.table(file="'."$tmpdir/Result.txt".'", pvalue.df[, setdiff(colnames(pvalue.df), analysis.final.tmp)], sep="\t", quote=FALSE, row.names=FALSE)');

open(IN, "<", "$tmpdir/Result.txt") or $LOGGER->logdie($!);
my $h = <IN>;
chomp($h);
my @hres = split(/\t/, $h);
my %mirnaq;
while(<IN>) {
    chomp;
    my %r;
    @r{ @hres} = split(/\t/, $_);
    foreach my $e (@hres) {
        if ($e =~ /\.qvalue$/) {
            $mirnaq{ $r{'mirna'} }->{ $e } = $r{$e};
        }            
    } 
}
close(IN);

print join("\t", '#miRNA', @sample, (map { $_.'.norm' } @sample), @DE_fold_header, @DE_pvalue, @DE_qvalue),"\n";

foreach my $art (@tmp) {
    my $mirna = $art->[0];
    print join("\t", @{$art}, map { sprintf("%.3f", $_) } @{ $mirnaq{$mirna} }{@DE_qvalue}),"\n";
}



#&R::eval('save(list = ls(all=TRUE,envir=globalenv()), file ="/tmp/debug.RData")');

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
        -i      --infile    Input file (output of miRDeep2 quantifier)
        -s      --samples   Sample list

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

sub log2 {
   my $n = shift;
   return log($n)/log(2);
}

