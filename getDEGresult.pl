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

my ($level, $infile, $samples, $excnotest, $lfc, $qvalue, $minfpkm);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=> \$infile,
            "s|samples=s"=>\$samples,
            "n|exnotest"=>\$excnotest,
            "g|lfc=f"=>\$lfc,
            "q|qvalue=f"=>\$qvalue,
            "m|minfpkm=f"=>\$minfpkm
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

my @data;

$LOGGER->logdie("Missing samples") unless ($samples);

my %valid;

my @sample = split(/\,/, $samples);
for (my $i = 0; $i <= $#sample; $i++) {
	$sample[$i]=~s/^\s+|\s+$//;
    $valid{$sample[$i]} = 1;
}

$LOGGER->logdie("Insufficient number of samples. Choose at least 2 samples.") if ( scalar( keys %valid ) < 2 );

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

open(IN, "<", $infile) or $LOGGER->logdie($!);

my @pvalue;

my @header;

my $first_line = <IN>;
chomp($first_line);

@header=split(/\t/, $first_line);

while(<IN>) {
    chomp;
    my %line;
    @line{@header} = split(/\t/, $_); 
    
    if ((exists $valid{ $line{'sample_1'} })&&(exists $valid{ $line{'sample_2'} } )) {

        $line{'value_1'}+=1;
        $line{'value_2'}+=1;
        
        my $log2fc = sprintf("%.6f", (log( $line{'value_2'}/$line{'value_1'} )/log(2)));
        $line{'log2(fold_change)'} = $log2fc;
        
        push(@data, \%line);
        push(@pvalue, [$line{'status'}, $line{'p_value'}]);
    }        
}

close(IN);

use File::Temp qw/ tempfile tempdir /;


use R;

use File::Basename;

&R::initR("--silent","--no-save");
&R::library("RSPerl");

$|=1;

my $tmpdir = tempdir( CLEANUP => 1 );
my ($tmpfh, $tmpfile) = tempfile( DIR => $tmpdir );

autoflush { $tmpfh } 1;

print { $tmpfh } join("\t", 'status', 'pvalue'),"\n";
foreach my $ar_p (@pvalue) {
    print { $tmpfh } join("\t", @{$ar_p}),"\n";
}

#print $tmpdir,"\t",$tmpfile,"\n";
&R::eval('pvalue.df <<- read.table(file="'.$tmpfile.'", header=TRUE, stringsAsFactors = FALSE)');
&R::eval('
   pvalue.df$qvalue <- NA
   ok <- which( pvalue.df$status == "OK" )
   pvalue.df[ok, "qvalue"] <<- round(p.adjust(pvalue.df[ok, "pvalue"], method="fdr"), 6)
');
&R::eval('write.table(file="'."$tmpdir/Result.txt".'", pvalue.df[,c("status", "pvalue", "qvalue")], sep="\t", quote=FALSE, row.names=FALSE)');


open(IN, "<", "$tmpdir/Result.txt") or $LOGGER->logdie($!);
my $h = <IN>;
chomp($h);
my @hres = split(/\t/, $h);
my $c = 0;
while(<IN>) {
    chomp;
    my %r;
    @r{ @hres} = split(/\t/, $_);
    if (($data[$c]->{'status'} eq 'OK')&&($r{'status'} eq 'OK' )) {
       $data[$c]->{'q_value'} = $r{'qvalue'};

# The code below was used to evaluate the qvalue calculation using R in comparison with cuffdiff results
#        if ($data[$c]->{'p_value'} == $r{'pvalue'}) {
#            if (sprintf("%.3f", $data[$c]->{'q_value'}) != sprintf("%.3f", $r{'qvalue'})) {
#                $LOGGER->logwarn("Check the qvalue for $data[$c]->{'test_id'}");
#            } else {
#                $data[$c]->{'q_value'} = $r{'qvalue'};
#                print "OK\n";
#            }
#        } else {
#            $LOGGER->logdie("Error found pvalue scrambled");
#        }
    }        
    $c++;
}
close(IN);

$lfc||=0;
$qvalue||=1;
$minfpkm||=0;

print join("\t", @header[0..($#header-1)], 'significant (cuffdiff)'),"\n";

foreach my $hr_data ( @data ) {
    if (($excnotest)&&($hr_data->{'status'} ne 'OK')) {
        next;
    }
    if (($lfc)&&( abs($hr_data->{'log2(fold_change)'}) < $lfc)) {
        next;
    }
    if (($qvalue)&&( $hr_data->{'q_value'} > $qvalue)) {
        next;
    }
    if ( ($hr_data->{'value_1'} >= $minfpkm) || ($hr_data->{'value_2'} >= $minfpkm) ) {
        print join("\t", @{$hr_data}{@header}),"\n";
    }
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
        -i      --infile    Input file (output of cuffdiff)
        -s      --samples   Sample list
        -n      --excnotest Exclude NOTEST/FAIL results
        -g      --lfc       Absolute log2-fold-change threshold [Default: 0]
        -q      --qvalue    Adjusted p-value threshold [Default: 1]
        -m      --minfpkm   Minimum FPKM in one or other sample [Default: 0]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

