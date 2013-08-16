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

use File::Temp qw/tempdir/;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $input, $reffileq, $reffiles, $threshold, $exactfile, $beddir, $evalue_threshold);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|input=s"=>\$input,
            "rq|reffileq=s"=>\$reffileq,
            "rs|reffiles=s"=>\$reffiles,
            "t|threshold=f"=>\$threshold,
            "e|exactfile=s"=>\$exactfile,
            "b|beddir=s"=>\$beddir,
            "v|evalue=f"=>\$evalue_threshold
    ) or &Usage();

$evalue_threshold||=0.001;    

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

use Bio::SearchIO; 

my %refscoreq;
my %refscores;

$LOGGER->logdie("Missing exact file") unless ($exactfile);
$LOGGER->logdie("Wrong exact file ($exactfile)") unless (-e $exactfile);

$LOGGER->logdie("Missing reference file (query)") unless ($reffileq);
$LOGGER->logdie("Wrong reference file (query) ($reffileq)") unless (-e $reffileq);

$LOGGER->logdie("Missing reference file (subject)") unless ($reffiles);
$LOGGER->logdie("Wrong reference file (subject) ($reffiles)") unless (-e $reffiles);


$LOGGER->logdie("Missing threshold") unless ($threshold);
$LOGGER->logdie("Wrong threshold ($threshold)") if ( ($threshold > 1) || ($threshold < 0) );

my %exact;
open(EXACT, "<", $exactfile) or $LOGGER->logdie($!);
while(<EXACT>) {
    chomp;
    my ($q, $h) = split(/\t/, $_);
    foreach my $eh (split(';', $h)) {
        $exact{$q}->{$eh} = undef;
    }
    print $q,"\t",$h,"\t",'==',"\n";
}
close(EXACT);

open(REF, "<", $reffileq) or $LOGGER->logdie($!);

while(<REF>) {
    chomp;
    if ($_=~/Query=(\S+) Hit=(\S+) Score=(\S+) /) {
        my ($q, $h, $s) = ($1, $2, $3);
        if ($q ne $h) {
            $LOGGER->logdie("Wrong association found $_");
        }
        $refscoreq{$q} = $s;
    }
}

close(REF);

open(REF, "<", $reffiles) or $LOGGER->logdie($!);

while(<REF>) {
    chomp;
    if ($_=~/Query=(\S+) Hit=(\S+) Score=(\S+) /) {
        my ($q, $h, $s) = ($1, $2, $3);
        if ($q ne $h) {
            $LOGGER->logdie("Wrong association found $_");
        }
        $refscores{$q} = $s;
    }
}

close(REF);

my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $input);

use File::Basename;

$LOGGER->logdie("Missing beddir") unless ($beddir);

my ($b1, $b2) = split('_x_', basename($input, '.bls'));

my $bed1f = $beddir.'/'.$b1.'.bed';
my $bed2f = $beddir.'/'.$b2.'.bed';

$LOGGER->logdie("Not found input1 bed file ($bed1f)") unless (-e $bed1f);
$LOGGER->logdie("Not found input2 bed file ($bed2f)") unless (-e $bed2f);

while( my $result = $in->next_result ) {
  ## $result is a Bio::Search::Result::ResultI compliant object
  my $best_hit_score;
  my @bh_s;
  my @bh_q;
  
  while( my $hit = $result->next_hit ) {
      last if ( $hit->significance() > $evalue_threshold );

      my $set = undef;

      $LOGGER->logdie("Not found refscores for ".$hit->name())  unless ($refscores{$hit->name()});
      
      if ( (scalar(@bh_q) == 0)  && (scalar(@bh_s) == 0) ) {
          push(@bh_s, $hit);
          next;
      }          

      if ( abs($hit->bits()-($best_hit_score||$hit->bits())) <= 2) {
              
            if (($hit->bits()/$refscores{$hit->name()}) >= $threshold) {
                push(@bh_s, $hit);
                $best_hit_score||= $hit->bits();
                $set = 1;
            } elsif (scalar(@bh_s)==0) {
                $LOGGER->logdie("Not found refscoreq for ".$result->query_name()) unless ($refscoreq{$result->query_name()});
                if ( ($hit->bits()/$refscoreq{$result->query_name()}) >= $threshold) {
                    push(@bh_q, $hit);
                    $best_hit_score||= $hit->bits();
                    $set = 1;
                }
            }

      }

      unless ( $set ) {
              print STDERR "INSUFFICIENT","\t",
                           "Query=",   $result->query_name,
                           " Hit=",    $hit->name,
                           " Score=",  $hit->bits,
                           " Evalue=", $hit->significance,
                           " Threshold=", $threshold,
                           " Reference(s)=", $refscores{$hit->name()},
                           " Reference(q)=", $refscoreq{$result->query_name()},
                           "\n";
      }                           
  }



    my $class;
    if (scalar(@bh_s)>0) {
        $class='=~';
    }
    elsif (scalar(@bh_q)>0) {
        $class='~=';
    }
    else {
        $class='';
    }
    my @best_hit = (@bh_s, @bh_q);

    #print join(';', map {$_->name()} @best_hit),"\n";

    if ( scalar(@best_hit)==1 ) {
        ## $hit is a Bio::Search::Hit::HitI compliant object
        if (exists $exact{ $result->query_name() }) {
            foreach my $bh (@best_hit) {
               unless ( exists $exact{ $result->query_name() }->{$bh->name()} ) {
                    $LOGGER->logdie("Inconsistent match ! Distinct between analysis: ".join(';', keys %{ $exact{ $result->query_name() } })." and ".$bh->name());
               }
            }               
        }
        else {        
            print $result->query_name,"\t",$best_hit[0]->name(),"\t",$class,"\n";
        }           
    }
    elsif ( scalar(@best_hit)>1 ) {
        
        if (exists $exact{ $result->query_name() }) {
            my @inconsistent;
            foreach my $bh (@best_hit) {
               unless ( exists $exact{ $result->query_name() }->{$bh->name()} ) {
                   push(@inconsistent, $bh->name());
               }
            }
            if (scalar(@inconsistent) == scalar(@best_hit)) {
                $LOGGER->logdie("Inconsistent match ! Distinct between analysis: ".join(';', keys %{ $exact{ $result->query_name() } })." and ".join(';', @inconsistent ));
            }                    
        }
        else {
            my @r = map { $_->name() } @best_hit;
            my $x = &disambiguation($result->query_name(), $bed1f, $bed2f, @r);
            print $result->query_name,"\t",$x,"\t",$class,"\n";
        }            
    }
    else {
        print STDERR "NO MATCH","\t",
                     "Query=",   $result->query_name,
                     "\n";
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
        -i      --input     Input file
        -rq     --reffileq  Reference file (query)
        -rs     --reffiles  Reference file (subject)
        -e      --exactfile Exact file
        -b      --beddir    BED dir
	    -t      --threshold Threshold
        -v      --evalue    Evalue threshold [Default: 0.001]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub disambiguation {
    my ($s1, $bed1f, $bed2f, @r) = @_;

    my $tmpdir = tempdir( CLEANUP => 1 );

    if (scalar(@r) > 1) {

        my ($b1) = $s1=~/\|?([^\|]+)$/;

        $LOGGER->logdie("Not found ACC for $s1") unless ($b1);

        $b1=~s/_Group.*//;
        $b1=~s/P/R/;

        `grep '$b1' $bed1f | cut -f 1,2,3,3,4,5,6 > $tmpdir/BED1.bed`;
        my @x;
        if (-z "$tmpdir/BED1.bed") {
            $LOGGER->logwarn("Difficult to map $s1 to any of ".join(';', @r).' because '.$b1.' not found in '.$bed1f);
        }
        else {

            foreach my $p (@r) {
                my ($b2) = $p=~/\|?([^\|]+)$/;
                $LOGGER->logdie("Not found ACC for $p") unless ($b2);

                $b2=~s/_Group.*//;
                $b2=~s/P/R/;
                `grep '$b2' $bed2f | cut -f 1,2,3,4,5,6 > $tmpdir/BED2.bed`;
                if (-z "$tmpdir/BED2.bed") {
                    $LOGGER->logwarn("Difficult to map $s1 to any of ".join(';', @r).' because '.$b2.' not found in '.$bed2f);
                }
                else {
                    `intersectBed -s -b $tmpdir/BED1.bed -a $tmpdir/BED2.bed -f 0.7 -r -wo > $tmpdir/INTER.txt`;
                    open(INTER, "<", "$tmpdir/INTER.txt") or $LOGGER->logdie($!);
                    while(<INTER>){
                        chomp;
                        my (@d) = split(/\t/, $_);
                        $LOGGER->logdie("Not found overlap length") unless (defined $d[12]);
                        push(@x, [ $p, $d[12] ]);
                    }
                    close(INTER);
                }
            }
        }
        if (scalar(@x) == 1) {
            return($x[0]->[0]);
        }
        elsif (scalar(@x) > 1) {

            my @sel;
            my $best=0;
            foreach my $y (sort {$b->[1] <=> $a->[1]} @x) {
                if ($best <= $y->[1]) {
                    push(@sel, $y);
                    $best = $y->[1];
                }
                else {
                    last;
                }
            }

            if (scalar(@sel)==1) {
                return($sel[0]->[0]);
            }
            else {
                my %selhash;
                @selhash{ (map { $_->[0] } @sel) } = ( map { $_->[1] } @sel );

                $LOGGER->logwarn("Difficult to map $s1 to any of ".join(';', @r).' because '.join(';', map { $_."[".$selhash{$_}."]"  } keys %selhash).' have the same intersection length so we are considering all of them.');

                return (join(';', keys %selhash));
            }
 
            return(undef);
        }
        else {
            $LOGGER->logwarn('Not found pieces or sufficient intersection @r=('.join(',', @r).') b1='.$b1.' bed1f='.$bed1f.' bed2f='.$bed2f);
            return(join(';', @r));
        }
    }
    else {
        return($r[0]);
    }
}

