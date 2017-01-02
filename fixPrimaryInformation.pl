#!/usr/bin/perl

use strict;
use warnings;

use POSIX 'isatty';
use FileHandle;

my $infile = $ARGV[0];

my $fh;

if ($infile) {

    die "Wrong input file ($infile)" unless (-e $infile);

    $fh = FileHandle->new;
    $fh->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fh = \*STDIN;
    } else {
        die "Missing input file (-i/--infile) or STDIN data";
    }
}



my $id='';
my $pacount_first = 0;
my $pacount_second = 0;
my @aln_first;
my @aln_second;

while(<$fh>) {
    chomp;
    if ($_=~/^@/) {
        print $_,"\n";
        next;
    }

    my @F=split(/\t/, $_);

    if ($id ne $F[0]) {
        if ($id ne '') {
            print STDERR $id,"\t",$pacount_first,"\t",$pacount_second,"\n";
            &consume(\@aln_first, \@aln_second);
        }

        $id=$F[0];

        $pacount_first=0;
        $pacount_second=0;

        @aln_first=();
        @aln_second=();
    }
    
    my $revflag=reverse(sprintf("%012b",$F[1]));
    
    push(@aln_first, [substr($revflag,8,1), \@F]) if ( substr($revflag,6,1)==1 );
    push(@aln_second,[substr($revflag,8,1), \@F]) if ( substr($revflag,7,1)==1 );

    if ( substr($revflag,8,1)==0 ) {
        #print "\t",$F[0],"\t",$F[1],"\n";
        $pacount_first++ if ( substr($revflag,6,1)==1 );
        $pacount_second++ if ( substr($revflag,7,1)==1 );
    }
}

$fh->close();

print STDERR $id,"\t",$pacount_first,"\t",$pacount_second,"\n";
&consume(\@aln_first, \@aln_second);

sub consume {
    my ($ar_first, $ar_second) = @_;
    
    my @first = sort { $a->[0] <=> $b->[0] } @{ $ar_first };
    
    &change_print($first[0]->[1], 8, 0);

    foreach (my $i=1; $i<=$#first; $i++) {
        &change_print($first[$i]->[1], 8, 1);
    }

    my @second = sort { $a->[0] <=> $b->[0] } @{ $ar_second };
    
    &change_print($second[0]->[1], 8, 0);

    foreach (my $i=1; $i<=$#aln_second; $i++) {
        &change_print($second[0]->[1], 8, 1);
    }
} 
    
sub change_print {
    my ($ar_line,$bit,$binv) = @_;
    
    my $revflag=reverse(sprintf("%012b", $ar_line->[1]));
    substr($revflag, $bit, 1, $binv);
    my $adjflag=reverse($revflag);
    $ar_line->[1] = oct('0b'.$adjflag);
    
    print join("\t", @{ $ar_line }),"\n";
}    
