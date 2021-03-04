#!/usr/bin/env perl
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
#  Copyright (C) 2019  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
#  Faculdade de Ciências Agrárias e Veterinárias (FCAV)
#  Laboratório de Bioinformática (LB)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://www.fcav.unesp.br
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

Copyright (C) 2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

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

use POSIX 'isatty';
use FileHandle;

my ($level, $nthreads, $infile);

GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "t|threads=i"=>\$nthreads
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

$nthreads||=2;

my $fh;

if ($infile) {

    die "Wrong input file ($infile)" unless (-e $infile);

    $fh = FileHandle->new;
    $fh->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fh = \*STDIN;
    } else {
        die "Missing input file (-i/--infile) or STDIN data [taxonomy assignment matrix]";
    }
}

while (<$fh>) {
    chomp;
    next if ($_=~/^#/);

    my ($seqid, @taxid) = split(/\t/, $_);
    my @vtaxid;
    foreach my $tid (@taxid) {
        push(@vtaxid, $tid) if ($tid);
    }

    my $btid= &recLCA(\@vtaxid, undef);
    if (defined $btid) {
        my $taxinfo=`echo $btid | taxonkit lineage -j $nthreads -t`;
        chomp($taxinfo);
    
        print $seqid,"\t",$taxinfo,"\t",@taxid,"\n";
    } else {
        print STDERR "WARNING: Not found $seqid (".join(';', @taxid).")\n";
    }
}

$fh->close();

sub recLCA {
    my ($ar_tid, $last)=@_;
    
    if (scalar(@{$ar_tid}) > 1) {
        #        print "[".join(";", @{$ar_tid})."]\n";
        my $lca_cmd='echo "'.join(' ',@{$ar_tid}).'" | GetLCA.pl -s';
        my $lca_taxid=`$lca_cmd`;
        chomp($lca_taxid);
        #        print ">>>$lca_taxid\n";
        
        if (($last)&&($lca_taxid == $last)) {
            return $lca_taxid;
        } else {
            my @valid;
            foreach my $id (@{$ar_tid}) {
                next if ($lca_taxid==$id);
                push(@valid, $id);
            }
            if (scalar(@valid)) {
                #                print ">>>>".join(";", @valid)," [$lca_taxid]\n";
                return &recLCA(\@valid, $lca_taxid);
            } else {
                return $lca_taxid;
            }                
        } 
    } else {
        #        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>$ar_tid->[0]\n";
        return $ar_tid->[0];
    }
}


# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file
        -t      --threads   Number of threads [Default: 2]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

