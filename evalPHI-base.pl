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

use Storable;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

use FileHandle;

my ($level,$infile,$outfile,$dumpdir,$phifile,$taxon_id);
# Command line named arguments (see function Usage() )
Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "o|outfile=s"=>\$outfile,
            "d|dumpdir=s"=>\$dumpdir,
            "p|phifile=s"=>\$phifile,
            "t|taxon_id=i"=>\$taxon_id
    ) or &Usage();


#Log4perl levels
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

$LOGGER->logdie("Missing input file with PHI-base search results") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

my %inh;
open( IN, "<", $infile ) or $LOGGER->logdie($!);
while (<IN>) {
    chomp;
    next if ( $_ =~ /^#/ );
    my ( $id, @F ) = split( /\t/, $_ );
    $inh{$id} = \@F;
}
close(IN);

$dumpdir||=".";

my $hr_phibase;
    
if ($phifile) {
    
    unless (-e $phifile) {
        $LOGGER->logdie("Missing PHI-base csv file");
    }

    open(PHI, "<", $phifile) or $LOGGER->logdie($!);
    { 
        local $/ = "";
        my $header_line=<PHI>;
        chomp($header_line);
        my @header=split(/,/, $header_line);
        foreach my $h (@header) {
            $h=~s/\s+//g;
        }
        while(my $l=<PHI>) {
            $l=~s/\n//g;
            chomp($l);
            my %data;
            my @m;
            $l=~s/""//g;
            while ($l=~ /("[^"]+")/g) {
                push(@m, $1);
            }
            foreach my $e (@m) {
                my $m=$e;
                $m=~s/,/;/g;
                my $qe=quotemeta($e);
                $l=~s/$qe/$m/;
            }
            if ($l) {            
                @data{@header} = split(/,/, $l);
                $LOGGER->logdie("Not found PHIMolConnID in this record ($.): $l") unless ($data{'PHIMolConnID'});
                $LOGGER->logdie("Not found the right PHIMolConnID (found: $data{'PHIMolConnID'}) in this record ($.): $l") unless ($data{'PHIMolConnID'}=~/^PHI:\d+$/);
                $LOGGER->logdie("Not found HostID in this record ($.): $l") unless ($data{'HostID'});
                if (($data{'HostID'} eq 'No host tests done')||($data{'HostID'} eq 'Lethal pathogen phenotype')) {
                    $data{'HostID'}=0;
                }
                $data{'HostID'}=~s/^[^0-9]+|[^0-9]+$//;
                $LOGGER->logdie("Not found the right HostID (found: $data{'HostID'}) in this record ($.): $l") unless ($data{'HostID'}=~/^\d+$/);
                $hr_phibase->{ $data{'PHIMolConnID'} }->{'HostID'} = $data{'HostID'};
            }                
        }
    }        
    close(PHI);
    store $hr_phibase, "$dumpdir/evalPHI-base.dump";
} else {

    unless (-e "$dumpdir/evalPHI-base.dump") {
        $LOGGER->logdie("Missing PHI-base dump or csv file");
    }
   
}

if (-e "$dumpdir/evalPHI-base.dump") {
    # Retrieve the hash reference from a file
    $hr_phibase = retrieve("$dumpdir/evalPHI-base.dump");
}

#foreach my $p (keys %{ $hr_phibase }) {        
#    print $p,"\t",$hr_phibase->{$p}->{'HostID'}||'',"\n";
#}

open(OUT, ">", $outfile) or $LOGGER->logdie($!);

foreach my $id (sort { $inh{$b}->[1]<=>$inh{$a}->[1] } keys %inh) {
    my @r_info=split(/#/, $inh{$id}->[0]);
    my $set=undef;
    foreach my $phi_id ( split(/__/, $r_info[1]) ) {
        $LOGGER->logdie("Not found PHI identifier ($phi_id) on database") unless (exists $hr_phibase->{$phi_id});
        my $r_taxon_id=$hr_phibase->{$phi_id}->{'HostID'};
        $LOGGER->logdie("Not found HostID for PHI identifier ($phi_id) on database") unless (defined $r_taxon_id);
        my $test_cmd = 'echo -e "'.$r_taxon_id.'\n'.$taxon_id.'" | GetLCA.pl 2>/dev/null';
        my $test_id = `$test_cmd`;
        chomp($test_id);
        if (($test_id==$taxon_id)||
            ($test_id==$r_taxon_id)) {
            $set=1;
            last;
        }            
    }
    print OUT $id,"\t",join("\t", @{$inh{$id}}),"\n" if ($set);
}

close(OUT);

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
        -i      --in        Input file of PHI-base search results 
        -e      --excludes  Exclude lineages (comma separated names)
        -o      --outfile   Output file with filtered PHI-base search results considering Host taxonomy clade
        -d      --dumpdir   Dump directory [Default: ./]
        -p      --phifile   PHI-base csv file
        -t      --taxon_id  NCBI taxonomy for host taxonomy clade filtering

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

