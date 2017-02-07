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

my ($level, $iprfile, $tabfile, $colname, $newname, $choosen_analysis, $kegg_db, $cross_clstr);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|ipr=s"=>\$iprfile,
            "t|tab=s"=>\$tabfile,
            "c|column=s"=>\$colname,
            "a|analysis=s"=>\$choosen_analysis,
            "n|nanme=s"=>\$newname,
            "k|kegg_db=s"=>\$kegg_db,
            "x|cross_clstr=s"=>\$cross_clstr
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


$LOGGER->logdie("Missing interproscan tsv file") unless ($iprfile);
$LOGGER->logdie("Wrong interproscan tsv file ($iprfile)") unless (-e $iprfile);

$LOGGER->logdie("Missing tabular file") unless ($tabfile);
$LOGGER->logdie("Wrong tabular file ($tabfile)") unless (-e $tabfile);

$LOGGER->logdie("Missing column name") unless ($colname);


my %clstrmember;

if ($cross_clstr)  {
    $LOGGER->logdie("Wrong cluster file ($cross_clstr)") unless (-e $cross_clstr);
    my $cltid;
    open(CLT, "<", $cross_clstr) or $LOGGER->logdie($!);
    while(<CLT>) {
        chomp;
        if ( $_ =~ /^>(.+)/ ) {
            $cltid = $1;
            $cltid =~ s/\s/_/g;
            $cltid =~ s/Cluster/cluster/;
        }
        else {
            if ($cltid) {
                if ( $_ =~ />([^\.\s]+)/ ) {
                    my $member=$1;
                    push(@{$clstrmember{$cltid}},$member);
                }
            }
        }
    }
    close(CLT);

}

my $dbh;

if ($kegg_db) {
    require DBI;
    import DBI;
    $LOGGER->logdie("Wrong KEGG db file ($kegg_db)") unless (-e $kegg_db);
    $dbh = DBI->connect("dbi:SQLite:dbname=$kegg_db","","");
}

$newname||=$choosen_analysis;
$choosen_analysis||='IPR';

open(TAB, "<", $tabfile) or $LOGGER->logdie($!);
my $header_line=<TAB>;
chomp($header_line);
my @header=split(/\t/, $header_line);

{
    my %check;
    @check{ @header } = (1..$#header);

    $LOGGER->logdie("The column name already exists in tabular file") if (exists $check{ $newname });
    $LOGGER->logdie("The column name already exists in tabular file") if (exists $check{ $newname.'.id' });
    $LOGGER->logdie("The column name already exists in tabular file") if (exists $check{ $newname.'.desc' });

    $LOGGER->logdie("Not found column $colname in tabular file") unless (exists $check{ $colname });
}

my @newheader;
foreach my $h (@header) {
    push(@newheader, $h);
    if ( $h eq $colname ) {
        push(@newheader, $newname.'.id');
        push(@newheader, $newname.'.desc');
    }
}

my %annot;

my @res;
while(<TAB>) {
    chomp;
    my %data;
    @data{ @header } = split(/\t/, $_);
    push(@res, \%data);
    $annot{ $data{ $colname } } = undef;

    if ($cross_clstr) {
        my ($q) =  $data{$colname}=~/(cluster_\d+)/;
        if ($q) {
            if ($clstrmember{ $q }) {
                for my $x (@{$clstrmember{ $q }}) {
                    $annot{ $x } = undef;
                }
            } else {
                $LOGGER->logdie("Error: $colname ($data{$colname})");
            }
        }
    }

}
close(TAB);


open(IPR, "<", $iprfile) or $LOGGER->logdie($!);
while(<IPR>) {
    chomp;
    my ($acc, $md5, $seqlen, $analysis_name, $analysis_id, $analysis_desc, $start, $stop, $score, $status, $date, $ipr_id, $ipr_desc, $gos, $paths) = split(/\t/, $_);
    if (exists $annot{ $acc }) {
        if ($choosen_analysis eq 'IPR') {
            if ($ipr_id) {
                $annot{ $acc }->{ $ipr_id } = $ipr_desc;
            }                
        } elsif ($analysis_name eq $choosen_analysis) {
            $annot{ $acc }->{ $analysis_id } = $analysis_desc;
        } elsif ($choosen_analysis eq 'GO') {
            if ($gos) {
                foreach my $go_id ( split(/\|/, $gos) ) {
                   $annot{ $acc }->{ $go_id } = undef;
                }
            }
        } elsif ($choosen_analysis eq 'KEGG') {
            if ($paths) {
                foreach my $path_id ( split(/\|/, $paths) ) {
                   if ($path_id =~/^KEGG: ([^\+]+)(\+\S+)?/) {
                       my $kegg_id='ko'.$1;
                       if ($dbh) {
                           ($annot{ $acc }->{ $kegg_id }) = $dbh->selectrow_array("SELECT name FROM Pathways WHERE pid=?", undef, $kegg_id);
                       } else {
                            $annot{ $acc }->{ $kegg_id }='';
                       }
                   }
                }
            }
        }
    }
}
close(IPR);

print join("\t", @newheader),"\n";
foreach my $hr_data ( @res ) {
    
    my %info;

    my @x = $hr_data->{ $colname };
    if ($cross_clstr) {
        my ($q) =  $hr_data->{$colname}=~/(cluster_\d+)/;
        if ($q) {
            if ($clstrmember{ $q }) {
                push(@x, @{$clstrmember{ $q }});
            } else {
                $LOGGER->logdie("Error: $colname ($hr_data->{$colname})");
            }
        }
    }
    foreach my $q (@x) {
        foreach my $id (keys %{ $annot{ $q } }) {
            $info{$id} = $annot{$q}->{$id};
        }
    }
    
    $hr_data->{ $newname.'.id' } = join(';', map {$_||''} keys %info);
    $hr_data->{ $newname.'.desc' } = join(';', map {$_||''} values %info);
    
    print join("\t", @{$hr_data}{ @newheader }),"\n";
}

# Subroutines

sub Usage {
    my ($msg) = @_;
	my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -i      --ipr           InterProScan .tsv file
        -t      --tab           Tabular file (count matrix)
        -c      --col           ID column name
        -n      --name          Annotation column basename [Default: choosen analysis -a/--analysis]
        -a      --analysis      InterProScan Analysis [Default: IPR]
        -k      --kegg_db       KEGG db (ko.db from KOBAS)
        -x      --cross_clstr   Cross clstr data

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

