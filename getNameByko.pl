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

use File::Temp qw/ tempfile tempdir /;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level,$koid,$infile,$dbdir);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "k|ko=s"=>\$koid,
            "i|infile=s"=>\$infile,
            "d|db=s"=>\$dbdir
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

my $dir;
if ($dbdir) {
    
    $LOGGER->logdie("Missing DB directory ($dbdir)") unless (-d $dbdir);

    $dir=$dbdir;

} else {

    $dir = tempdir('getNameByko_XXXXX', DIR => '.' , CLEANUP => 1 );

}

my %ko;

unless ($koid) {
    unless ($infile) {
        $LOGGER->logdie("Missing ko info: id (-k) or input file (-i).");
    } else {
        open(IN, "<", $infile) or $LOGGER->logdie($!);
        while(<IN>) {
            chomp;
            my ($id) = $_;
            $ko{lc($id)} = undef;
        }
        close(IN);
    }
} else {
    $ko{lc($koid)} = undef;
}


# Create a user agent object
use LWP::UserAgent;
my $ua = LWP::UserAgent->new;
$ua->agent("MyApp/0.1 ");


foreach my $kid (keys %ko) {
    my @content_line;
    if ( ( ! -e "$dir/KO/$kid.txt") || ( -z "$dir/KO/$kid.txt") ) {
        # Create a request
        my $req = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$kid);
        $req->content_type('application/x-www-form-urlencoded');
        $req->content('query=libwww-perl&mode=dist');
        
        # Pass request to the user agent and get a response back
        my $res = $ua->request($req);
        
        # Check the outcome of the response
        if ($res->is_success()) {
            my $content = $res->content();

            @content_line = split(/\n/, $content);

            mkdir("$dir/KO") unless (-d "$dir/KO");

            open(KO, ">", "$dir/KO/$kid.txt") or $LOGGER->logdie($!);
            print KO $content;
            close(KO);
        }
        else {
            print $res->status_line, "\n";
        }
        sleep(3);

    } else {
        open(KO, "<", "$dir/KO/$kid.txt") or $LOGGER->logdie($!);
        while(<KO>) {
            chomp;
            push(@content_line, $_);
        }
        close(KO);
    }        
     
    my $set = undef;
    foreach my $line ( @content_line ) {
        if ( ( $line =~ /^NAME/ ) || ($set) ) {
            last if ( ($set) && ( $line =~ /^\S+/ ) );
            $set = 1;
            my ($kdesc) = $line =~ /^(?:NAME)?\s+(.+)/;
            $ko{$kid} .= $kdesc;

        }
    }
}


foreach my $kid (keys %ko) {
    print $kid,"\t",($ko{$kid}||''),"\n";
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
        -k      --ko        KEGG pathway (ko)
        -i      --infile    Input file
        -d      --db        Database directory of KEGG files

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

