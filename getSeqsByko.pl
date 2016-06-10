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

my ($level,$koid);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "k|ko=s"=>\$koid
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


$LOGGER->logdie("Missing ko id") unless ($koid);

# Create a user agent object
use LWP::UserAgent;
my $ua = LWP::UserAgent->new;
$ua->agent("MyApp/0.1 ");

# Create a request
my $req = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$koid);
$req->content_type('application/x-www-form-urlencoded');
$req->content('query=libwww-perl&mode=dist');

# Pass request to the user agent and get a response back
my $res = $ua->request($req);

# Check the outcome of the response
if ($res->is_success()) {
    my $set = undef;
    my $content = $res->content();
    foreach my $line ( split(/\n/, $content) ) {
        if (($line =~ /^ORTHOLOGY/)||($set)) {
            last if (($set)&&($line=~/^\S+/));
            $set = 1;
            my ($kid, $kdesc) = $line=~/^(?:ORTHOLOGY)?\s+(\S+)\s+(.+)/;
            #print $kid,";",$kdesc,";\n";

            my $req_k = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$kid);
            $req_k->content_type('application/x-www-form-urlencoded');
            $req_k->content('query=libwww-perl&mode=dist');

            my $res_k = $ua->request($req_k);

            if ($res_k->is_success()) {
                my $set_k = undef;
                my $content_k = $res_k->content();
                foreach my $line_k ( split(/\n/, $content_k) ) {
                    if (($line_k =~ /^GENES/)||($set_k)) {
                        last if (($set_k)&&($line_k=~/^\S+/));
                        $set_k = 1;
                        my ($db, $acc_list) = $line_k=~/^(?:GENES)?\s+(\S+):\s+(.+)/;
                        #print "\t$db\n";
                        foreach my $acc ( split(/ /, $acc_list) ) {
                            $acc=~s/\([^\)]+\)//;
                            #print "\t\t$acc\n";

                            #http://rest.kegg.jp/get/ag:BAA36421/aaseq
                            
                            my $req_s = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.lc($db).':'.$acc.'/aaseq');
                            $req_s->content_type('application/x-www-form-urlencoded');
                            $req_s->content('query=libwww-perl&mode=dist');
                            
                            my $res_s= $ua->request($req_s);
                            
                            if ($res_s->is_success()) {
                                my $set_s = undef;
                                my $content_s = $res_s->content();
                                
                                my ($fasta_id, $fasta_desc, $fasta_seq) = (undef, undef, '');
                                foreach my $line_s ( split(/\n/, $content_s) ) {
                                    if ($line_s=~/^>(\S+)\s+(.*)/) {
                                        ($fasta_id, $fasta_desc) = ($1, $2);
                                    } else {
                                        $fasta_seq.=$line_s;
                                    }
                                }
                                
                                $fasta_seq=~s/(.{50})/$1\n/g;
                                print '>'.$kid.':'.$acc."\t".$fasta_desc."\n".$fasta_seq."\n";
                                
                            }

                        }
                    }
                }
            }
        }
    }
}
else {
    print $res->status_line, "\n";
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

        -h      --help  Help
        -l      --level Log level [Default: FATAL]
        -k      --ko    KEGG pathway (ko)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

