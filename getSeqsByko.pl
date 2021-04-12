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
use File::Temp qw/ tempfile tempdir /;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level,$koid,$dbdir);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "k|ko=s"=>\$koid,
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


$LOGGER->logdie("Missing ko id") unless ($koid);

my $dir;
if ($dbdir) {
    
    $LOGGER->logdie("Missing DB directory ($dbdir)") unless (-d $dbdir);

    $dir=$dbdir;

} else {

    $dir = tempdir('getSeqsByko_XXXXX', DIR => '.' , CLEANUP => 1 );

}


# Create a user agent object
use LWP::UserAgent;
my $ua = LWP::UserAgent->new;
$ua->agent("MyApp/0.1 ");

my @content_line;
if ( ( ! -e "$dir/ko/$koid.txt") || ( -z "$dir/ko/$koid.txt") ) {

    # Create a request
    my $req = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$koid);
    $req->content_type('application/x-www-form-urlencoded');
    $req->content('query=libwww-perl&mode=dist');

    # Pass request to the user agent and get a response back
    my $res = $ua->request($req);

    # Check the outcome of the response
    if ($res->is_success()) {
        my $content = $res->content();

        @content_line = split(/\n/, $content);
        
        mkdir("$dir/ko") unless (-d "$dir/ko");
        
        open(KO, ">", "$dir/ko/$koid.txt") or $LOGGER->logdie($!);
        print KO $content;
        close(KO);

    } else {
        print $res->status_line, "\n";
    }
    sleep(3);
} else {
    open(KO, "<", "$dir/ko/$koid.txt") or $LOGGER->logdie($!);
    while(<KO>) {
        chomp;
        push(@content_line, $_);
    }
    close(KO);
}   

my $set = undef;
foreach my $line ( @content_line ) {
    if (($line =~ /^ORTHOLOGY/)||($set)) {
        last if (($set)&&($line=~/^\S+/));
        $set = 1;
        my ($kid, $kdesc) = $line=~/^(?:ORTHOLOGY)?\s+(\S+)\s+(.+)/;
        #print $kid,";",$kdesc,";\n";
        
        my @content_k;

        if ( ( ! -e "$dir/KO/$kid.txt") || ( -z "$dir/KO/$kid.txt") ) {
            my $req_k = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$kid);
            $req_k->content_type('application/x-www-form-urlencoded');
            $req_k->content('query=libwww-perl&mode=dist');

            my $res_k = $ua->request($req_k);
            
            if ($res_k->is_success()) {
                my $content_k = $res_k->content();
                @content_k = split(/\n/, $content_k);

                mkdir("$dir/KO") unless (-d "$dir/KO");
                
                open(KO, ">", "$dir/KO/$koid.txt") or $LOGGER->logdie($!);
                print KO $content_k;
                close(KO);

            } else {
                print $res_k->status_line, "\n";
            }
        } else {
            open(KO, "<", "$dir/KO/$kid.txt") or $LOGGER->logdie($!);
            while(<KO>) {
                chomp;
                push(@content_k, $_);
            }
            close(KO);
        }
                
        my $set_k = undef;
        
        foreach my $line_k ( @content_k ) {
            if (($line_k =~ /^GENES/)||($set_k)) {
                last if (($set_k)&&($line_k=~/^\S+/));
                $set_k = 1;
                my ($db, $acc_list) = $line_k=~/^(?:GENES)?\s+(\S+):\s+(.+)/;
                #print "\t$db\n";
                foreach my $acc ( split(/ /, $acc_list) ) {
                    $acc=~s/\([^\)]+\)//;
                    #print "\t\t$acc\n";

                    #http://rest.kegg.jp/get/ag:BAA36421/aaseq
                    my $cat = lc($db);

                    my $set_s = undef;
                    my $content_s;

                    if ( ( ! -e "$dir/$cat/$acc.fa") || ( -z "$dir/$cat/$acc.fa") ) {
                    
                        my $req_s = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.lc($db).':'.$acc.'/aaseq');
                        $req_s->content_type('application/x-www-form-urlencoded');
                        $req_s->content('query=libwww-perl&mode=dist');
                        
                        my $res_s= $ua->request($req_s);
                        
                        if ($res_s->is_success()) {
                            $content_s = $res_s->content();
                        } else {
                            print $res_s->status_line(),"\n";
                        }                            
                    } else {
                        open(FASTA, "<", "$dir/$cat/$acc.fa") or $LOGGER->logdie($!);
                        while(<FASTA>) {
                            $content_s.=$_;
                        }
                        close(FASTA);
                    }

                    my ($fasta_id, $fasta_desc, $fasta_seq) = (undef, undef, '');
                    foreach my $line_s ( split(/\n/, $content_s) ) {
                        if ($line_s=~/^>(\S+)\s+(.*)/) {
                            ($fasta_id, $fasta_desc) = ($1, $2);
                        } else {
                            $fasta_seq.=$line_s;
                        }
                    }
                    
                    $fasta_seq=~s/(.{50})/$1\n/g;
                    print '>'.$kid.':'.lc($db).':'.$acc."\t".$fasta_desc."\n".$fasta_seq."\n";
                        
                }
            }
        }
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

        -h      --help  Help
        -l      --level Log level [Default: FATAL]
        -k      --ko    KEGG pathway (ko)
        -d      --db    Database directory of KEGG files

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

