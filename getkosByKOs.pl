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
use File::Basename;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $KOid, $infile, $outprefix, $keep, $dbdir, $ignorehuman);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "K|KO=s"=>\$KOid,
            "i|infile=s"=>\$infile,
            "o|outprefix=s"=>\$outprefix,
            "p|keep"=>\$keep,
            "d|db=s"=>\$dbdir,
            "ih|ignore_human"=>\$ignorehuman
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

$outprefix||='./getkosByKOs';

if ($dbdir) {
    $keep=1;
}

my %KO;

unless ($KOid) {
    unless ($infile) {
        $LOGGER->logdie("Missing KO info: id (-K) or input file (-i).");
    } else {
        open(IN, "<", $infile) or $LOGGER->logdie($!);
        while(<IN>) {
            chomp;
            my ($id) = $_;
            $KO{uc($id)} = undef;
        }
        close(IN);
    }
} else {
    $KO{uc($KOid)} = undef;
}


# Create a user agent object
use LWP::UserAgent;
my $ua = LWP::UserAgent->new;
$ua->agent("MyApp/0.1 ");

my $dir;
if ($keep) {
    if ($dbdir) {
       $dir=$dbdir;
    } else {
       $dir = dirname($outprefix).'/'.basename($outprefix);
    }
    mkdir($dir);
} else {
    $dir = tempdir(basename($outprefix).'_XXXXX', DIR => dirname($outprefix), CLEANUP => (($keep) ? 0 : 1) );
}

foreach my $kid (keys %KO) {

    my @content_line;
    
    if ( ( ! -e "$dir/KO/$kid.txt") || ( -z "$dir/KO/$kid.txt") ) {

        # Create a request
        my $req = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/ko:'.$kid);
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
            sleep(3);
        } else {
            print $kid,"\t",$res->status_line, "\n";
            next;
        }
    } else {
        open(KO, "<", "$dir/KO/$kid.txt") or $LOGGER->logdie($!);
        while(<KO>) {
            chomp;
            push(@content_line, $_);
        }
        close(KO);
    }
        
    if (scalar(@content_line)) {
        my $set_name = undef;
        my $set_definition = undef;
        my $set_pathway = undef;
        foreach my $line ( @content_line ) {
            if (($line =~ /^NAME/)||($set_name)) {
                if (($set_name)&&($line=~/^[A-Z_]+\s+/)&&($line!~/^NAME/)) {
                    $set_name=undef;
                } else {
                    $set_name = 1;
                    my ($kname) = $line=~/^(?:NAME\s+)?(.+)/;
                    $KO{$kid}->{'name'}.=$kname if ($kname);
                }
            }
            if (($line =~ /^DEFINITION/)||($set_definition)) {
                if (($set_definition)&&($line=~/^[A-Z_]+\s+/)&&($line!~/^DEFINITION/)) {
                    $set_definition=undef;
                } else {
                    $set_definition = 1;
                    my ($kdesc) = $line=~/^(?:DEFINITION\s+)?(.+)/;
                    $KO{$kid}->{'definition'}.=$kdesc if ($kdesc);
                }
            }
            if (($line =~ /^(?:MODULE|PATHWAY)/)||($set_pathway)) {
                if (($set_pathway)&&($line=~/^[A-Z_]+\s+/)&&($line!~/^(?:MODULE|PATHWAY)/)) {
                    $set_pathway=undef;
                    last;
                } else {
                    $set_pathway = 1;
                    if ($line=~/\s+((?:ko|map)\d{5})\s+.*/) {
                        my $ko = $1;
                        $KO{$kid}->{'ko'}->{$ko} = undef;
                    } elsif ($line=~/\s+((?:M)\d{5})\s+.*/) {
                        # Not used
                        my $mo = $1;
                        $KO{$kid}->{'M'}->{$mo} = undef;
                    }
                }
            }
            
            
        }
    }
}

my %pathway;
open(OUT, ">", $outprefix.'-KO-annotations.txt') or $LOGGER->logdie($!);

foreach my $kid (keys %KO) {
    print OUT join("\t", $kid, map { $_||'' } @{ $KO{$kid} }{'name', 'definition'} ),"\n";
    
    foreach my $pid (keys %{ $KO{$kid}->{'ko'} }) {
        
        next if (exists $pathway{$pid});

        my @content_line;

        if ( ( ! -e "$dir/ko/$pid.txt") || ( -z "$dir/ko/$pid.txt") ) {

            # Create a request
            my $req = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$pid);
            $req->content_type('application/x-www-form-urlencoded');
            $req->content('query=libwww-perl&mode=dist');
            
            # Pass request to the user agent and get a response back
            my $res = $ua->request($req);
            # Check the outcome of the response
            if ($res->is_success()) {
                my $content = $res->content();
                @content_line = split(/\n/, $content);

                mkdir("$dir/ko") unless (-d "$dir/ko");

                open(KO, ">", "$dir/ko/$pid.txt") or $LOGGER->logdie($!);
                print KO $content;
                close(KO);
                sleep(3);
            } else {
                print $pid,"\t",$res->status_line, "\n";
                next;
            }
        } else {
            open(KO, "<", "$dir/ko/$pid.txt") or $LOGGER->logdie($!);
            while(<KO>) {
                chomp;
                push(@content_line, $_);
            }
            close(KO);
        }

        if (scalar(@content_line)) {
            my $set_name = undef;
            my $set_orthology = undef;
            my $set_module = undef;
            my $set_class = undef;
            foreach my $line ( @content_line ) {
                if (($line =~ /^NAME/)||($set_name)) {
                    if (($set_name)&&($line=~/^\S+/)&&($line!~/^NAME/)) {
                        $set_name=undef;
                    } else {
                        $set_name = 1;
                        my ($kname) = $line=~/^(?:NAME\s+)?(.+)/;
                        $pathway{$pid}->{'name'}.=$kname if ($kname);
                    }
                }
                if ($ignorehuman) {
                    if (($line =~ /^CLASS/)||($set_class)) {
                        if (($set_class)&&($line=~/^[A-Z_]+\s+/)&&($line!~/^CLASS/)) {
                            $set_class=undef;
                        } else {
                            $set_class = 1;
                            my ($kname) = $line=~/^(?:CLASS\s+)?(.+)/;
                            my $set_human_diseases=undef;
                            foreach my $c (split(/; /, $kname)) {
                                if ($c eq 'Human Diseases') {
                                    $set_human_diseases=1;
                                    delete($pathway{$pid});
                                    delete($KO{$kid}->{'ko'}->{$pid});
                                    last;
                                }
                            }
                            if ($set_human_diseases) {
                                last;
                            }
                        }
                    }
                } 
                if (($line =~ /^ORTHOLOGY/)||($set_orthology)) {
                    if (($set_orthology)&&($line=~/^[A-Z_]+\s+/)&&($line!~/^ORTHOLOGY/)) {
                        $set_orthology=undef;
                        last;
                    } else {
                        $set_orthology = 1;
                        if ($line=~/(K\d+)\s+(.*)/) {
                            my $K = $1;
                            my $name = $2;
                            $pathway{$pid}->{'orthology'}->{$K} = $name;
                        }
                    }
                }
                if (($line =~ /^MODULE/)||($set_module)) {
                    if (($set_module)&&($line=~/^[A-Z_]+\s+/)&&($line!~/^MODULE/)) {
                        $set_module=undef;
                        last;
                    } else {
                        $set_module = 1;
                        if ($line=~/(M\d+)\s+(.*)/) {
                            my $mid = $1;
                            my @m_content_line;
                            
                            if ( ( ! -e "$dir/M/$mid.txt") || ( -z "$dir/M/$mid.txt") ) {
                                # Create a request
                                my $m_req = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$mid);
                                $m_req->content_type('application/x-www-form-urlencoded');
                                $m_req->content('query=libwww-perl&mode=dist');
                                
                                # Pass request to the user agent and get a response back
                                my $m_res = $ua->request($m_req);
                                # Check the outcome of the response
                                if ($m_res->is_success()) {
                                    my $m_content = $m_res->content();
                                    @m_content_line = split(/\n/, $m_content);
                                    
                                    mkdir("$dir/M") unless (-d "$dir/M");
                                    
                                    open(MO, ">", "$dir/M/$mid.txt") or $LOGGER->logdie($!);
                                    print MO $m_content;
                                    close(MO);
                                    sleep(3);
                                } else {
                                    print $mid,"\t",$m_res->status_line, "\n";
                                    next;
                                }
                            } else {
                                    
                                open(MO, "<", "$dir/M/$mid.txt") or $LOGGER->logdie($!);
                                while(<MO>) {
                                    chomp;
                                    push(@m_content_line, $_);
                                }
                                close(MO);
                            }

                            if (scalar(@m_content_line)) {
                                my $set_m_name = undef;
                                my $set_m_orthology = undef;
                                foreach my $m_line ( @m_content_line ) {
                                    if (($line =~ /^NAME/)||($set_m_name)) {
                                        if (($set_m_name)&&($m_line=~/^\S+/)&&($m_line!~/^NAME/)) {
                                            $set_m_name=undef;
                                        } else {
                                            $set_m_name = 1;
                                            my ($mname) = $m_line=~/^(?:NAME\s+)?(.+)/;
                                        }
                                    }
                                    if (($m_line =~ /^ORTHOLOGY/)||($set_m_orthology)) {
                                        if (($set_m_orthology)&&($m_line=~/^[A-Z_]+\s+/)&&($m_line!~/^ORTHOLOGY/)) {
                                            $set_m_orthology=undef;
                                            last;
                                        } else {
                                            $set_m_orthology = 1;
                                            if ($m_line=~/(K\d+)\s+(.*)/) {
                                                my $K = $1;
                                                my $name = $2;
                                                $pathway{$pid}->{'orthology'}->{$K} = $name;
                                            }
                                        }
                                    }
                                }
                            }
                        }                            
                    }
                }
                
                
            }
        }
    
    }
    
}
close(OUT);

open(OUT, ">", $outprefix.'-ko-annotations.txt') or $LOGGER->logdie($!);
open(GMT, ">", $outprefix.'-ko-annotations.gmt') or $LOGGER->logdie($!);
foreach my $pid (keys %pathway) {
    print OUT join("\t", $pid, $pathway{$pid}->{'name'} ),"\n";
    print GMT $pid,"\t";
    if (($pathway{$pid}->{'orthology'})&&(scalar(keys %{ $pathway{$pid}->{'orthology'} }))) {
        print GMT join(";", map { $_||'' } keys %{ $pathway{$pid}->{'orthology'} } );
    }
    print GMT "\n";
}
close(OUT);
close(GMT);

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -K      --KO            KEGG Orthology ID (KO)
        -i      --infile        Input file
        -o      --outprefix     Output prefix
        -p      --keep          Keep temporary directory [True if -d has value]
        -d      --db            Database directory, i.e. with "KO/" (KEGG Orthology), "ko/" (KEGG Reference Pathway), and "M/" (KEGG Modules) directories
        -ih     --ignore_human  Ignore KEGG pathways specific from "Human Diseases" category

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

