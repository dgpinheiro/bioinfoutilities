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

my ($level, $input);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|input=i"=> \$input
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

$LOGGER->logdie("Missing input JGI taxon_oid") unless ($input);
    
use JSON::Parse "parse_json"; 

#use WWW::Curl::Easy;

#my $curl = WWW::Curl::Easy->new;

#$curl->setopt(CURLOPT_HEADER, 0);
#$curl->setopt(CURLOPT_URL, 'https://taxonomy.jgi-psf.org/img/'.$input);


# A filehandle, reference to a scalar or reference to a typeglob can be used here.
#my $response_body;
#$curl->setopt(CURLOPT_WRITEDATA,\$response_body);

### Example: https://taxonomy.jgi.doe.gov/img/2724679250

use LWP::UserAgent;

my $ua = LWP::UserAgent->new;

$ua->ssl_opts(
    verify_hostname => 0
);

my $response = $ua->get('https://img.jgi.doe.gov/cgi-bin/genomesMetadata.cgi?ids='.$input);
#my $response = $ua->get('https://taxonomy.jgi-psf.org/img/'.$input);
 
my $response_body='';
if ($response->is_success) {
    $response_body=$response->decoded_content;
}
else {
    die $response->status_line;
}

my $perl = parse_json($response_body); 
        

print $perl->{$input}->{"ncbi_taxon_id"},"\n" if (($perl)&&($perl->{$input})&&($perl->{$input}->{"ncbi_taxon_id"}));


# Starts the actual request
#my $retcode = $curl->perform;
#
# Looking at the results...
#if ($retcode == 0) {
#       #print("Transfer went ok\n");
#	my $response_code = $curl->getinfo(CURLINFO_HTTP_CODE);
#        # judge result and next action based on $response_code
#	print("Received response: $response_body\n");
#
#	my $perl = parse_json($response_body); 
        
#	print $perl->{$input}->{"tax_id"},"\n" if (($perl)&&($perl->{$input})&&($perl->{$input}->{"tax_id"}));
#
#} else {
#        # Error code, type of error, error message
#	print("An error happened: $retcode ".$curl->strerror($retcode)." ".$curl->errbuf."\n");
#}




# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help  Help
        -l      --level Log level [Default: FATAL]
        -i      --input JGI IMG taxon_oid

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

