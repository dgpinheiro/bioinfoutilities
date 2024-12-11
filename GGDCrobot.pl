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

my ($level, $qfile, $rfiles, $email, $blast, $threshold, $wait, $gfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "q|qfile=s"=>\$qfile,
            "r|rfiles=s"=>\$rfiles,
	    "g|gfile=s"=>\$gfile,
            "e|email=s"=>\$email,
            "b|blast=s"=>\$blast,
            "t|threshold=i"=>\$threshold,
            "w|wait=i"=>\$wait
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

$wait||=5;

$threshold||=50;

$LOGGER->logdie("Threshold must be in this range: 0 - 100") if (($threshold>100)||($threshold<0));

$LOGGER->logdie("Missing query file") unless ($qfile);
$LOGGER->logdie("Wrong query file ($qfile)") unless (-e $qfile);

my @multipleRefGenomes;

my @rfile;
my $gcontent;

if ($rfiles) {
	@rfile=split(/,/, $rfiles);

	foreach my $rf (@rfile) {
	    $LOGGER->logdie("Wrong reference file ($rf) in reference file list") unless (-e $rf);
	    my $rcontent='';
	    open(REF, "<", $rf);
	    while(<REF>) {
		$rcontent.=$_;
	    }
	    close(REF);

	    push(@multipleRefGenomes, {
		    filename => basename($rf),
		    content  => $rcontent,
		    content_type  => 'text/plain'
		}
	    );
	    
	}
} elsif ($gfile) {

	$gcontent='';
	open(REFACCS, "<", $gfile);
	while(<REFACCS>) {
	    if ($_=~/^\S+/) {
		    $gcontent.=$_;
	    }
	}
	close(REFACCS);
	
} else {
	$LOGGER->logdie("Missing at least one reference fasta file or one text file with Genbank accs");
}


$blast||='GBDP2_BLASTPLUS';

$LOGGER->logdie("Missing email") unless ($email);
$LOGGER->logdie("Wrong email ($email)") unless ($email=~/^[^\@ :;,\]\[\{\}\?\\\/\<\>\|\#]+\@[^\@ :;,\]\[\{\}\?\\\/\<\>\|\#]+$/);

my %avail;
@avail{'GBDP2_BLASTPLUS', 'GBDP2_BLAT', 'GBDP2_BLASTZ', 'GBDP2_WU-BLAST', 'GBDP2_MUMMER'} = ();

$LOGGER->logdie("Wrong program $blast. Please select one of these: GBDP2_BLASTPLUS,GBDP2_BLAT,GBDP2_BLASTZ,GBDP2_WU-BLAST,GBDP2_MUMMER") unless (exists $avail{$blast});

use strict;
use warnings;

use File::Basename;

use HTTP::Tiny;
use HTTP::Tiny::Multipart;

while (&try()>$threshold) {
    print STDERR "WAIT! Current slot usage is above $threshold ...\n";
    &waitmin($wait);
}

my $url = 'https://ggdc.dsmz.de/submit_ggdc_job.php';

my $http = HTTP::Tiny->new();
my $qcontent='';
open(QUERY, "<", $qfile);
while(<QUERY>) {
    $qcontent.=$_;
}
close(QUERY);

my $hr_input_params = {
    targetName => '',
    targetGenome => {
        filename => basename($qfile),
        content  => $qcontent,
        content_type  => 'text/plain'
        },
    'email' => $email,
    'blastVariant' => $blast
 };

if (scalar(@rfile)) {
	$LOGGER->info("Job: ".basename($qfile)." X ".scalar(@rfile)." genome files\n");
    	$hr_input_params->{'multipleRefGenomes[]'} = \@multipleRefGenomes;
} elsif ($gcontent) {
	my $nlines = () = $gcontent =~ /\n/g;
	# Add 1 if the last line doesn't end with a newline
	$nlines++ if $gcontent !~ /\n\z/;

	$LOGGER->info("Job: ".basename($qfile)." X ".$nlines." genome accessions\n");
    	$hr_input_params->{'refGenbank'} = $gcontent;
} else {
	$LOGGER->logdie("Job sent failed. Wrong input parameters");
}

my $response = $http->post_multipart( $url, $hr_input_params);

if (length $response->{content}) {
	if ($response->{content}=~/your job with ID <b>([^<]+)<\/b> has been successfully submitted/) {
		print "Your job with ID $1 has been successfully submitted and the results should be found at $email email box\n";
	} else {
		my $randnum=rand(10000);
		open(OUT,">", "/tmp/GGDCrobot_error_$randnum.html") or $LOGGER->logdie("$!");
		print OUT $response->{content};
		close(OUT);
		$LOGGER->logdie("Something wrong with this submission. See the content of /tmp/GGDCrobot_error_$randnum.html");
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
        -q      --qfile     Query fasta file
        -r      --rfiles    Reference fasta files separated by comma
	-g	--gfile	    Text file with Genbank references. One line per genome. If genome consists of multiple accessions provide these either separated by blanks (e.g. 'AE000782 AE000783 AE000784') or as a range (e.g. 'AE000782-AE000784')
        -e      --email     Email
        -b      --blast     blastVariant {GBDP2_BLASTPLUS,GBDP2_BLAT,GBDP2_BLASTZ,GBDP2_WU-BLAST,GBDP2_MUMMER} [Default: GBDP2_BLASTPLUS]
        -t      --threshold Threshold for current slot usage (wait if current slot usage > than this threshold to send a job) Range: 0-100 [Default: 50]
        -w      --wait      Wait time for check current slot usage in minutes [Default: 5]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub try {
    my $resp = HTTP::Tiny->new->get('https://ggdc.dsmz.de/ggdc.php');
     
    die "Failed!\n" unless $resp->{success};
     
    if (length $resp->{content}) {
        my @content=split(/\n/, $resp->{content});
        foreach my $l (@content) {
            if ($l=~/div class="progress-bar progress-bar-success" role="progressbar" aria-valuenow="3" aria-valuemin="0" aria-valuemax="100" style="width: (\d+)%;">/) {
                my $n = $1;
                return $n;
            }
        } 
    } else {
        die "Error. Cannot get content";       
    }
}

sub waitmin {
    my ($n) = $_;
    $n||=1;
    sleep($n*60);
}
