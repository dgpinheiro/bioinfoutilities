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


use vars qw/$LOGGER $MERGED_DMP_FILE_PATH/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
    
    use File::Basename;
    use Cwd 'abs_path';

    $MERGED_DMP_FILE_PATH=dirname(abs_path(`readlink -f \$(which GetLCA.pl)`)).'/../taxonomy/merged.dmp';
}

my ($level, $mfile, $infile);

#Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "m|mfile=s"=>\$mfile,
            "i|infile=s"=>\$infile
    ) or &Usage();

use FileHandle;    
use POSIX 'isatty';
use Storable;

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

$mfile||=$MERGED_DMP_FILE_PATH;

$LOGGER->logdie("Not found Krona GetLCA.pl executable. Please, check if Krona GetLCA.pl path is on the PATH environment variable. Or put the value of NCBI taxonomy merged.dmp path file in --mfile/-m parameter.") unless ($mfile);

my $hr_data;

if (-e "$mfile.dump") {
    $hr_data = retrieve("$mfile.dump");
} else {
    open(MFILE, "<", $mfile) or $LOGGER->logdie($!);
    while(<MFILE>) {
        chomp;
        #12      |       74109   |
        #30      |       29      |
        $_=~s/\s+//g;
        my ($old, $new) = split(/\|/, $_);
        $hr_data->{ $old } = $new;
    }
    close(MFILE);
    store $hr_data, "$mfile.dump";
}

my $fh;

if ($infile) {

    $LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

    $fh = FileHandle->new;
    $fh->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fh = \*STDIN;
    } else {
        $LOGGER->logdie("Missing input file or STDIN data");
    }
}

while(<$fh>) {
    chomp;
    my $taxids=$_;
    $taxids=~s/^\s+//;
    $taxids=~s/^\s+$//;
    my @taxids=split(/\s+/, $taxids);
    foreach my $taxid (@taxids) {
        $taxid+=0;
        if (exists $hr_data->{$taxid}) {
            print $hr_data->{$taxid},"\n";
        } else {
            print "$taxid\n";
        }
    }        
}

$fh->close();

# Subroutines

sub Usage {
    my ($msg) = @_;

    my $mfile||=$MERGED_DMP_FILE_PATH;

	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file
        -m      --merged    merged.dmp file [Default: $mfile]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

