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
#  Copyright (C) 2020  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
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

Copyright (c) 2020 Universidade Estadual Paulista "Júlio de Mesquita Filho"

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Getopt::Long;


use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}


my (	$level, 
	$r1file, 
	$r2file,
	$nbases,
	$outfile 
	);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions(	"h|?|help" => sub { &Usage(); },
		"l|level=s"=> \$level,
		"r1|r1file=s"=>\$r1file,
		"r2|r2file=s"=>\$r2file,
		"o|outfile=s"=>\$outfile,
		"nb|nbases=i"=>\$nbases
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

$LOGGER->logdie("Missing R1 and/or R2 input file(s)") if (!(defined $r1file)&&(defined $r2file));
$LOGGER->logdie("Wrong R1 input file ($r1file)") if (! -e $r1file);
$LOGGER->logdie("Wrong R2 input file ($r2file)") if (! -e $r2file);

$nbases||=0;

$LOGGER->logdie("Missing output file") if (! defined $outfile);

`paste $r1file $r2file > $outfile.tmp`;

my $N='N' x $nbases;
my $Q='!' x $nbases;

if ( -e "$outfile.tmp" ) {
	open(IN, "<", "$outfile.tmp") or $LOGGER->logdie($!);
    open(OUT, ">", "$outfile") or $LOGGER->logdie($!);
    while(<IN>) {
        chomp;
        if ( ($. % 4) == 1) {
            my @twoseqs;
            while ($_=~/\@(\S+)/g) {
                push(@twoseqs, $1);
            }
            if ($twoseqs[0] ne $twoseqs[1]) {
                $LOGGER->logdie("Seq R1 ($twoseqs[0]) differs of ($twoseqs[1])");
            } else {
                print OUT '@'.$twoseqs[0]."\n";
            }
        } elsif ( ($. % 4) == 2) {
            $_=~s/\s+/$N/;
            print OUT $_,"\n",
        } elsif ( ($. % 4) == 3) {
            print OUT '+',"\n";
        } elsif ( ($. % 4) == 0) {
            $_=~s/\s+/$Q/;
            print OUT $_,"\n",
        }
	}
    close(IN);
    close(OUT);
    unlink("$outfile.tmp");
}

# Subroutines

sub Usage {
    my ($msg) = @_;

    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2020 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

	-h	--help			Help
	-l	--level			Log level [Default: FATAL]
	-r1	--r1file		Read 1 (forward) file
	-r2	--r2file		Read 2 (reverse) file
	-o	--outfile		Output file
    -n  --nbases        Number of Ns between reads

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

