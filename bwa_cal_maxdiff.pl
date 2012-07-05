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
#  Copyright (C) 2010  Fundação Hemocentro de Ribeirão Preto
#
#  Laboratório de Genética Molecular e Bioinformática
#  Núcleo de Bioinformática
#  BiT -  Bioinformatics Team
#  Fundação Hemocentro de Ribeirão Preto
#  Rua Tenente Catão Roxo, 2501
#  Ribeirão Preto - São Paulo
#  Brasil
#  CEP 14051-140
#  Fone: 55 16 39639300 Ramal 9603
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://lgmb.fmrp.usp.br
#
# $Id$

=head1 NAME

    bwa_cal_maxdiff.pl - Calculate the k (edit distance) based on parameters.

=head1 SYNOPSIS
    
    # uniform error rate 2%
    # sequence length: 36bp
    # probability threshold: 4%

    $ bwa_cal_maxdiff.pl 36 0.02 0.04

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Calculate the k (edit distance) based on poisson distribution, given the uniform error rate, the threshold probability and the sequence length.

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2010 Regional Blood Center of Ribeirão Preto

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;

my $length = $ARGV[0];
die "Missing length" unless ($length);
my $err = $ARGV[1];
die "Missing error" unless ($err);
my $thres = $ARGV[2];
die "Missing threshold" unless ($thres);

print "R: ",&bwa_cal_maxdiff($length, $err, $thres),"\n";

sub bwa_cal_maxdiff {
    my ($l, $err, $thres) = @_;
    my $elambda = exp(-$l * $err);
    my $k = 1;
    my $x = 1;
    my $sum = $elambda;
    my $y = 1;
    print "k\ty\tx\tsum\t1-sum\n";
    for (;$k<1000;++$k) {

        $y*=$l*$err;
        $x*=$k;
        # sum = soma das probabilidades de que existam "exatamente" k ocorrências (erros) = "no máximo k" ocorrências (erros)
        # P(X<=k) = P(X=1)+...+P(X=k)
        $sum+=($elambda*$y/$x);
        print $k,"\t",$y,"\t",$x,"\t",$sum,"\t",(1-$sum),"\n";
        # 1-sum = probabilidade de que existam "no mínimo" k ocorrências (erros)
        # P(X>k) = 1-P(X<=k)
        if ((1-$sum) < $thres) {
            return $k;
        }
    }
    return 2;
}
