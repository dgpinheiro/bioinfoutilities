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
use Storable;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $freqthreshold, $hetthreshold, $ord, $poplvl, $restlvl);

my %order = ('f'=>'frequency', 'd'=>'data');
my %alllvls = ('one'=>1, 'two'=>2, 'three'=>3);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=> \$infile,
            "t|freqthreshold=s"=>\$freqthreshold,
            "h|hetthreshold=s"=>\$hetthreshold,
            "e|level=s"=>\$poplvl,
            "o|order=s"=>\$ord,
            "r|restlevel=s"=>\$restlvl
    ) or &Usage();

$ord||='d';
$poplvl||='one';


if (defined $hetthreshold) {
    $LOGGER->logdie("Heterozygous threshold must be greather than 1") if ($hetthreshold < 2);
}
if (defined $restlvl) {
    $LOGGER->logdie("Wrong restrict level ($restlvl). Choose one of: ".join(', ', keys %alllvls)) unless (exists $alllvls{$restlvl});
}

$LOGGER->logdie("Wrong order ($ord). Possible values are: ".join(',', map { $_.' ('.$order{$_}.')' } keys %order ).".") unless (exists $order{$ord});    

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

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

#######################
# FIRST STEP: Load data and split them by number of heterozygous
#######################
$LOGGER->info("FIRST STEP...");

my %DATA;

#open(POP, "<", 'pop.txt') or $LOGGER->logdie($!);
#while(<POP>) {
while(<DATA>){
    chomp;
    next if ($_=~/^#/);
    my ($lv_one, $lv_two, $hmap, $desc) = split(/\t/, $_);
    $DATA{'POPULATION'}->{$hmap} = {    'one'   => $lv_one,
                                        'two'   => $lv_two,
                                        'three' => $hmap,
                                        'desc'  => $desc
                                    };
}

#close(POP);

open(IN, "<", $infile) or $LOGGER->logdie($!);

my %qtd;

my @header;
my @rs;

my %haplotype;

#my @full;

my %view;

my %subject;
## TEST DATA FOR FIRST STEP

#my @in= ( "hmap\tpaciente\tsex\trs1\trs2\trs3\trs4\n", 
#        "P1\tS1\tMale\tA/A\tC/C\tC/C\tC/N\n",
#        "P1\tS2\tMale\tA/T\tG/G\tC/C\tG/C\n",
#        "P1\tS3\tMale\tA/A\tG/G\tC/C\tG/G\n",
#        "P1\tS4\tFemale\tN/A\tN/N\tC/C\tG/G\n",
#        "P2\tS5\tFemale\tT/G\tC/C\tC/C\tG/G\n",
#        "P2\tS6\tFemale\tN/T\tC/C\tC/N\tG/G\n",
#        "P2\tS7\tFemale\tA/A\tC/C\tC/C\tG/G\n",
#        "P2\tS8\tFemale\tG/A\tG/C\tG/C\tG/G\n",
#        "P3\tS9\tMale\tA/A\tC/C\tT/C\tG/C\n"
#        );
#my $x = 1;
#foreach(@in) {

while(<IN>) {
    chomp;
    if ($.==1) {
#    if ($x==1) {
#        $x++;
        @header = split(/\t/, $_);
        foreach (@header) {
            push(@rs, $_) if ($_=~/^rs\d+/);
        }
#        print join("\t", 'LV1', 'LV2', 'LV3', 'ID', @rs, 'HAPLOTYPE'),"\n";
#        print join("\t", 'LV1', 'LV2', 'LV3', 'ID', 'HAPLOTYPE'),"\n";
    }
    else {
        my %v;
        @v{@header} = split(/\t/, $_);
        
        @rs = @rs;
        for my $m (@rs) {
            $v{$m} = [ sort {$a cmp $b} split('/', $v{$m} ) ];
        }

        $v{'sex'} = 'Unknow' if (($v{'sex'} ne 'Male')&&($v{'sex'} ne 'Female'));
        
        $subject{$v{'paciente'}} = \%v;
        
        my $hap = join(';', (map { (($_->[0] eq $_->[1]) ? $_->[0] : $_->[0].'/'.$_->[1]) } @v{@rs}));
        $qtd{ 'one' }->{ $DATA{'POPULATION'}->{ $v{'hmap'} }->{'one'} }++;
        $qtd{ 'two' }->{ $DATA{'POPULATION'}->{ $v{'hmap'} }->{'two'} }++;
        $qtd{ 'three' }->{ $DATA{'POPULATION'}->{ $v{'hmap'} }->{'three'} }++;
#        print join("\t",    $DATA{'POPULATION'}->{ $v{'hmap'} }->{'one'}, 
#                            $DATA{'POPULATION'}->{ $v{'hmap'} }->{'two'},
#                            $v{'hmap'},
#                            $v{'paciente'},
#                            $hap
#                        ),"\n";
        

#        push(@full,     [   $DATA{'POPULATION'}->{ $v{'hmap'} }->{'one'},
#                            $DATA{'POPULATION'}->{ $v{'hmap'} }->{'two'},
#                            $v{'hmap'},
#                            $v{'paciente'},
#                            $hap
#                        ]);

         my $nhet=0;
         $nhet++ while($hap=~m/\//g);

         push(@{ $view{ $DATA{'POPULATION'}->{ $v{'hmap'} }->{$poplvl} }->{ $nhet } }, [$hap, $v{'paciente'}]);


#        print join("\t",    $DATA{'POPULATION'}->{ $v{'hmap'} }->{'one'}, 
#                            $DATA{'POPULATION'}->{ $v{'hmap'} }->{'two'},
#                            $v{'hmap'},
#                            $v{'paciente'},
#                            (map { $_->[0] } @v{@rs}),
#                            join('', (map { $_->[0] } @v{@rs}))
#                        ),"\n";
#
#        print join("\t",    $DATA{'POPULATION'}->{ $v{'hmap'} }->{'one'}, 
#                            $DATA{'POPULATION'}->{ $v{'hmap'} }->{'two'},
#                            $v{'hmap'},
#                            $v{'paciente'},
#                            (map { $_->[1] } @v{@rs}),
#                            join('', (map { $_->[1] } @v{@rs}))
#                        ),"\n";
    }
}

close(IN);

my %general;
my %initial;

foreach my $pel (keys %view) {

    my %freq;

    #######################
    # SECOND STEP: Count frequency for 0 and 1 heterozygous in the haplotype
    #######################
    $LOGGER->info("Population: $pel # SECOND STEP...");

    foreach my $vid (0..1) {
        
        $LOGGER->info("Heterozygosity number: $vid");

        foreach my $view_ar ( @{ $view{$pel}->{$vid} } ) {

            my ($hap, $sbjct) = @{ $view_ar };
            if ($vid == 0) {
                my $f = 2;
                $freq{$hap}->{'data'}->{$sbjct} = $f;
                $freq{$hap}->{'f'}+=$f;
                $LOGGER->debug("$pel:$sbjct:$hap=$f");
            }
            else {
                my $orig = $hap;
                my @alter;
                &recpossibleshet($hap, \@alter, 0);
                $LOGGER->logdie("Wrong number of alleles: ".scalar(@alter)) if (scalar(@alter)!=2);
                foreach my $alt (@alter) {
                    my $f = 2/scalar(@alter);
                    $freq{$alt}->{'data'}->{$sbjct} = $f;
                    $freq{$alt}->{'f'}+=$f;
                    $LOGGER->debug("$pel:$sbjct:$alt=$f");
                }

            }
        }
    }    
    
    ## TEST DATA FOR THIRD STEP
    #%freq = ( 'A;C;G;N;T;N' =>  {   'data'  =>  [ ['S1',2], ['S2',2] ] ,
    #                                'f'     =>  4
    #                            },
    #          'A;C;G;T;T;T' =>  {   'data'  =>  [ ['S3',2] ],
    #                                'f'     =>  2
    #                            },
    #          'A;C;G;T;T;N' =>  {   'data'  =>  [ ['S4',2] ],
    #                                'f'     =>  2
    #                            },
    #          'A;A;T;A;A;A' =>  {   'data'  =>  [ ['S5',2] ],
    #                                'f'     =>  2
    #                            }           
    #);

    #######################
    # THIRD STEP: Count Ns
    #######################
    $LOGGER->info("Population: $pel # THIRD STEP...");

    my @tmp = keys %freq;
    foreach my $h ( @tmp ) {
        if ($h =~/N/) {
            my $orig = $h;
            #print ">",$orig,"\t",scalar(@{$freq{$orig}}),"\n";
            my @alter;
            &recpossibles($h, \@alter, 0);
            my %aux;
            my $qglobal = 0;

            # SE EXISTE GANHA UMA PROBABILIDADE MAIOR

            foreach my $alt (@alter) {
                if ( exists $freq{$alt} ) {
                    $aux{$alt} =  1+$freq{$alt}->{'f'};
                    $qglobal+=$aux{$alt};
                }
                else {
                    $aux{$alt}=1;
                    $qglobal+=1;
                }
            }
            
            
            foreach my $alt (@alter) {
                #print "\t",$alt,"\n";
                foreach my $sbjct (keys %{ $freq{ $orig }->{'data'} }) {
                    my $f = ($freq{ $orig }->{'data'}->{ $sbjct }*$aux{$alt})/$qglobal;
                    $freq{$alt}->{'data'}->{$sbjct} = $f;
                    $freq{$alt}->{'f'}+=$f;
                    $LOGGER->debug("$pel:$sbjct:$alt=$f");
                }                
            }
            delete($freq{$orig});
        }
    }
    
    foreach my $h ( keys %freq ) {
        foreach my $sbjct(keys %{ $freq{ $h }->{'data'} } ) {
            $initial{ $h }->{'data'}->{ $sbjct }=$freq{ $h }->{'data'}->{ $sbjct };
            $initial{ $h }->{'f'}+=$initial{ $h }->{'data'}->{ $sbjct };
        }
    }

}

%general = %initial;

foreach my $pel (keys %view) {
    
    my %freq = %initial;

    #######################
    # FOURTH STEP: Count frequency for 2 or more heterozygous in the haplotype
    #######################
    $LOGGER->info("Population: $pel # FOURTH STEP...");

    foreach my $vid (2..($hetthreshold||scalar(@rs))) {
        $LOGGER->info("Heterozygosity number: $vid");

        foreach my $view_ar ( @{ $view{$pel}->{$vid} } ) {
            my ($hap, $sbjct) = @{ $view_ar };
            my @haps = split(/\;/, $hap);

            my %tmp;
            foreach my $h (keys %freq) {
                my @hs = split(/\;/, $h);
                my $match = 0;
                for (my $i=0; $i<=$#hs; $i++) {
                    if ($haps[$i]=~/$hs[$i]|N/) {
                        $match++;
                    }
                    else {
                        last;
                    }
                }
                if ($match==scalar(@hs)) {
                    $tmp{ $h } = undef;
                }
            }
            my @alter = keys %tmp;

            my @betteralter;
            if (scalar(@alter) == 0 ) {
                $LOGGER->logwarn("Not found an aproximate haplotype for $hap. Trying to found all possibles...");
                &recpossibleshet($hap,\@betteralter,0);
            }
            else {
                my $c = 0;
                my $last_freq = 0;
                for my $x (sort { (($freq{$b}) ? $freq{$b}->{'f'} : 0) <=> (($freq{$a}) ? $freq{$a}->{'f'} : 0) } @alter) {
                    if ($c<1) {
                        push(@betteralter,$x);
                        my $f = $freq{$x}->{'f'}||0;
                        if ($last_freq<$f) {
                            $c++;
                        }
                        $last_freq = $f;
                    }
                }
                $betteralter[1] = join(';', &otherhap(\@haps,$betteralter[0]));
            }

                foreach my $balt (@betteralter) {
                    my $f = 2/scalar(@betteralter);
#                    print ">>>>>>>>>>>>>>>>>>$f (".scalar(@betteralter).")\n";
                    $freq{$balt}->{'data'}->{$sbjct} = $f;
                    $freq{$balt}->{'f'}+=$f;
                    $LOGGER->debug("$pel:$sbjct:$balt=$f");
                }
        } 
    }

    #######################
    # FIFTH STEP: Merge haplotypes and frequencies
    #######################
    $LOGGER->info("Population: $pel # FIFTH STEP...");
    
    foreach my $h ( keys %freq ) {
        foreach my $sbjct(keys %{ $freq{ $h }->{'data'} } ) {
            my $aux = ($freq{ $h }->{'data'}->{ $sbjct }-($initial{ $h }->{'data'}->{$sbjct}||0));
            $general{ $h }->{'data'}->{ $sbjct }+=$aux;
            $general{ $h }->{'f'}+=$aux;
        }
    }
    
}

# FINAL RESULT
my @arres;
if ($ord eq 'f') {
    @arres = (sort { $general{$b}->{'f'} <=> $general{$a}->{'f'} } keys %general);
}    
else {    
    @arres = (sort {$a cmp $b} keys %general);
}

my @sellvls;
if (defined $restlvl) {
    @sellvls = $restlvl;
}
else {
    @sellvls = sort { $alllvls{$a} <=> $alllvls{$b} } keys %alllvls;
}

print join(';', @rs),"\n";
foreach my $h ( @arres ) {
    if ((! defined $freqthreshold) || ($general{$h}->{'f'} >= $freqthreshold)) {
        print $h,"\t",  sprintf("%.2f",$general{$h}->{'f'}).'/'.(scalar(keys %subject)*2),"\t",
                        sprintf("%.2f",100*( $general{$h}->{'f'}/(scalar(keys %subject)*2))),"%\n";
        # FREQUENCIES                        
        my %r;
        foreach my $sbjct (keys %{ $general{ $h }->{'data'} }) {
            foreach my $lv ( @sellvls ) {
                $r{$lv}->{'data'}->{ 
                         $DATA{'POPULATION'}->{ $subject{ $sbjct }->{'hmap'}  }->{$lv}
                         }+= $general{ $h }->{'data'}->{$sbjct};
                $r{$lv}->{'f'}+= $general{ $h }->{'data'}->{$sbjct};
            }
        }
        foreach my $lv ( @sellvls ) {
            print "\tLEVEL: $lv\t(",sprintf("%.2f", ($r{$lv}->{'f'})),"/",(scalar(keys %subject)*2)," = ",sprintf("%.2f", 100*($r{$lv}->{'f'}/(scalar(keys %subject)*2))),"%)\n";
            foreach my $el (keys %{ $r{$lv}->{'data'} }) {
                print "\t\t",$el,"\t",sprintf("%.2f",$r{$lv}->{'data'}->{$el}),"/",($qtd{$lv}->{$el}*2),"\t",sprintf("%.2f",100*($r{$lv}->{'data'}->{$el}/($qtd{$lv}->{$el}*2))),"%\n";
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

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -i      --infile        Input file
        -f      --freqthreshold Frequency threshold
        -h      --hetthreshold  Heterozygotes threshold
        -o      --order         Order 
        -e      --level         Population level [Default: one]
        -r      --restlevel     Restrict to this level 

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

#print "\n\n###\n\n";
#my @ar;
#&recpossibles('N;N;N;N;N',\@ar,0);
#foreach my $x (@ar) {
#    print $x,"\n";
#}
#
#print "\n\n###\n\n";
#my @ar;
#&recpossibleshet('A;T;C/A;G;T/A',\@ar,0);
#foreach my $x (@ar) {
#    print $x,"\n";
#}

sub recpossibles {
    my ($h,$arr,$lv) = @_;
    if ($h=~m/N/g) {
        my $p = pos($h);
#        print "$lv BEGIN($p): $h\n";
        my $tmp = $h;
        foreach my $b ('A','C','G','T') {
#            print "$lv SUB($b): $h\n";
           substr($tmp, $p-1, 1, $b);
#            print "$lv SENT: \t\t\t",$h,"\n";
            &recpossibles($tmp, $arr, $lv+1);
        }
    }
    else {
#    print "$lv STORED: $h\n";
        push(@{$arr}, $h);
    }
}

sub recpossibleshet {
    my ($h,$arr,$lv,$hr_freq) = @_;
#    print "HAP: $h\n";
    if ($h=~m/\//g) {
#        print "FOUND HET: $h\n";
        my @allm = split(/\;/, $h);
        #print join('-', @allm),"\n";
        my $i = 0;
        while ( $allm[$i]!~/\// ) {
#            print $i,"\t",$allm[$i],"\n";
            $i++;
        }
        my ($left, $right) = split(/\//, $allm[$i]);
        my $tmp;
        #print ">>>$i\t$allm[$i]\n";
        $allm[$i] = $left;
        $tmp=join(';', @allm);
        &recpossibleshet($tmp, $arr, $lv+1, $hr_freq);
        $allm[$i] = $right;
        $tmp=join(';', @allm);
        &recpossibleshet($tmp, $arr, $lv+1, $hr_freq);
    }
    else {
#    print "$lv STORED: $h\n";
        if ((!defined $hr_freq) || (exists $hr_freq->{$h})) {
             push(@{$arr}, $h);
        }
    }
}

sub otherhap {
    my ($ar_haps, $h) = @_;
    my @r;
    my @h = split(';', $h);
    for (my $i=0;$i<=$#{$ar_haps};$i++) {
        if ($ar_haps->[$i] eq $h[$i]) {
            push(@r, $ar_haps->[$i]);
        }
        else {
            foreach my $eb (split(/\//, $ar_haps->[$i])) {
                unless ($eb eq  $h[$i]) {
                    push(@r, $eb);
                    last;
                }
            }
        }
    }
    return @r;
}


__DATA__
European	European	CEU	Utah residents (CEPH) with Northern and Western European ancestry (CEU)
European	European	TSI	Toscani in Italia (TSI)
European	European	GBR	British from England and Scotland (GBR)
European	European	FIN	Finnish from Finland (FIN)
European	European	IBS	Iberian populations in Spain (IBS)
Asian	East-Asian	CHB	Han Chinese in Beijing, China (CHB)
Asian	East-Asian	JPT	Japanese in Toyko, Japan (JPT)
Asian	East-Asian	CHS	Han Chinese South (CHS)
Asian	East-Asian	CHD	Chinese in Denver, Colorado (CHD) (pilot 3 only)
Asian	South-Asian	GIH	Gujarati Indian in Houston, TX (GIH)
African	West-African	YRI	Yoruba in Ibadan, Nigeria (YRI)
African	West-African	LWK	Luhya in Webuye, Kenya (LWK)
African	West-African	MKK	Maasai in Kinyawa, Kenya (MKK)
African	Afro-American	ASW	African Ancestry in Southwest US (ASW)
American	North-American	MEX	Mexican Ancestry in Los Angeles, CA (MXL)
American	North-American	PUR	Puerto Rican in Puerto Rico (PUR)
American	South-American	CLM	Colombian in Medellin, Colombia (CLM)
#Asian	East-Asian	CDX	Chinese Dai in Xishuangbanna (CDX)
#Asian	East-Asian	KHV	Kinh in Ho Chi Minh City, Vietnam (KHV)
#Asian	South-Asian	PJL	Punjabi in Lahore, Pakistan (PJL)
#Asian	South-Asian	BEB	Bengali in Bangladesh (BEB)
#Asian	South-Asian	STU	Sri Lankan Tamil in the UK (STU)
#Asian	South-Asian	ITU	Indian Telegu in the UK (ITU)
#African	West-African	GWD	Gambian in Western Division, The Gambia (GWD)
#African	West-African	MSL	Mende in Sierra Leono (MSL)
#African	West-African	ESN	Esan in Nigeria (ESN)
#African	Afro-American	ACB	African Caribbean in Barbados (ACB)
#American	South-American	PEL	Peruvian in Lima, Peru (PEL)
