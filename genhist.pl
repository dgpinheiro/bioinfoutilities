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

use File::Spec;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $fqfile, $bamfilea, $gtfile, $genhisto, $div, $minreads, $reps, $gffile, $minrpkm);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "f|fqfile=s"=>\$fqfile,
            "a|bamfilea=s"=>\$bamfilea,
            "t|gtfile=s"=>\$gtfile,
            "r|gffile=s"=>\$gffile,
            "g|genhisto=s"=>\$genhisto,
            "m|minreads=i"=>\$minreads,
            "k|minrpkm=i"=>\$minrpkm,
            "d|division=i"=>\$div,
            "n|replicates=i"=>\$reps
    ) or &Usage();

$LOGGER->logdie("Missing fastq file") unless ($fqfile);
$LOGGER->logdie("Wrong fastq file ($fqfile)") unless (-e $fqfile);

$LOGGER->logdie("Missing bam file - alignments") unless ($bamfilea);
$LOGGER->logdie("Wrong bam file - alignments ($bamfilea)") unless (-e $bamfilea);

$LOGGER->logdie("Missing gtf file - cufflinks output") unless ($gtfile);
$LOGGER->logdie("Wrong gtf file - cufflinks output ($gtfile)") unless (-e $gtfile);

$LOGGER->logdie("Missing gff file - reference ") unless ($gffile);
$LOGGER->logdie("Wrong gff file - reference ($gffile)") unless (-e $gffile);

$div||=100000;
$minreads=1 unless (defined $minreads);
$minrpkm=0 unless (defined $minrpkm);
$reps||=1;
$genhisto||='.';

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

use File::Basename;

use Digest::MD5::File qw( file_md5_hex );

mkdir($genhisto) unless (-d $genhisto);

use R;
use Cwd;

use Storable;

use File::Basename;

&R::initR("--silent","--no-save");
&R::library("RSPerl");

$|=1;

my %checksum;
if (-e "$genhisto/checksum.txt") {
    open(CS, "<", "$genhisto/checksum.txt") or $LOGGER->logdie($!);
    while(<CS>) {
        chomp;
        my ($file, $checksum) = split(/\t/, $_);
        $checksum{$file} = $checksum;
    }
    close(CS); 
}

open(CS, ">", "$genhisto/checksum.txt") or $LOGGER->logdie($!);

my $accepted_hits_x_known_genes_md5='';

if (-e "$genhisto/accepted_hits_x_known_genes.txt") {
    $accepted_hits_x_known_genes_md5 = file_md5_hex("$genhisto/accepted_hits_x_known_genes.txt");
}

my %translenknown;
if ((! -e "$genhisto/accepted_hits_x_known_genes.txt")||
    ( ! exists $checksum{basename("$genhisto/accepted_hits_x_known_genes.txt")})||
    ( $checksum{basename("$genhisto/accepted_hits_x_known_genes.txt")} ne $accepted_hits_x_known_genes_md5 )) {
       
       $LOGGER->info("GFF to BED conversion ...");
       `grep 'exon' $gffile | gff2bed > $genhisto/reference_exons.bed`;
        
       open(EXON, "<", "$genhisto/reference_exons.bed") or $LOGGER->logdie($!);
       while(<EXON>){
           chomp;
           my ($chr, $start, $end, undef, undef, undef, undef, undef, undef, $attr) = split(/\t/, $_);
           my ($parent)= $attr=~/Parent=([^;]+)/;
           $translenknown{ $parent } = 0 unless (exists $translenknown{ $parent });
           $translenknown{ $parent }+=($end-$start);
       }
       close(EXON);
    
       store \%translenknown, "$genhisto/translength_known.dump";

       $LOGGER->info("intersectBed ...");
       `intersectBed -abam -bed -wo -a $bamfilea -b $genhisto/reference_exons.bed | nsort -k4,4 > $genhisto/accepted_hits_x_known_genes.txt`;
       `intersectBed -v -abam -bed -wo -a $bamfilea -b $genhisto/reference_exons.bed | nsort -k4,4 > $genhisto/remaining_accepted_hits.txt`;
       
       $accepted_hits_x_known_genes_md5 = file_md5_hex("$genhisto/accepted_hits_x_known_genes.txt");
       unlink("$genhisto/session.RData") if (-e "$genhisto/session.RData");
} else {
        $LOGGER->logdie("Not found $genhisto/translength_known.dump") unless (-e "$genhisto/translength_known.dump");
        $LOGGER->info("Retrieve length of known transcripts from dump ...");
        %translenknown = %{ retrieve("$genhisto/translength_known.dump") };
}

print CS basename("$genhisto/accepted_hits_x_known_genes.txt"), "\t", $accepted_hits_x_known_genes_md5, "\n";

my $accepted_hits_x_assembled_genes_md5='';

if (-e "$genhisto/accepted_hits_x_assembled_genes.txt") {
    $accepted_hits_x_assembled_genes_md5 = file_md5_hex("$genhisto/accepted_hits_x_assembled_genes.txt");
}
    
my %translenasm;
if ((! -e "$genhisto/accepted_hits_x_assembled_genes.txt")||
    ( ! exists $checksum{basename("$genhisto/accepted_hits_x_assembled_genes.txt")})||
    ( $checksum{basename("$genhisto/accepted_hits_x_assembled_genes.txt")} ne $accepted_hits_x_assembled_genes_md5 )) {
       
       $LOGGER->logdie("Not found $genhisto/remaining_accepted_hits.txt") unless (-e "$genhisto/remaining_accepted_hits.txt");
        
       #$LOGGER->info("GFF reference comparison (cuffcompare)");
       #`cuffcompare -r $gffile $gtfile -o $genhisto/transcripts`;
       $LOGGER->info("GTF to BED conversion ...");
       #`grep 'exon' $genhisto/transcripts.combined.gtf | gtf2bed > $genhisto/assembled_exons.bed`;
       `grep 'exon' $gtfile | gtf2bed > $genhisto/assembled_exons.bed`;
        
       open(EXON, "<", "$genhisto/assembled_exons.bed") or $LOGGER->logdie($!);
       while(<EXON>){
           chomp;
           my ($chr, $start, $end, undef, undef, undef, undef, undef, undef, $attr) = split(/\t/, $_);
           my ($parent)= $attr=~/transcript_id "([^"]+)"/;
           $translenasm{ $parent } = 0 unless (exists $translenasm{ $parent });
           $translenasm{ $parent }+=($end-$start);
       }
       close(EXON);

       store \%translenasm, "$genhisto/translength_asm.dump";

       $LOGGER->info("intersectBed ...");
       `intersectBed -bed -wao -a $genhisto/remaining_accepted_hits.txt -b $genhisto/assembled_exons.bed | nsort -k4,4 > $genhisto/accepted_hits_x_assembled_genes.txt`;
    
       $accepted_hits_x_assembled_genes_md5 = file_md5_hex("$genhisto/accepted_hits_x_assembled_genes.txt");
       unlink("$genhisto/session.RData") if (-e "$genhisto/session.RData");
} else {
        $LOGGER->logdie("Not found $genhisto/translength_asm.dump") unless (-e "$genhisto/translength_asm.dump");
        $LOGGER->info("Retrieve length of assembled transcripts from dump ...");
        %translenasm = %{ retrieve("$genhisto/translength_asm.dump") };
}

print CS basename("$genhisto/accepted_hits_x_assembled_genes.txt"), "\t", $accepted_hits_x_assembled_genes_md5, "\n";

my $genes_md5='';
if (-e "$genhisto/genes.txt") {
    $genes_md5 = file_md5_hex("$genhisto/genes.txt");
}

my %genetranscript;

my %genelen;
if ((!-e "$genhisto/genes.txt")||
    (! exists $checksum{basename("$genhisto/genes.txt")})||
    ( $checksum{basename("$genhisto/genes.txt")} ne $genes_md5 )){
    
    $LOGGER->info("Collecting genes's info from gff ...");
    #GroupUn4146     LBDAv0.4        mRNA    66      1123    .       -       .       ID=XM_001123196.2;Parent=LOC727487/GB48040
    my %genes;
    open(GFF, "<", $gffile) or $LOGGER->logdie($!);
    while(<GFF>) {
        next if ($_=~/^#/);
        chomp;
        my ($seq, $src, $type, $start, $end, $score, $strand, $phase, $extra)= split(/\t/, $_);
        if ($type =~/RNA|transcript/) {
            my ($id)=$extra=~/ID=([^;]+)/;
            my ($parent)=$extra=~/Parent=([^;]+)/;
            $genetranscript{$id} = $parent;
            $genes{$parent} = 1;
        }
    }
    close(GFF);

    $LOGGER->info("Collecting genes's info from gtf (identifying novel loci) ...");
    #Group1  Cufflinks       transcript      234561  234772  1000    .       .       gene_id "CUFF.2"; transcript_id "CUFF.2.1"; FPKM "9.2736794454"; frac "1.000000"; conf_lo "0.884430";  conf_hi "3.390314"; cov "14.837959"; full_read_support "yes";
    my %test;
    open(GTF, "<", $gtfile) or $LOGGER->logdie($!);
    while(<GTF>) {
        next if ($_=~/^#/);
        chomp;
        my ($seq, $src, $type, $start, $end, $score, $strand, $phase, $extra)= split(/\t/, $_);
        if ($type eq 'transcript') {
            my ($id)=$extra=~/transcript_id "([^"]+)"/;
            my ($parent)=$extra=~/gene_id "([^"]+)"/;
            $test{$parent}->{$id} = undef;
        }
    }
    close(GTF);
    # Há casos em que uma isoforma nova aparece como gene novo mas há isoformas com transcritos identificados
    # neste caso não é um novo loci
    my %similar;
    my %novel;
    foreach my $g (keys %test) {
        $novel{$g} = undef;
        foreach my $t (keys %{ $test{$g} }) {
            if ($t !~ /CUFF/) {
                delete($novel{$g});
                $similar{$g} = $t;
                last;
            } 
        }
    }
    $LOGGER->info("   Novel loci: ".scalar(keys %novel));
    $LOGGER->info("Collecting genes's info from gtf ...");
    #Group1  Cufflinks       transcript      234561  234772  1000    .       .       gene_id "CUFF.2"; transcript_id "CUFF.2.1"; FPKM "9.2736794454"; frac "1.000000"; conf_lo "0.884430";  conf_hi "3.390314"; cov "14.837959"; full_read_support "yes";
    open(GTF, "<", $gtfile) or $LOGGER->logdie($!);
    while(<GTF>) {
        next if ($_=~/^#/);
        chomp;
        my ($seq, $src, $type, $start, $end, $score, $strand, $phase, $extra)= split(/\t/, $_);
        if ($type eq 'transcript') {
            my ($id)=$extra=~/transcript_id "([^"]+)"/;
            my ($parent)=$extra=~/gene_id "([^"]+)"/;
            unless (exists $genetranscript{$id}) {
                if (exists $similar{$parent}) {
                    $genetranscript{$id} = $genetranscript{ $similar{$parent} };
                    $LOGGER->logdie("Not found gene for $similar{$parent}") unless ($genetranscript{$id});
                } else {
                    $genetranscript{$id}  = $parent;
                    if ($id =~/CUFF/) {
                        $genes{$parent} = 0;
                    } else {
                        $genes{$parent} = 1;
                    }
                }

            }
        }
    }
    close(GTF);
    
    my $zero = 0;
    $LOGGER->info("Calculating gene's length ...");
    foreach my $trans (keys %genetranscript) {
        # asm
        if ($genes{ $genetranscript{ $trans } } == 0) {
            $LOGGER->logdie("Not found length for $trans (asm)") unless (exists $translenasm{ $trans });
        
            if ($genelen{ $genetranscript{ $trans } }) {
                if ($genelen{ $genetranscript{ $trans } } < $translenasm{ $trans }) {
                    $genelen{ $genetranscript{ $trans } } = $translenasm{ $trans };
                }
            } else {
                $genelen{ $genetranscript{ $trans } } = $translenasm{ $trans };
            }
        } else {
            my $len;
            unless (exists $translenknown{ $trans }) {
                unless (exists $translenasm{ $trans }) {
                    $LOGGER->logdie("Not found length for $trans (known or asm)");
                }
                $len = $translenasm{ $trans };
            } else {
                $len = $translenknown{ $trans };
            }
        
            if ( $genelen{ $genetranscript{$trans} } ) {
                if ($genelen{ $genetranscript{$trans} } < $len ) {
                    $genelen{ $genetranscript{$trans} } = $len;
                }
            } else {
                $genelen{ $genetranscript{ $trans } } = $len;
            }
        }
    }
    open(GENES, ">", "$genhisto/genes.txt") or $LOGGER->logdie($!);
    foreach my $gene_id (keys %genes) {
        print GENES join("\t", $gene_id, $genes{ $gene_id }),"\n";
        if ($genes{ $gene_id } == 0) {
            $zero++;
        }
    }
    close(GENES);
    
    $LOGGER->logdie("Novel genes divergence: ".scalar(keys %novel)." != ".$zero) if ($zero != scalar(keys %novel));

    store \%genetranscript, "$genhisto/transcripts.dump";
    store \%genelen, "$genhisto/genelength.dump";

    $genes_md5 = file_md5_hex("$genhisto/genes.txt");
    unlink("$genhisto/session.RData") if (-e "$genhisto/session.RData");
} else {
        $LOGGER->logdie("Not found $genhisto/genelength.dump") unless (-e "$genhisto/genelength.dump");
        $LOGGER->info("Retrieve gene's length from dump ...");
        %genelen = %{ retrieve("$genhisto/genelength.dump") };
}

print CS basename("$genhisto/genes.txt"), "\t", $genes_md5, "\n";


my $reads_md5 = '';
if (-e "$genhisto/reads.txt") {
    $reads_md5 = file_md5_hex("$genhisto/reads.txt");
}

if ((!-e "$genhisto/reads.txt")||
    (! exists $checksum{basename("$genhisto/reads.txt")})||
    ( $checksum{basename("$genhisto/reads.txt")} ne $reads_md5 ) ) {
    
    if (scalar(keys %genetranscript) == 0) {
        $LOGGER->logdie("Not found $genhisto/transcripts.dump") unless (-e "$genhisto/transcripts.dump");
        $LOGGER->info("Retrieve gene transcript relationship from dump ...");
        %genetranscript = %{ retrieve("$genhisto/transcripts.dump") };
    }


    $LOGGER->info("Collecting reads's info ...");
    
    my %reads;
    
    $LOGGER->info(" ... from fastq ...");
    open(FASTQ, "<", $fqfile) or $LOGGER->logdie($!);
    while(<FASTQ>) {
        if ( ($. % 4) == 1 ) {
            if ($_=~/^@([^\/]+)/) {
                $reads{$1}=undef;
            }
        } 
    }
    close(FASTQ);
        
    $LOGGER->info(" ... from intersectBed output (known genes) ...");
    open(INTER, "<", "$genhisto/accepted_hits_x_known_genes.txt") or $LOGGER->logdie($!);
    
    while(<INTER>) {
        chomp;
        #Group4  4017326 4017376 OFQ-1:166:D15UPACXX:2:1101:10000:100498 50      -       Group4  4017263 4017489 GB49543-RA#exon3        .       +       LBDAv0.4        exon    .       ID=GB49543-RA#exon3;Parent=GB49543-RA        50
        my ($read, $attrs, $overlap, $transcript_id);

        (undef, undef, undef, $read, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, $attrs, $overlap) = split(/\t/, $_);
        if ($attrs=~/Parent=([^;]+)/) {
            $transcript_id=$1;
        } else {
            $LOGGER->logdie("Cannot found transcript_id attribute ($attrs)");
        }
        my $gene_id;
        if (exists $genetranscript{ $transcript_id } ) {
            $gene_id=$genetranscript{ $transcript_id };
        } else {
            $LOGGER->logdie("Cannot found gene_id for $transcript_id");
        }

        $reads{$read}->{$gene_id} = undef;
    }
    close(INTER);
    
    $LOGGER->info(" ... from intersectBed output (assembled genes) ...");
    open(INTER, "<", "$genhisto/accepted_hits_x_assembled_genes.txt") or $LOGGER->logdie($!);
   
    while(<INTER>) {
        chomp;
        #Group11 8683594 8683644 OFQ-1:166:D15UPACXX:2:1101:10001:37284  50      +       Group11 8681261 8687153 CUFF.3236       593     -       Cufflinks       exon    .       gene_id "CUFF.3236"; transcript_id "CUFF.3236.5"; exon_number "1"; FPKM "24.3786885491"; frac "0.406885"; conf_lo "20.811306"; conf_hi "22.871057"; cov "38.910667"; 50
        my ($read, $attrs, $overlap, $transcript_id);

        (undef, undef, undef, $read, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, $attrs, $overlap) = split(/\t/, $_);
        if (($attrs) && ($attrs ne '.')) {
            if ($attrs=~/transcript_id "([^"]+)"/) {
                $transcript_id=$1;
            } else {
                $LOGGER->logdie("Cannot found transcript_id attribute ($attrs): $_");
            }
            my $gene_id;
            if (exists $genetranscript{ $transcript_id } ) {
                $gene_id=$genetranscript{ $transcript_id };
            } else {
                $LOGGER->logdie("Cannot found gene_id for $transcript_id");
            }
            $reads{$read}->{$gene_id} = undef;
        } else {
            $reads{$read}->{''} = undef;
        }
    }
    close(INTER);
    
    my $nreads=0;
    $LOGGER->info(" Writting results (reads.txt) ...");
    open(OUT, ">", "$genhisto/reads.txt") or $LOGGER->logdie($!);
    foreach my $read (keys %reads) {
        $nreads++;
        if ($reads{ $read }) {
            my $n = scalar(keys %{ $reads{ $read } });
            foreach my $gene_id (keys %{ $reads{ $read } }) {
                if ($gene_id ne '') { 
                    print OUT join("\t", $nreads, $read, $gene_id, (1/$n), ((1/$n)/($genelen{ $gene_id }/1000))),"\n";
                } else {
                    print OUT join("\t", $nreads, $read, $gene_id, (1/$n), 0),"\n";
                }                    
            }
        } else {
                print OUT join("\t", $nreads, $read, '', 0, 0),"\n";
        }
    }
    close(OUT);
    &R::eval('n.reads <<- '.$nreads);
    
    $reads_md5 = file_md5_hex("$genhisto/reads.txt");
    unlink("$genhisto/session.RData") if (-e "$genhisto/session.RData");
}

print CS basename("$genhisto/reads.txt"), "\t", $reads_md5, "\n";

close(CS);

&R::eval('library(package="ggplot2")');
&R::eval('library(package="scales")');
&R::eval('library(package="data.table")');
&R::eval('library(package="plyr")');
&R::eval('
# based on findgr() function (http://www.princeton.edu/genomics/botstein/protocols/Growth-Rate-Using-R.pdf)
find.line <<- function(x, y, int, r2.cutoff) {
    n = length(x)
    mat = NULL
    fit = lm(x~y)
    m = abs(coefficients(fit)[[2]])
    if (m < 0.00001) {
        return( c(start=y[1], end=length(y), slope=0, r2=0, intercept=0  ) )
    } else {

        x[which(x <= 0)] = 0.001 #transform values < 0
        for (i in 1:(n-int-1)) {
          fit = lm(x[i:(i+int-1)]~y[i:(i+int-1)])               #linear regression on log transformed data.
          m = coefficients(fit)[[2]]
          b = coefficients(fit)[[1]]
          r2 = summary(fit)$r.squared
          mat = rbind(mat, c(i, b, m, r2))
        }
        mat = mat[which(mat[,4] > r2.cutoff),]  #only include slopes greater than the R2 cutoff.
        min = mat[which.min(mat[,3]),]
        if ( length(min) > 0 ) {
            fit.line = sapply(y, function(x) min[3]*x+min[2])
            resid = fit.line-x
            residper = round(abs(resid/fit.line), 2)
            resid.mat = rbind(y, fit.line, x, resid, residper)
            resid.mat = resid.mat[,which(resid.mat[5,] <= 0.01)]
            
            par(las=1, mar=c(5, 4, 4, 4) + 0.1)
            plot(y,x, type="n", pch=20, xlab="x-values", ylab="logn(y-values)", main="TEST")
            usr.old = par("usr")
            mtext(paste("lag =",round( resid.mat[1,1] ,2)),side=3, line=-3, at=0, cex=0.8, adj=0)
            abline(v=resid.mat[1,1], col="cadet blue", lty=2)
            abline(v=resid.mat[1,ncol(resid.mat)], col="cadet blue", lty=2)
            par(usr=usr.old)
            points(y,x,pch=20)
            abline(a=min[2], b=min[3], col="red", lty=2, lwd=2)
            points(y[min[1]:(min[1]+int-1)], x[min[1]:(min[1]+int-1)], col="red")
            mtext(paste("m =",round(min[3],3)), side=3, line=-1, at=0, cex=0.8, adj=0)
            mtext(paste("r2 =",round(min[4],4)), side=3, line=-2, at=0, cex=0.8, adj=0)
            dx = diff(x)/diff(y)
            #par(usr=c(par("usr")[1:2],min(dx)*1.05, min(dx)*1.05))
            points(y[1:(length(y)-1)],dx, pch=18, type="o", col="dark grey", lty=1)
            axis(4, col.axis="dark grey", col.ticks="dark grey")
            mtext("delta(x)/delta(y)", side=4, line=3, col="dark grey", las=3)
            #print(c(start=as.numeric(resid.mat[1,1]), end=as.numeric(resid.mat[1,ncol(resid.mat)])))
               
            return( c(start=as.numeric(resid.mat[1,1]), end=as.numeric(resid.mat[1,ncol(resid.mat)]), slope=round(min[3],2), r2=round(min[4],2), intercept=min[2] ) )
        } else {
            return( c(start=y[1], end=length(y), slope=0, r2=0, intercept=0 ) )
        }
    }
}
');

unless (-e "$genhisto/session.RData") {

    $LOGGER->info("Loading genes ...");
    &R::eval('genes.df <<- read.delim(file="'.File::Spec->rel2abs( "$genhisto/genes.txt" ).'", header=FALSE, col.names=c("gene", "known"), stringsAsFactors = FALSE)
              ngenes <<- dim(genes.df)[1]
              known.genes <<- as.character(genes.df[ which(genes.df$known == 1)  ,"gene"])
              unknown.genes <<- as.character(genes.df[ which(genes.df$known == 0)  ,"gene"])
             ');

    # Número total de reads = leituras não alinhadas (count==0) + soma das contagens das leituras alinhadas
    # length(unique(exp.reads[which(exp.reads$count==0),"id"]))+sum(exp.reads$count)

    $LOGGER->info("Reading reads.txt file ...");
    &R::eval('exp.reads <<- read.delim(file="'.File::Spec->rel2abs( "$genhisto/reads.txt" ).'", header=FALSE, stringsAsFactors = FALSE)
        colnames(exp.reads) <<- c("id", "name", "gene", "count", "rpk")
        ');
    
    $LOGGER->info("Sampling ...");
    &R::eval('rvar <<- data.frame( "X1"=sample(1:n.reads, n.reads) )');

    $LOGGER->info("Saving session ...");
    &R::eval('save(list = ls(all=TRUE,envir=globalenv()), file ="'.File::Spec->rel2abs( "$genhisto/session.RData" ).'" )');
} else {
    $LOGGER->info("Loading saved session ...");
    &R::eval('load(file="'.File::Spec->rel2abs( "$genhisto/session.RData" ).'", envir = globalenv())');
}

&R::eval('
n.reps <<- '.$reps.'
n.div <<- '.$div.'
n.min <<- '.$minreads.'
n.minrpkm <<- '.$minrpkm.'

m <<- data.frame( n = c(1:(n.reads/n.div)))');

if ($reps > 1) {
    for (my $r=2; $r<=$reps; $r++) {
        &R::eval('rvar$X'.$r.' <<- sample(1:n.reads, n.reads)');
    } 
}
 
$LOGGER->info("Processing data ...");

&R::eval('
mapped.pe <<- NULL

for (r in 1:n.reps) {
    # selecionados
    sel <<- list()
    sel[[as.character(0)]] <<- c()

    # sobras
    left <<- list()
    left[[as.character(0)]] <<- c()
    
    last.mapped.reads <<- 0
    last.mapped.reads.unique <<- 0

    for  (i in 1:dim(m)[1]) {

        # conjunto de dados corrente - iteração i
        cur <<- subset(exp.reads, id %in% rvar[[paste("X",r,sep="")]][(((i-1)*n.div)+1):((i)*n.div)] )[, c("gene", "count", "rpk")]

        # mapped reads
        mapped.total <<- last.mapped.reads + sum(cur$count)

        # somente irá filtrar com "cur" pois contém os mapeados fora de regiões gênicas
        mapped.unique <<- last.mapped.reads.unique + sum((subset(cur, count >= 1))$count)
        
        mapped.pe <<- rbind(mapped.pe, c( i, "total", mapped.total ))
        mapped.pe <<- rbind(mapped.pe, c( i, "unique", mapped.unique ))

        last.mapped.reads <<- mapped.total
        last.mapped.reads.unique <<- mapped.unique
            
        # conjunto de dados corrente - iteração i - exceto os contabilizados anteriormente
        cur.sel <<- subset(cur, ! gene %in% c(sel[[as.character(i-1)]]$gene,"") )
        cur.sel$rpkm <<- 0
            
        # os que poderão ser contabilizados (cur.sel) mais os que não foram contabilizados anteriormente, i.e., as sobras.
        cur.union <<- rbind( cur.sel, left[[as.character(i-1)]] )
        if ( dim(cur.union)[1] > 0) {
            
            cur.union.aggreg <<- aggregate(. ~ gene, data=cur.union, sum)

            # somente irá filtrar com "cur" pois contém os mapeados fora de regiões gênicas
            cur.union.aggreg.unique <<- aggregate(. ~ gene, data=subset(cur, count >= 1), sum)
            

            # rpkm
            cur.union.aggreg$rpkm <<- (cur.union.aggreg$rpk/( mapped.total/1000000 ))

            # selecionados nesta iteração
            cur.union.aggreg.sel <<- subset(cur.union.aggreg, count >= n.min & rpk >= n.minrpkm)   

            # os selecionados nesta iteração junto com os selecionados na anterior
            sel[[as.character(i)]] <<- rbind(sel[[as.character(i-1)]] , cur.union.aggreg.sel)

            # sobras desta iteração
            left[[as.character(i)]] <<- subset(cur.union.aggreg, ! gene %in% c( cur.union.aggreg.sel$gene, "")  )
        } else {
            sel[[as.character(i)]] <<- cur.union
            left[[as.character(i)]] <<- cur.union
        }

        m[i, paste("genes", r,sep="")] <<- dim( sel[[as.character(i)]] )[1]
        m[i, paste("known", r,sep="")] <<- length( grep(TRUE, sel[[as.character(i)]]$gene %in% known.genes) ) 
        m[i, paste("unknown", r,sep="")] <<- m[i, paste("genes", r,sep="")]-m[i, paste("known", r,sep="")] 
    }
}

colnames(mapped.pe) <<- c("n", "class", "count")
mapped.pe <<- as.data.frame(mapped.pe)
mapped.pe$n <<- as.numeric(as.character(mapped.pe$n))
mapped.pe$class <<- as.factor(as.character(mapped.pe$class))
mapped.pe$count <<- as.numeric(as.character(mapped.pe$count))

mapped.stats <<- ddply(mapped.pe, c("class","n"), function(df) return(c(hwy.avg=mean(df$count), hwy.sd=sd(df$count))))
mapped.stats$hwy.se <- (mapped.stats$hwy.sd/sqrt(n.reps))

mapped.stats$avg.perc <- apply( mapped.stats, 1, function(x,d) { return(((as.numeric(x[["hwy.avg"]]))/(as.numeric(x[["n"]])*d))*100) }, n.div)
mapped.stats$avg.up.perc <- apply( mapped.stats, 1, function(x,d,v) { return(((as.numeric(x[["hwy.avg"]])+(as.numeric(x[["hwy.se"]])*v) )/(as.numeric(x[["n"]])*d))*100) }, n.div, 1.96)
mapped.stats$avg.down.perc <- apply( mapped.stats, 1, function(x,d,v) { return(((as.numeric(x[["hwy.avg"]])-(as.numeric(x[["hwy.se"]])*v) )/(as.numeric(x[["n"]])*d))*100) }, n.div, 1.96)


ggplot(data =mapped.stats, aes(x=n,y=avg.perc, fill=class)) + geom_bar(position="dodge", stat="identity", width=0.9)+theme_bw()+theme(axis.text.x=element_text(angle=60, hjust=1) )+scale_fill_discrete(name="Type", breaks=c("unique", "total"), labels=c("Unique", "Total"))+ylab("Percentage of aligned reads (%)")+xlab("Number of sequenced reads")+ggtitle("Alignments of Reads")+scale_y_continuous(breaks = round(seq(min(mapped.stats$avg.down.perc), 100, by = 5),2))+scale_x_continuous(breaks = round(seq(0, max(mapped.stats$n)),1))+geom_errorbar(aes(ymin=avg.down.perc, ymax=avg.up.perc, width=0.25), position = position_dodge(width = 0.90))
ggsave(file="'."$genhisto/genhistreads.png".'")


');
    
#&R::eval('save(list = ls(all=TRUE,envir=globalenv()), file ="'.File::Spec->rel2abs( "$genhisto/DEBUG.RData" ).'" )');

&R::eval('


sampled.df <<- data.frame(n=m[,"n"])
sampled.df$mean <<- apply(as.data.frame(m[, grep("^genes", colnames(m)) ]), 1, mean)
sampled.df$se <<- apply(as.data.frame(m[, grep("^genes", colnames(m)) ]), 1, sd)/sqrt(n.reps)
sampled.df[["se"]][is.na(sampled.df[["se"]])] <<- 0
sampled.df$type <<- "sampled"
sampled.df$status <<- "all"

lo.model <<- lm(mean ~ log(n), data=sampled.df)

extra.df <<- data.frame(n=seq(1, (max(sampled.df$n)+as.integer((max(sampled.df$n)/3)))   ))

extra.pred <<- predict(lo.model, extra.df, se=TRUE)
extra.df$mean <<- extra.pred$fit
extra.df$se <<- extra.pred$se.fit
extra.df$type <<- "predicted"
extra.df$status <<- "all"

known.df <<- data.frame(n=m[,"n"])
known.df$mean <<- apply(as.data.frame(m[, grep("^known", colnames(m)) ]), 1, mean)
known.df$se <<- apply(as.data.frame(m[, grep("^known", colnames(m)) ]), 1, sd)/sqrt(n.reps)
known.df[["se"]][is.na(known.df[["se"]])] <<- 0
known.df$type <<- "sampled"
known.df$status <<- "known"

lo.known.model <<- lm(mean ~ log(n), data=known.df)

known.extra.df <<- data.frame(n=seq(1, (max(known.df$n)+as.integer((max(known.df$n)/3)))   ))

known.extra.pred <<- predict(lo.known.model, known.extra.df, se=TRUE)
known.extra.df$mean <<- known.extra.pred$fit
known.extra.df$se <<- known.extra.pred$se.fit
known.extra.df$type <<- "predicted"
known.extra.df$status <<- "known"

unknown.df <<- data.frame(n=m[,"n"])
unknown.df$mean <<- apply(as.data.frame(m[, grep("^unknown", colnames(m)) ]), 1, mean)
unknown.df$se <<- apply(as.data.frame(m[, grep("^unknown", colnames(m)) ]), 1, sd)/sqrt(n.reps)
unknown.df[["se"]][is.na(unknown.df[["se"]])] <<- 0
unknown.df$type <<- "sampled"
unknown.df$status <<- "novel"

lo.unknown.model <<- lm(mean ~ log(n), data=unknown.df)

unknown.extra.df <<- data.frame(n=seq(1, (max(unknown.df$n)+as.integer((max(unknown.df$n)/3)))   ))

unknown.extra.pred <<- predict(lo.unknown.model, unknown.extra.df, se=TRUE)
unknown.extra.df$mean <<- unknown.extra.pred$fit
unknown.extra.df$se <<- unknown.extra.pred$se.fit
unknown.extra.df$type <<- "predicted"
unknown.extra.df$status <<- "novel"

all.df <<- rbind( sampled.df, extra.df, known.df, known.extra.df, unknown.df, unknown.extra.df)

#save(list = ls(all=TRUE,envir=globalenv()), file="/tmp/DEBUG.RData");

png(file="'."$genhisto/genhistaux.png".'", bg="transparent", width=800, height=600)
line.limits <<- find.line( as.numeric(c(0,(subset(all.df, type=="sampled" & status=="all"))$mean)), 
                         as.numeric(c(0,(subset(all.df, type=="sampled" & status=="all"))$n)), 5, 0.98)
dev.off()


#  1.96 is the .975 quantile of the normal distribution = upper and lower 95% 

mainplot <<- ggplot(all.df, aes(x=n, y=mean, linetype=type, colour=status))+geom_line(stat="identity", position="identity", lwd=0.25)+geom_ribbon(aes(ymax=mean+(se*1.96), ymin=mean-(se*1.96), fill=status), colour=NA, alpha=0.2)+xlab(paste("Number of reads ( x ", n.div, " )", sep="") )+ylab(paste("Number of genes (more than ", n.min, " read(s) and ", n.minrpkm, " rpkm)",sep=""))+ggtitle("Sequencing Depths")+geom_hline(yintercept=dim(genes.df)[1], colour="red", alpha=0.2)+geom_hline(yintercept=length(known.genes), colour="green", alpha=0.2)+geom_hline(yintercept=length(unknown.genes), colour="blue", alpha=0.2)+scale_x_continuous(breaks = round(seq(0, max(all.df$n)),1))+theme_bw()+theme(axis.text.x=element_text(angle=60, hjust=1) )+scale_y_continuous()+scale_linetype_manual(values = c("sampled" = "solid", "predicted" = "dashed"))+scale_colour_manual(values = c("all" = "red", "known" = "green", "novel"="blue" ))

if (line.limits[["intercept"]] > 0) {
    mainplot+
    geom_vline(xintercept=line.limits[["start"]], colour="red", linetype="twodash", alpha=0.2)+
    geom_vline(xintercept=line.limits[["end"]], colour="red", linetype="twodash", alpha=0.2)+
    annotate("text", hjust=0, label = paste("r2=",line.limits[["r2"]], sep=""), x = 0, y = 21000, size = 4, colour = "gray10")+
    annotate("text", hjust=0, label = paste("m=",line.limits[["slope"]], sep=""), x = 0, y = 20000, size = 4, colour = "gray10")+
    annotate("text", hjust=0, label = paste("b=",line.limits[["intercept"]], sep=""), x = 0, y = 19000, size = 4, colour = "gray10")+
    geom_abline(intercept=line.limits[["intercept"]], slope=line.limits[["slope"]], colour="gray10",linetype="dotted")
} else {
    mainplot
}

ggsave(file="'."$genhisto/genhist.png".'")

predicted.all.diff <<- data.frame(n=as.numeric(c(0,(subset(all.df, type=="predicted" & status=="all"))$n) ))

predicted.all.diff$diff <<- c(1, diff(as.numeric(c(0,(subset(all.df, type=="predicted" & status=="all"))$mean)))/diff(as.numeric(c(0,(subset(all.df, type=="predicted" & status=="all"))$n))))

ggplot(predicted.all.diff,aes(y=diff,x=n))+ geom_bar(stat="identity",fill = I("grey50"), width=0.8)+ geom_line()+ scale_x_continuous(breaks = round(seq(0, max(predicted.all.diff$n)),1))+ scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+ theme_bw()+theme(axis.text.x=element_text(angle=60, hjust=1), panel.background = element_rect(fill = NA, colour="gray"))+ xlab(paste("Number of reads ( x ", n.div, " )", sep="") )+ ylab(paste("Number of increased genes (more than ", n.min, " read(s) and ", n.minrpkm, " rpkm)",sep=""))+ ggtitle("Sequencing Depths")

ggsave(file="'."$genhisto/genhistdiff.png".'")

');

close(CS);

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help       Help
        -l      --level      Log level [Default: FATAL]
        -f      --fqfile     fastq file
        -a      --bamfilea   bam file (alignment)
        -t      --gtfile     gtf (cufflinks output: transcripts.gtf)
        -r      --gffile     gff (gene model reference for cuffcompare)
        -g      --genhisto   genhist output
        -d      --division   Divided by n [Default: 100000]
        -m      --minreads   Minimum reads to be accounted [Default: 1]
        -k      --minrpkm    Minimum rpkm to be accounted [Default: 0]
        -n      --replicates Number of replicates (sampling)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

