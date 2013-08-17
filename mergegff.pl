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

use Bio::Tools::GFF;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $gff1, $gff2, $genome, $auxfile, $gffout, $source_tag);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "g1|gff1=s"=> \$gff1,
            "g2|gff2=s"=> \$gff2,
            "g|genome=s"=>\$genome,
            "a|auxfile=s"=>\$auxfile,
            "o|gffout=s"=>\$gffout,
            "s|source=s"=>\$source_tag
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

$LOGGER->logdie("Missing gff1") unless ($gff1);
$LOGGER->logdie("Wrong gff1 ($gff1)") unless (-e $gff1);

$LOGGER->logdie("Missing gff2") unless ($gff2);
$LOGGER->logdie("Wrong gff2 ($gff2)") unless (-e $gff2);

$LOGGER->logdie("Missing genome") unless ($genome);
$LOGGER->logdie("Wrong genome ($genome)") unless (-e $genome);

$LOGGER->logdie("Missing auxfile") unless ($auxfile);
open(AUX, ">", $auxfile) or $LOGGER->logdie($!);

$LOGGER->logdie("Missing gffout") unless ($gffout);
open(GFFOUT, ">", $gffout) or $LOGGER->logdie($!);
print GFFOUT '##gff-version 3',"\n";

$LOGGER->logdie("Missing source") unless ($source_tag);


my %chromosome;

my %ptag =  (   'gene'  =>1, 
                'mRNA'  =>1, 
                'exon'  =>1, 
                'CDS'   =>1 
            );

my %feature;
my %gene;
my %gff1counter;
my %group;
my $gffio1 = Bio::Tools::GFF->new(-file =>$gff1, -gff_version => 3);
while (my $feat = $gffio1->next_feature()) {
    my ($id,$parent);
    
    next unless (exists $ptag{ $feat->primary_tag() });
    
    foreach my $tag ( $feat->get_all_tags() ) {
        unless ($tag =~ /ID|Parent/) {
            $feat->remove_tag($tag);
        }
    }

    if ($feat->has_tag('ID')) {
        ($id) = $feat->get_tag_values('ID');
    }
    else {
        if ($feat->has_tag('Parent')) {
            ($id) = $feat->get_tag_values('Parent');
            if (($feat->primary_tag() eq 'CDS')&&($id =~/^GB\d+-R/)) {
                $id=~s/-R([ABC])/-P$1/;
            } 
        }
        else {
            $id = '';
        }    
        $gff1counter{ $feat->primary_tag() }->{ $id } = 0 unless ( exists $gff1counter{ $feat->primary_tag() }->{ $id } );
        $gff1counter{ $feat->primary_tag() }->{ $id }++;
        
        $id.='#'.$feat->primary_tag().$gff1counter{ $feat->primary_tag() }->{$id};

        $feat->add_tag_value('ID', $id);
    }
    

    unless (exists $feature{$id}) {
        $feature{$id} = {   'primary_tag'   => $feat->primary_tag(),
                            'child'         => undef };
        push(@{$chromosome{ $feat->seq_id() }}, $id);
    }                            
    

    if ($feat->has_tag('Parent')) {
        ($parent) = $feat->get_tag_values('Parent');

        push(@{$feature{$parent}->{'child'}->{ $feat->primary_tag() }}, $feat);
    }
    else {
        $gene{ $feat->seq_id() }->{$id} = $feat;
        $group{ $feat->seq_id() }->{'A'}->{$id} = undef;
    }

}

my %gff2counter;

my $gffio2 = Bio::Tools::GFF->new(-file =>$gff2, -gff_version => 3);
while (my $feat = $gffio2->next_feature()) {
    my ($id,$parent);
    
    next unless (exists $ptag{ $feat->primary_tag() });
    
    foreach my $tag ( $feat->get_all_tags() ) {
        unless ($tag =~ /ID|Parent/) {
            $feat->remove_tag($tag);
        }
    }

    if ($feat->has_tag('ID')) {
        ($id) = $feat->get_tag_values('ID');
    }
    else {
        if ($feat->has_tag('Parent')) {
            ($id) = $feat->get_tag_values('Parent');
            if (($feat->primary_tag() eq 'CDS')&&($id =~/^GB\d+-R/)) {
                $id=~s/-R([ABC])/-P$1/;
            } 
        }
        else {
            $id = '';
        }    
        $gff2counter{ $feat->primary_tag() }->{ $id } = 0 unless ( exists $gff2counter{ $feat->primary_tag() }->{ $id } );
        $gff2counter{ $feat->primary_tag() }->{ $id }++;
        
        $id.='#'.$feat->primary_tag().$gff2counter{ $feat->primary_tag() }->{$id};

        $feat->add_tag_value('ID', $id);

    }
    
    unless (exists $feature{$id}) {
        $feature{$id} = {   'primary_tag'   => $feat->primary_tag(),
                            'child'         => undef };
        push(@{$chromosome{ $feat->seq_id() }}, $id);
    }                            
    

    if ($feat->has_tag('Parent')) {
        ($parent) = $feat->get_tag_values('Parent');

        push(@{$feature{$parent}->{'child'}->{ $feat->primary_tag() }}, $feat);
    }
    else {
        $gene{ $feat->seq_id() }->{$id} = $feat;
        $group{ $feat->seq_id() }->{'B'}->{$id} = undef;
    }

    
}

my $tmpdir = tempdir( CLEANUP => 0 );


print ">$tmpdir\n";

foreach my $chr (keys %chromosome) {
        
    my $chrdir = $tmpdir.'/'.$chr;

    mkdir($chrdir);

    my ($filename1, $filename2);

    {
        my $fh1;

        ($fh1, $filename1) = tempfile( DIR => $chrdir );
    #    print "\t$filename1\n";

        foreach my $g (keys %{ $group{$chr}->{'A'} }) {
            my $f = $gene{ $chr }->{$g};
            print { $fh1 } join("\t", $chr, $f->start(), $f->end(), $g, 0, (($f->strand() > 0) ? '+' : '-')),"\n";
        }
        
        my $fh2;

        ($fh2, $filename2) = tempfile( DIR => $chrdir );
    #    print "\t$filename2\n";

        foreach my $g (keys %{ $group{ $chr }->{'B'} }) {
            my $f = $gene{ $chr }->{$g};
            print { $fh2 } join("\t", $chr, $f->start(), $f->end(), $g, 0, (($f->strand() > 0) ? '+' : '-')),"\n";
        }
    }

    my $i = 1;
    my %mitem;
    

    if (( ! -z $filename1)&&( ! -z $filename2)) {
#        print "RUNNING BEDTools -a $filename1 -b $filename2 \n";
        `intersectBed -bed -wo -s -a $filename1 -b $filename2 > $chrdir/MERGE.txt`;
        my %merge;

        open(MERGE, "<", "$chrdir/MERGE.txt") or $LOGGER->logdie($!);
        while(<MERGE>) {
    #        Group2  4190154 4191124 GB50778 0       +       Group2  4188631 4293258 igl     0       +       970
            my ($chra, $starta, $enda, $namea, $scorea, $stranda, $chrb, $startb, $endb, $nameb, $scoreb, $strandb, $mcount) = split(/\t/, $_);

                $merge{$namea}->{$nameb} = undef;
                $merge{$nameb}->{$namea} = undef;

                $mitem{$namea}->{'data'} = [$chra, $starta, $enda, $namea, $scorea, $stranda];
                $mitem{$namea}->{'type'} = 'A';
                $mitem{$nameb}->{'data'} = [$chrb, $startb, $endb, $nameb, $scoreb, $strandb];
                $mitem{$nameb}->{'type'} = 'B';
        }
        close(MERGE);
        

#        foreach my $item (keys %merge) {
#            foreach my $content (keys %{$merge{$item}}) {
#                print $item,"\t",$content,"\n";
#            }
#        }

        my %RNA;
        foreach my $item (keys %merge) {
            next unless ( exists $merge{$item} );

            my %hitem;
            &recitem(\%merge, \%hitem, $item);
            my $fA = "$chrdir/ITEM_$i".'_A'.".txt";
            my $fB = "$chrdir/ITEM_$i".'_B'.".txt";
            open(ITEMA, ">", $fA) or $LOGGER->logdie($!);
            open(ITEMB, ">", $fB) or $LOGGER->logdie($!);
            foreach my $hi (keys %hitem) {
                if ($mitem{$hi}->{'type'} eq 'A') {
                    &printall(\*ITEMA, $gene{$chr}, \%feature, $hi, $gffio1, \%RNA);
                }
                elsif ($mitem{$hi}->{'type'} eq 'B') {
                    &printall(\*ITEMB, $gene{$chr}, \%feature, $hi, $gffio2, \%RNA);
                }
                else {
                    $LOGGER->logdie("Problem found");
                }
            }
            close(ITEMA);
            close(ITEMB);

            my $ffA = $fA;
            $ffA=~s/\.txt$/\.fa/;
            `cat $fA | gffread -g $genome -w $ffA -y $ffA.aa`;
            my @seqA;
            # RNAs
            {
                my $in=Bio::SeqIO->new(-file=>$ffA); 
                while(my $seq=$in->next_seq()) {
                    my $seqfilename = $seq->display_id();
                    $seqfilename=~s/\|/_/g;
                    my $out=Bio::SeqIO->new(-file=>">$chrdir/".$seqfilename.".fa");
                    $out->write_seq($seq);
                    push(@seqA, $seq->display_id());
                }
            }
            # proteínas
            {
                my $in=Bio::SeqIO->new(-file=>"$ffA.aa"); 
                while(my $seq=$in->next_seq()) {
                    my $seqfilename = $seq->display_id();
                    $seqfilename=~s/\|/_/g;
                    my $out=Bio::SeqIO->new(-file=>">$chrdir/".$seqfilename.".aa.fa");
                    $out->write_seq($seq);
                }
            }

            my $ffB = $fB;
            $ffB=~s/\.txt$/\.fa/;
            `cat $fB | gffread -g $genome -w $ffB -y $ffB.aa`;
            my @seqB;
            # RNAs
            {
                my $in=Bio::SeqIO->new(-file=>$ffB); 
                while(my $seq=$in->next_seq()) {
                    my $seqfilename = $seq->display_id();
                    $seqfilename=~s/\|/_/g;
                    my $out=Bio::SeqIO->new(-file=>">$chrdir/".$seqfilename.".fa");
                    $out->write_seq($seq);   
                    push(@seqB, $seq->display_id());
                }
            }
            # proteínas
            {
                my $in=Bio::SeqIO->new(-file=>"$ffB.aa"); 
                while(my $seq=$in->next_seq()) {
                    my $seqfilename = $seq->display_id();
                    $seqfilename=~s/\|/_/g;
                    my $out=Bio::SeqIO->new(-file=>">$chrdir/".$seqfilename.".aa.fa");
                    $out->write_seq($seq);   
                }
            }
            
            my %result;
            my %RNAequal;
            my %CDSequal;
            foreach my $sA (@seqA) {
                foreach my $sB (@seqB) {
                    my $seqfilenameA = $sA;
                    $seqfilenameA=~s/\|/_/g;
                    my $seqfilenameB = $sB;
                    $seqfilenameB=~s/\|/_/g;

                    `bl2seq -i $chrdir/$seqfilenameA.fa -j $chrdir/$seqfilenameB.fa -p blastn -F f -D 0 -o $chrdir/$seqfilenameA-$seqfilenameB.bl2`;
                        my $searchio = new Bio::SearchIO(-format => 'blast',
                                                         -file   => "$chrdir/$seqfilenameA-$seqfilenameB.bl2");
                        while( my $bl2res = $searchio->next_result() ) {
                                my $query_name=$bl2res->query_name();

                                if ($bl2res->hits()) {
                                        while( my $hit = $bl2res->next_hit() ) {
                                            my @hsp;
                                            while( my $h = $hit->next_hsp() ) {
                                                push(@hsp, $h);
                                            }
                                                
                                            my ($geneA) = $RNA{$sA}->{'gene'};
                                            my ($geneB) = $RNA{$sB}->{'gene'};

                                            
                                            if (($hsp[0]->score() > 100)&&($hsp[0]->evalue() < 1e-10)) {
                                                print AUX join("\t", $chr, $source_tag.'_'.$chr.'_'.sprintf("%05d",$i), $geneA, $sA, $geneB, $sB, $hsp[0]->percent_identity(), $hsp[0]->score(), $hsp[0]->evalue() ),"\n";
                                                $result{$geneA}->{$geneB} = 0 unless (exists $result{$geneA}->{$geneB});
                                                $result{$geneB}->{$geneA} = 0 unless (exists $result{$geneB}->{$geneA});
                                                if (($hit->query_length() == $hit->length()) &&
                                                    ($hsp[0]->percent_identity() == 100) &&
                                                    ($hsp[0]->length('hit') == $hit->query_length())
                                                ) { # exactly equal sequences

                                                    $RNAequal{$sA}->{$sB} = undef;
                                                    $RNAequal{$sB}->{$sA} = undef;
                                                    if ((-e "$chrdir/$seqfilenameA.aa.fa")&&(-e "$chrdir/$seqfilenameB.aa.fa")) {
                                                        `bl2seq -i $chrdir/$seqfilenameA.aa.fa -j $chrdir/$seqfilenameB.aa.fa -p blastp -F f -D 0 -o $chrdir/$seqfilenameA-$seqfilenameB.aa.bl2`;
                                                        my $aasearchio = new Bio::SearchIO(-format => 'blast',
                                                                                           -file   => "$chrdir/$seqfilenameA-$seqfilenameB.aa.bl2");

                                                        while( my $aabl2res = $aasearchio->next_result() ) {
                                                            my $query_name=$aabl2res->query_name();
                                                            if ($aabl2res->hits()) {
                                                                while( my $aahit = $aabl2res->next_hit() ) {
                                                                    my @aahsp;
                                                                    while( my $h = $aahit->next_hsp() ) {
                                                                        push(@aahsp, $h);
                                                                    }
                                                
                                                                    if (($aahit->query_length() == $aahit->length()) &&
                                                                        ($aahsp[0]->percent_identity() == 100) &&
                                                                        ($hsp[0]->length('hit') == $hit->query_length())
                                                                    ) {
                                                                        $CDSequal{$sA}->{$sB} = undef;
                                                                        $CDSequal{$sB}->{$sA} = undef;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }

                                                }                                                    
                                            }
                                            last;
                                        }
                                }
                                last;
                        }
                }
            }
            
            my %groupgene;
            my %genegroup;
            foreach my $k (keys %result) {
                if (exists $result{$k}) {
                    my %groupname;
                    &recgroupname(\%result, $k, \%groupname);
                    my $genegroupname = join('/', keys %groupname);
                    push(@{$groupgene{ $genegroupname }}, keys %groupname);
                    foreach my $gn (keys %groupname) {
                        $genegroup{$gn} = $genegroupname;
                    }
                }
            }
            
#            foreach my $gn (keys %groupgene) {
#                print $gn,"\t",join(";", @{ $groupgene{$gn} }),"\n";
#            }
            
            foreach my $gn (keys %groupgene) {
                my ($gstart, $gend) = (9999999999999999999, 0);
                my $strand = undef;
                foreach my $g (@{ $groupgene{ $gn } }) {
                    my $f = $gene{ $chr }->{$g};
                    unless ($strand) {
                        $strand = $f->strand();
                    }
                    else {
                        $LOGGER->logdie("The strand cannot be different here") if ($strand != $f->strand());
                    }
                    if ($f->start() < $gstart) {
                        $gstart = $f->start();
                    }
                    if ($f->end() > $gend) {
                        $gend = $f->end();
                    }

                }
                
                my $feat = new Bio::SeqFeature::Generic ( -start => $gstart, -end => $gend,
                                -strand => $strand, -primary => 'gene',
                                -source_tag   => $source_tag,
                                -seq_id => $chr,
                                -score  => '.',
                                -tag    => {
                                    ID => $gn,
                                    Note => 'LOCUS:'.$source_tag.'_'.$chr.'_'.sprintf("%05d",$i) } );

                print GFFOUT $gffio1->gff_string($feat),"\n";
                
                foreach my $g (@{ $groupgene{ $gn } }) {
                    &printall(\*GFFOUT, $gene{$chr}, \%feature, $g, $gffio1, {}, { mRNA=>{'Parent'=>$gn} },1,\%RNAequal,\%CDSequal);
                }
                
            }

            foreach my $hi (keys %hitem) {
                unless (exists $genegroup{$hi}) {
                    &printall(\*GFFOUT, $gene{$chr}, \%feature, $hi, $gffio1, {}, {gene=>{'Note'=>'LOCUS:'.$source_tag.'_'.$chr.'_'.sprintf("%05d",$i)} });
                }
            }

            $i++;
        }
    }

    my %RNA;
        open(F1, "<", $filename1) or $LOGGER->logdie($!);
        while(<F1>) {
            chomp;
            my ($chr, $start, $end, $name, $score, $strand) = split(/\t/, $_);
            next if (exists $mitem{$name});

#            open(ITEM, ">", "$chrdir/ITEM_$i.txt") or $LOGGER->logdie($!);
#            &printall(\*ITEM, $gene{$chr}, \%feature, $name, $gffio1, \%RNA, {});
#            close(ITEM);
                    
            &printall(\*GFFOUT, $gene{$chr}, \%feature, $name, $gffio1, {}, {gene=>{'Note'=>'LOCUS:'.$source_tag.'_'.$chr.'_'.sprintf("%05d",$i)} });
            $i++;
        }
        close(F1);

        open(F2, "<", $filename2) or $LOGGER->logdie($!);
        while(<F2>) {
            chomp;
            my ($chr, $start, $end, $name, $score, $strand) = split(/\t/, $_);
            next if (exists $mitem{$name});
#            open(ITEM, ">", "$chrdir/ITEM_$i.txt") or $LOGGER->logdie($!);
#            &printall(\*ITEM, $gene{$chr}, \%feature, $name, $gffio1, \%RNA, {});
#            close(ITEM);
            
            &printall(\*GFFOUT, $gene{$chr}, \%feature, $name, $gffio1, {}, {gene=>{'Note'=>'LOCUS:'.$source_tag.'_'.$chr.'_'.sprintf("%05d",$i)} });
            $i++;
        }
        close(F1);
}
close(AUX);
close(GFFOUT);

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -g1     --gff1      GFF 1
        -g2     --gff2      GFF 2
        -g      --genome    Genome multi-fasta file
        -a      --auxfile   Auxiliary file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub recitem {
    my ($hr_merge, $hr_item, $item) = @_;
    unless (exists $hr_item->{$item}) {
        $hr_item->{$item} = undef;
        foreach my $it (keys %{ $hr_merge->{$item} }) {
            &recitem($hr_merge, $hr_item, $it, $hr_merge->{$item}->{$it} );
        }
        delete($hr_merge->{$item});
    }
}


sub printall {

    my ($FH, $hr_chr_gene, $hr_feature, $geneid, $gffio, $hr_RNA, $hr_extra, $nogene, $hr_RNAequal, $hr_CDSequal) = @_;
    

    unless ($nogene) {
        if (($hr_extra)&&($hr_extra->{'gene'})) {
            foreach my $k (keys %{ $hr_extra->{'gene'} }) {
                $hr_chr_gene->{ $geneid }->remove_tag($k) if ( $hr_chr_gene->{ $geneid }->has_tag($k) );
                $hr_chr_gene->{ $geneid }->add_tag_value($k, $hr_extra->{'gene'}->{$k});
            }            
        }

        $hr_chr_gene->{ $geneid }->source_tag($source_tag);
        print { $FH } $gffio1->gff_string( $hr_chr_gene->{ $geneid } ),"\n";
    }

    foreach my $mRNAfeat (@{ $hr_feature->{ $geneid }->{'child'}->{ 'mRNA' } } ) {

        my ($mRNAid) = $mRNAfeat->get_tag_values('ID');
        
        
        my $newmRNAid = $mRNAid;

        my $nomrna;
        my $nocds;

        if (($hr_RNAequal)&&(exists $hr_RNAequal->{$mRNAid})) {
            $LOGGER->logdie("Weird! More than one RNA equal ($mRNAid): ".join(';', keys %{ $hr_RNAequal->{$mRNAid} })) if ( scalar(keys %{ $hr_RNAequal->{$mRNAid} }) > 1 );
            my ($Requal) = keys %{ $hr_RNAequal->{$mRNAid} };
            # já foi mesclado anteriormente
            $nomrna = 1 if ($hr_RNAequal->{$mRNAid}->{$Requal});
            
            if (($hr_extra)&&($hr_extra->{'mRNA'})&&($hr_extra->{'mRNA'}->{'Note'})) {
                $hr_extra->{'mRNA'}->{'Note'}||='';
                $hr_extra->{'mRNA'}->{'Note'}.= $hr_extra->{'mRNA'}->{'Note'}.',';
            }
            $hr_extra->{'mRNA'}->{'Note'}.='Identical_RNA_from:'.$Requal;

            $mRNAfeat->remove_tag('ID');
            $newmRNAid = join('/',(sort ($mRNAid,$Requal))); # NEW mRNA id
            $mRNAfeat->add_tag_value('ID', $newmRNAid);
            
            $hr_RNAequal->{$mRNAid}->{$Requal} = 1;
            $hr_RNAequal->{$Requal}->{$mRNAid} = 1;
        }

        if (($hr_CDSequal)&&(exists $hr_CDSequal->{$mRNAid})) {
            $LOGGER->logdie("Weird! More than one CDS equal ($mRNAid): ".join(';', keys %{ $hr_CDSequal->{$mRNAid} })) if ( scalar(keys %{ $hr_CDSequal->{$mRNAid} }) > 1 );
            my ($Requal) = keys %{ $hr_CDSequal->{$mRNAid} };
            
            if (($hr_extra)&&($hr_extra->{'mRNA'})&&($hr_extra->{'mRNA'}->{'Note'})) {
                $hr_extra->{'mRNA'}->{'Note'}||='';
                $hr_extra->{'mRNA'}->{'Note'}.= $hr_extra->{'mRNA'}->{'Note'}.',';
            }
            $hr_extra->{'mRNA'}->{'Note'}.='Identical_CDS_from:'.$Requal;
            
            # já foi mesclado anteriormente
            $nocds = 1 if ($hr_CDSequal->{$mRNAid}->{$Requal});

            $hr_CDSequal->{$mRNAid}->{$Requal} = 1;
            $hr_CDSequal->{$Requal}->{$mRNAid} = 1;
        }

        $hr_RNA->{$newmRNAid}->{'gene'} = $geneid;

        if (($hr_extra)&&($hr_extra->{'mRNA'})) {
            foreach my $k (keys %{ $hr_extra->{'mRNA'} }) {
                $mRNAfeat->remove_tag($k) if ($mRNAfeat->has_tag($k));
                $mRNAfeat->add_tag_value($k, $hr_extra->{'mRNA'}->{$k});
            }            
        }
        
        $mRNAfeat->source_tag($source_tag);

        unless ($nomrna) {
            print { $FH } $gffio->gff_string( $mRNAfeat ),"\n";
        }

        my @fe;
        $LOGGER->logdie("Unable to find exon for $mRNAid") unless (exists $hr_feature->{$mRNAid}->{'child'}->{'exon'});
        if ($mRNAfeat->strand() == -1) {
            @fe = sort { $b->start() <=> $a->start() } @{ $hr_feature->{$mRNAid}->{'child'}->{'exon'} };
        }
        elsif ($mRNAfeat->strand() == 1) {
            @fe = sort { $a->start() <=> $b->start() } @{ $hr_feature->{$mRNAid}->{'child'}->{'exon'} };
        }
        foreach my $exonfeat ( @fe ) {
            $exonfeat->source_tag($source_tag);
            $exonfeat->remove_tag('Parent');
            $exonfeat->add_tag_value('Parent', $newmRNAid);
            unless ($nomrna) {
                print { $FH } $gffio->gff_string( $exonfeat ),"\n";
            }
        }
        
        $hr_RNA->{$mRNAid}->{'gene'} = $geneid;

        if (exists $hr_feature->{$mRNAid}->{'child'}->{'CDS'}) {
            my %eachCDS;
            {
                my @allcds = @{ $hr_feature->{$mRNAid}->{'child'}->{'CDS'} };
                foreach my $cf (@allcds) {
                    my ($CDSid);
                    # OGS and NCBI (each mRNA has only one CDS)
                    if ($cf->has_tag('ID')) {
                        ($CDSid) = $cf->get_tag_values('ID');
                    }
                    else {
                        $LOGGER->logdie("Not found an ID for CDS of $mRNAid");
                    }
                    my ($CDSgroup) = $CDSid=~/^([^#]+)/;
                    push(@{$eachCDS{ $CDSgroup }}, $cf);
                }
            }
            

            foreach my $CDSgroup (keys %eachCDS) {
                
                my @fc;
                if ($mRNAfeat->strand() == -1) {
                    @fc = sort { $b->start() <=> $a->start() } @{ $eachCDS{ $CDSgroup } };
                }
                elsif ($mRNAfeat->strand() == 1) {
                    @fc = sort { $a->start() <=> $b->start() } @{ $eachCDS{ $CDSgroup } };
                }
                else {
                    $LOGGER->logdie("Cannot found strand for $mRNAid");
                }
                foreach my $cdsfeat ( @fc ) {
                    $cdsfeat->source_tag($source_tag);
                    $cdsfeat->remove_tag('Parent');
                    $cdsfeat->add_tag_value('Parent', $newmRNAid);
                    unless ($nocds) {
                        print { $FH } $gffio->gff_string( $cdsfeat ),"\n";
                    }                        
                }
            }
        }
    }
}
 
sub recgroupname {
    my ($hr_result, $gene, $hr_groupname) = @_;
    unless (exists $hr_groupname->{$gene}) {
        $hr_groupname->{$gene} = undef;
        foreach my $w (keys %{$hr_result->{$gene}}) {
            &recgroupname($hr_result, $w, $hr_groupname);
        }
        delete($hr_result->{$gene});
    }                    
}

