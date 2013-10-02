#!/usr/bin/perl

use strict;
use warnings; 

use Bio::SeqIO; 
use Bio::Tools::GFF; 
use Bio::Seq; 

die "Usage: $0 <genome file> <promoter size> " if ($#ARGV==-1);

my $genomefile = $ARGV[0];
my $promotersize = $ARGV[1];

die "Missing genome multi fasta file" unless ($genomefile);
die "Wrong genome multi fasta file ($genomefile)" unless (-e $genomefile);

die "Missing promoter size" unless (defined $promotersize);
die "Wrong promoter size ($promotersize)" unless ($promotersize =~ /^\d+$/);

my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 3); 

my %gene; 

while (my $feat = $gffio->next_feature()) {
	die "Sorry... only gene features are accepted!" if ($feat->primary_tag() ne 'gene');
	push(@{$gene{$feat->seq_id()}}, $feat); 
} 

$gffio->close(); 

my $promoter=Bio::SeqIO->new(-fh=>\*STDOUT,-format=>"FASTA"); 

my $genome=Bio::SeqIO->new(-file=>$genomefile,-format=>"FASTA"); 

while(my $seq=$genome->next_seq()){ 
	if (exists $gene{ $seq->display_id() }) { 
		foreach my $feat (@{ $gene{$seq->display_id()} }) { 

			my ($geneid)=$feat->get_tag_values("ID"); 
			my ($start, $end);
			if ($feat->strand() == 1) { 
				next if ($feat->start() == 1);
				$start=1; 
				if ($feat->start() > $promotersize) { 
					$start=$feat->start()-$promotersize;  
				} 
				$end=($feat->start()-1);
			} else {
				next if ($feat->end() == $seq->length()); 
				$end=$seq->length(); 
				if ( ($seq->length()-$feat->end) > $promotersize ) { 
					$end=$feat->end()+$promotersize;
				}
				$start=$feat->end()+1;
			}   
			my $tmp = $seq->trunc($start, $end); 
			$tmp->display_id($geneid); 
			$promoter->write_seq($tmp);   
		}   
	}  
}
