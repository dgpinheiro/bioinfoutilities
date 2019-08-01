#!/usr/bin/env python2


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys, csv, StringIO, random, decimal, argparse

def printf(format, *args):
    sys.stdout.write(format % args)

parser = argparse.ArgumentParser(description='Generate variable-length single metagenomic fragments.')
parser.add_argument('-r', metavar='<reference_fasta>', dest="ref", help="Multi-FASTA file containing genomic sequences from which reads will be sampled.")
parser.add_argument('-a', metavar='<abundance_file>', dest="abund", help="Tab-delimited abundance file with an abundance value for each corre- sponding genome sequence in <reference fasta>")
parser.add_argument('-o', metavar='<output_file>', dest="output", help="Name for output file containing simulated uniform-length reads")
parser.add_argument('-t', metavar='<total_reads>', type=int, dest="total", help="The total number of fragments to sample from all genomes")
parser.add_argument('-i', metavar='<insert_mean_length>', type=int, dest="insert", default="0", help="Average length of inserts.")
parser.add_argument('-s', metavar='<insert_stddev>', type=int, dest="stddev", default="0", help="Standard deviation of insert length" )
parser.add_argument('-d', '--direction', action='store_true', dest="direction", help="Use this switch to generate reads in both forward and reverse orientations" )
args = parser.parse_args()


#Reference metagenome database file (FASTA)
f1 = open(args.ref);

#abundance file (tab-delimited .txt)
f2 = open(args.abund);

total_reads = args.total

insert_avg = args.insert
insert_stddev = args.stddev

if(insert_avg):
	f4 = open(args.output + '.1.fasta', 'w')
#	f5 = open(args.output + '.2.fasta', 'w')
else:
	f4 = open(args.output, 'w')

frags=[]

div_file = csv.reader(f2, delimiter='\t')
species=[]
diversity=[]

lengths=[]
freqs=[]
comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', 'W':'W', 'R':'Y', 'Y':'R', 'S':'S', 'K':'M', 'M':'K', 'B':'V', 'V':'B', 'D':'H', 'H':'D'}
for row in div_file:
	species.append(row[0][1:])
	diversity.append(decimal.Decimal(row[1]))


for i in SeqIO.parse(f1, 'fasta') :
	genome_num=0
	while(not(species[genome_num] in i.description)) :
		genome_num+=1
	if(species[genome_num] in i.description) :
		coverage=max(1, int((decimal.Decimal(diversity[genome_num])*total_reads)))

		limit=len(i.seq)
		for j in range(0, coverage) :
                	rand = random.random()
                	rand_length = 0
                	numLen = len(lengths)-1

			if( (insert_avg != 0) & (insert_stddev != 0)):
				cur_insert = int(random.gauss(insert_avg, insert_stddev))
				if(limit > cur_insert):
					start1 = random.randint(0, limit-cur_insert)
					end1 = start1 + cur_insert
				else:
					start1 = 0
					end1 = limit

				read1 = i.seq[start1:end1]
				if(args.direction):
					check = random.random()
					if(check < 0.5): #forward orientation
						f4.write(">%s\n" % i.description)
						f4.write("%s\n" % read1)
					else: #reverse orientation
						f4.write(">%s\n" % i.description)
						f4.write("%s\n" % read1[::-1])
				if(args.direction and random.random() < 0.5):
				        #reverse orientation
					f4.write(">%s\n" % i.description)
					f4.write("%s\n" % read1[::-1])
				else:
                                        #forward orientation
                                        f4.write(">%s\n" % i.description)
					f4.write("%s\n" % read1)

			else:
				if(limit > cur_insert) :
					start=random.randint(0, limit-max_read_length)
					end=start+max_read_length
				else:
					start=0
					end=limit
				read = i.seq[start:end]
				if(args.direction and random.random() < 0.5):
                                        #reverse orientation
					read = ''.join([comp[b] for b in i.seq[start:end][::-1]])
					f4.write(">%s\n" % i.description)
					f4.write("%s\n" % read)
					
				else:
                                        #forward orientation
					f4.write(">%s\n" % i.description)
					f4.write("%s\n" % read)

	if (genome_num >= len(species) ) :
		break;

f1.close()
f2.close()
f4.close()
