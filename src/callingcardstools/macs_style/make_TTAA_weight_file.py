"""
make_TTAA_weight_file.py 

written 11/17/16 by RDM
modified 7/27/17 to work with ccf format.
 

This file takes in a genome (e.g. mm10 or Hg19) and gives the locations of all of the 
TTAA's.  Since TTAA is palindromic, it will only report the positions on the sense strand.
It will report positions on all numbered chromosomes and X, Y, M, but not 1_GL456210_random
for example.  

Usage 
python make_TTAA_weight_file.py -g <genome path and filename>  -o <output path and filename>

"""
from Bio import SeqIO
from Bio import Seq
import argparse
import re
import pandas as pd


def make_TTAA_weight_file(inputfile,outputfile):

	#open files

	in_filehandle = open(inputfile,'rU')
	out_filehandle = open(outputfile,'w')

	#read in first fasta record
	for seq_record in SeqIO.parse(in_filehandle,"fasta"):
		match_flag = 0
		#is it a numbered chromosome or X,Y,M?
		match_test = re.match(r"^chr\d+$",seq_record.name,re.I)
		if match_test:
			match_flag = 1
		else:
			match_test = re.match(r"^chrM$",seq_record.name,re.I)
			if match_test:
				match_flag = 1
			else:
				match_test = re.match(r"^chrX$",seq_record.name,re.I)
				if match_test:
					match_flag = 1
				else:
					match_test = re.match(r"^chrY$",seq_record.name,re.I)
					if match_test:
						match_flag = 1
		if match_flag:
			chr = seq_record.name
			#if so, find all TTAAs and print their locations to output files

			#convert to upper case
			seq_record.seq = seq_record.seq.upper()
			start_coor_zi = seq_record.seq.find("TTAA")
			while start_coor_zi != -1:
				#output coord, start
				print >> out_filehandle,chr+"\t"+str(start_coor_zi+1)+"\t"+str(start_coor_zi+4)
				start_coor_zi = seq_record.seq.find("TTAA",start_coor_zi+4)
	out_filehandle.close()
	TTAA_frame = pd.read_csv(outputfile,delimiter='\t',header=None,names=['Chr','Start','End'])
	chrlist = ['chr1','chr2','chr3','chr4','chr5',
		'chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
		'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21',
		'chr22','chrX','chrY','chrM']	
	TTAA_frame["Chr"] = pd.Categorical(TTAA_frame["Chr"],chrlist)
	TTAA_frame = TTAA_frame.sort_values(['Chr','Start','End'])
	TTAA_frame.to_csv(outputfile,sep="\t",header = False, index = False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='cc_filter_reads.py')
    parser.add_argument('-g','--genome',help='genome fasta file path and name',required=True)
    parser.add_argument('-o','--output',help='output filename (full path)',required=True)
    args = parser.parse_args()
    make_TTAA_weight_file(args.genome,args.output)


	

	


