#!/usr/bin/env python

import numpy as np
import subprocess
import sys
import os
import sys

inputOptions = sys.argv[1:]

#usage: python search_virus_domains_new.py ../results/blastx/RNA24586.results.tab ../contigs/RNA24586.coverage_selected_contigs.fasta RNA24586


def main():

	sequences=read_fasta(inputOptions)
	sequences=read_blast(inputOptions,sequences)
	sequences=evaluate_sequences_coverage(sequences)

	for sequence in sequences.keys():
		if sequences[sequence].coverage>=0.3 and sequences[sequence].length >= 500:
			seq_info=inputOptions[2]+"|"+sequences[sequence].name+"|"+str(sequences[sequence].length)+"|"+str(sequences[sequence].cov)+"|"+str(round(sequences[sequence].coverage,2))
			
			virus_string=""
			for virus in sorted(sequences[sequence].hits, key=lambda Hit: Hit.identities, reverse=True):
				virus_string+=virus.name.split("|")[4].split("virus")[0]+"virus:"+str(round(virus.identities,2))+","
		
			print ">"+seq_info+"\t"+virus_string	
		 	print sequences[sequence].sequence


def read_fasta(inputOptions):

	fasta_file = [n for n in open(inputOptions[1],'r').read().replace("\r","").split("\n") if len(n)>0]
	
	
	sequences={}

	for line in fasta_file:
		if line[0:1]==">":
			name=line[1:].split(" ")[0]
			sequence=Sequence(name)
			sequence.length=int(line.split(" len_")[1].split("|")[0])
			sequence.cov=float(line.split("|rpb_")[1])
			sequences[name]=sequence
			
		else:
			sequences[name].sequence+=line

	return sequences


def read_blast(inputOptions,sequences):
	# # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	input_file = [n for n in open(inputOptions[0],'r').read().replace("\r","").split("\n") if len(n)>0]

	for line in input_file:
		if line.lower().find("phage")==-1 and line[0:1]!="#":

			hit=Hit(line.split("\t")[1])	
			
			if len(line.split("|len:"))>1:
				hit.length_prot=int(line.split("|len:")[1].split("|")[0])
			else:
				hit.length_prot=10		

			hit.alignment_length= int(line.split("\t")[3])
			hit.similarity=float(line.split("\t")[2])
			hit.identities=float(hit.alignment_length)*(float(hit.similarity)/100)

			start=int(line.split("\t")[6])
			end=int(line.split("\t")[7])
			if start>end:
				hit.start=end
				hit.end=start
			else:
				hit.start=start
				hit.end=end

			if hit.similarity >= 30 and float(hit.alignment_length)/float(hit.length_prot) >= 0.3:
				sequences[line.split("\t")[0]].hits.append(hit)

	return sequences



def evaluate_sequences_coverage(sequences):
	for sequence in sequences.keys():	
		covered={}
			
		for hit in sequences[sequence].hits:

			for i in range(hit.start,hit.end):
				covered[i]=1
		
		sequences[sequence].coverage= float(len(covered.keys()))/float(sequences[sequence].length)
		

	return sequences

class Sequence:
	def __init__(self, name):
		self.name=name
		self.sequence=""
		self.hits=list()
		self.length=0
		self.coverage=0
		self.cov=0

class Hit:
	def __init__(self, name):
		self.name=name
		self.simple_name=""
		self.alignment_length=0
		self.similarity=0
		self.identities=0
		self.length_prot=0
		self.start=0
		self.end=0
			



main()
