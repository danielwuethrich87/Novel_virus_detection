#!/usr/bin/env python

import numpy as np
import subprocess
import sys
import os
import sys
import gzip

inputOptions = sys.argv[1:]

#usage: file1


def main():




	info={}
	info["sequence"]=""
	info["prot_ID"]=""
	info["acession_number"]=""
	info["occurance"]=""
	info["organism"]=""
	info["sequence_length"]=0
	info["protein_name"]=""

	with gzip.open(inputOptions[0],'r') as f:
    		for line in f:
			construct_fasta(line.replace("\n",""),info)
			

def construct_fasta(line,info):
	if line[0:2]=="DE" and info["protein_name"]=="":
		if line.find("Full=")!=-1:
			info["protein_name"]+=line.split("Full=")[1].split("{")[0]
	if line[0:2]=="ID":
		info["prot_ID"]+=line.split("   ")[1]
	if line[0:2]=="AC":
		info["acession_number"]+=line.split("   ")[1]
	if line[0:2]=="OC":
		info["occurance"]+=line.split("   ")[1]
	if line[0:2]=="OS":
		info["organism"]+=line.split("   ")[1]
	if line[0:2]=="SQ":
		info["sequence_length"]+=int(line.split("   ")[2].split(" ")[0])
	if line[0:5]=="     ":
		info["sequence"]+=line.replace(" ","")
	if line[0:2]=="//":
		assert(info["sequence_length"]==len(info["sequence"]))
			
		print ((">"+info["prot_ID"]+"|"+info["acession_number"]+"|len:"+str(info["sequence_length"])+"|"+info["protein_name"]+"|"+info["organism"]).replace(" ","_")+" "+info["occurance"])[0:999]

		i=0
		print_string=""
		while (i<=len(info["sequence"])):
			print_string+= info["sequence"][i:i+60]+"\n"
			i=i+60
		assert(info["sequence_length"]==len(print_string.replace("\n","")))

		print print_string

		info["sequence"]=""
		info["prot_ID"]=""
		info["acession_number"]=""
		info["occurance"]=""
		info["organism"]=""
		info["sequence_length"]=0
		info["protein_name"]=""
main()
