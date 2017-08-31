
import csv
import sys
import math
from Bio.Blast import NCBIXML


inputOptions = sys.argv[1:]

# usage: out_xml sequence.fasta  

def main():

	host_homology={}
	input_file = [n for n in open(inputOptions[1],'r').read().replace("\r","").split("\n") if len(n)>0]
	sequences={}
	info={}
	queries={}
	blast_result={}
	for line in input_file:
		if line[0:1]==">":
			name = line.split("\t")[0][1:]
			sequences[name]=""		
			info[name]=line[1:]
			queries[name]=0
			host_homology[name]=bool(0)

			blast_result[name]=list()


			hit=("none",0)
			blast_result[name].append(hit)

		else:
			sequences[name]+=line


		        

	blast_file = [n for n in open(inputOptions[0],'r').read().replace("\r","").split("\n") if len(n)>0]


	with open(inputOptions[0],'r') as f:
    		for raw_line in f:
			line=raw_line.replace("\n","").replace("\r","")

			query=line.split("\t")[0]

			align_length=int(line.split("\t")[5])
			hit_def=line.split("\t")[2]
			identities=int(line.split("\t")[7])
			title=line.split("\t")[1]
	
			if align_length>=80 and hit_def.lower().find(inputOptions[2].replace("_"," ").lower())!=-1:
				host_homology[query]=bool(1)
			
			hit=(title+hit_def, str(round(float(identities)/float(len(sequences[query])),4)))

			blast_result[query].append(hit)

			blast_result[query]=sorted(blast_result[query], key=lambda tup: tup[1],reverse=True)[0:1]






	for query in queries.keys():

		to_print=""
		for i in blast_result[query]:
			to_print+="\t"+str(i[0])+"\t"+str(i[1])

		print info[query].replace("|","\t") + to_print+"\t"+str(host_homology[query])+"\t"+sequences[query]





	

main()		
