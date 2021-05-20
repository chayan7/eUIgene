#Finding intergenic regions between query gene and the preceding gene
#Author: Chayan Kumar Saha
#This tool uses an output file from FlaGs (that have a suffix *_operon.tsv) as input file
#FlaGs can be downloaded from: https://github.com/GCA-VH-lab/FlaGs



from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC
import random
import argparse
import ftplib
import os, sys, os.path, math
import gzip
import subprocess
import glob
import time
import textwrap


usage= ''' Description:  Get intergenic regions between query gene and upstream genes using FlaGs outfile that has suffix "_operon.tsv" '''




parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-ot", "--operontsv", help="Operon tsv file")
parser.add_argument("-sp", "--sposition", help="300 or 700, Upstream is extracted based on the position input, but within the same contig, even if it encroaches in upstream genes. Default is intergenic nucleotides between preceding gene end to query gene start")
parser.add_argument("-ip", "--iposition", help="300 or 700, Upstream is extracted only from the intergenic region. Default is intergenic nucleotides between preceding gene end to query gene start")
parser.add_argument("-o", "--out_prefix", required= True, help="Any Keyword to define your output eg. MyQuery")
parser.add_argument("-k", "--keep", action="store_true", help="If you want to keep the intermediate files eg. gff3 use [-k]. By default it will remove.")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.1')


args = parser.parse_args()
parser.parse_args()




'''
input_file
WP_082602131.1#2_Phenylobacterium_sp._Root700   537     -       -       0       -2241   -1705   121099  121635  WP_056735442.1#2        NZ_LMHT01000030.1       GCF_001429025.1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/429/025/GCF_001429025.1_Root700
WP_082602131.1#2_Phenylobacterium_sp._Root700   576     -       -       0       -1663   -1088   120482  121057  WP_056735441.1#2        NZ_LMHT01000030.1       GCF_001429025.1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/429/025/GCF_001429025.1_Root700
WP_082602131.1#2_Phenylobacterium_sp._Root700   549     -       +       0       -1110   -562    119956  120504  WP_082602138.1#2        NZ_LMHT01000030.1       GCF_001429025.1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/429/025/GCF_001429025.1_Root700
'''

refDb='./refSeq.db'
genDb='./genBank.db'
if os.path.isfile(refDb):
	refDbSize=os.path.getsize(refDb)
else:
	refDbSize='0'
if os.path.isfile(genDb):
	genDbSize=os.path.getsize(genDb)
else:
	genDbSize='0'

ftp = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
ftp.cwd("/genomes/refseq") # move to refseq directory

filenames = ftp.nlst() # get file/directory names within the directory
if 'assembly_summary_refseq.txt' in filenames:
	ftp.sendcmd("TYPE i")
	if int(ftp.size('assembly_summary_refseq.txt'))!=int(refDbSize):#check if the previously downloaded db exists and if that's updated to recent one
		ftp.retrbinary('RETR ' + 'assembly_summary_refseq.txt', open('refSeq.db', 'wb').write) # get the assembly summary from refseq
	else:
		pass

ftp_gen = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
ftp_gen.cwd("/genomes/genbank") # move to refseq directory

filenames = ftp_gen.nlst() # get file/directory names within the directory
if 'assembly_summary_genbank.txt' in filenames:
	ftp_gen.sendcmd("TYPE i")
	if int(ftp_gen.size('assembly_summary_genbank.txt'))!=int(genDbSize):#check if the previously downloaded db exists and if that's updated to recent one
		ftp_gen.retrbinary('RETR ' + 'assembly_summary_genbank.txt', open('genBank.db', 'wb').write) # get the assembly summary from refseq
	else:
		pass

assemblyName={}
bioDict={} #bioproject as keys and assemble number (eg.GCF_000001765.1) as value
accnr_list_dict={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
with open('refSeq.db', 'r') as fileIn:
	for line in fileIn:
		if line[0]!='#':
			Line=line.rstrip().split('\t')
			accnr_list_dict[Line[0]]= Line[7]+'\t'+Line[19]
			bioDict[Line[1]]=Line[0]
			assemblyName[Line[0]]=Line[0]

ftp_gen = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
ftp_gen.cwd("/genomes/genbank") # move to refseq directory


assemblyName_GCA={}
bioDict_gen={}
accnr_list_dict_gen={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
with open('genBank.db', 'r') as fileIn:
	for line in fileIn:
		if line[0]!='#':
			Line=line.rstrip().split('\t')
			if len(Line)>19:
				if Line[18]=='identical':
					if Line[17] in accnr_list_dict:
						bioDict_gen[Line[1]]=Line[0]
						accnr_list_dict_gen[Line[0]]= accnr_list_dict[Line[17]]
						assemblyName_GCA[Line[0]]=Line[17]
					else:
						accnr_list_dict_gen[Line[0]]=Line[7]+'\t'+Line[19]
				else:
					accnr_list_dict_gen[Line[0]]=Line[7]+'\t'+Line[19]

bioDict.update(bioDict_gen)
accnr_list_dict.update(accnr_list_dict_gen)
assemblyName.update(assemblyName_GCA)

print ('\n'+ '>> Database Downloaded. Cross-checking of the accession list in progress ...'+ '\n')




acc_sp={}
acc_ftp={}
acc_chrm={}
acc_gcf={}


operontsvIn=open(args.operontsv, 'r').read().rstrip()
splittedOTsv=operontsvIn.split("\n\n")
positionDict={}
strandDict={}
count=set()
for item in splittedOTsv:
	if item!='':
		start=[]
		end=[]
		fgenesList=[]
		queryaccName=''
		strandOrder=''
		for items in item.split('\n'):
			accs=items.split('\t')[0].split('|')[0]#+'_'+items.split('\t')[0].split('_')[1]
			acc_sp[accs]=items.split('\t')[0].split('|')[1]
			if items.split('\t')[11] in accnr_list_dict: #if gcf in accnr_list_dict
				#print(items.split('\t')[11], accnr_list_dict[items.split('\t')[11]])
				acc_ftp[accs]=accnr_list_dict[items.split('\t')[11]].split('\t')[1] #ftp link
				acc_chrm[accs]=items.split('\t')[10]
				acc_gcf[accs]=items.split('\t')[11]
				startP=int(items.split('\t')[7])
				endP=int(items.split('\t')[8])
				fgenes=items.split('\t')[9]
				fgenesList.append(fgenes)
				queryaccName=items.split('\t')[0].split('|')[0]
				strandOrder=items.split('\t')[2]
				start.append(startP)
				end.append(endP)
			else:
				count.add(items.split('\t')[11])
				acc_ftp[accs]='not_found'

		#if strandOrder=='+':
		#print(accs, start, end)
		#print(accs,fgenesList, start, end)
		if args.sposition:
			if fgenesList:
				if fgenesList.index(queryaccName)!=0:
					if strandOrder=='+':
						positionDict[accs]=str(start[fgenesList.index(queryaccName)]-int(args.sposition))+'\t'+str(start[fgenesList.index(queryaccName)]-1)
						strandDict[accs]='+'
					if strandOrder=='-':
						positionDict[accs]=str(end[fgenesList.index(queryaccName)]+1)+'\t'+str(end[fgenesList.index(queryaccName)]+int(args.sposition))
						strandDict[accs]='-'
		else:
			if fgenesList:
				if strandOrder=='+':
					#print(accs, len(end), fgenesList.index(queryaccName)+1, fgenesList, queryaccName)
					if fgenesList.index(queryaccName)!=0:
					#if end[fgenesList.index(queryaccName)-1]:
						positionDict[accs]=str(end[fgenesList.index(queryaccName)-1]+1)+'\t'+str(start[fgenesList.index(queryaccName)]-1)
						strandDict[accs]='+'
				if strandOrder=='-':
					#print(accs, len(end), fgenesList.index(queryaccName), fgenesList, queryaccName)
					if fgenesList.index(queryaccName)!=0:
					#if end[fgenesList.index(queryaccName)+1]:
						positionDict[accs]=str(end[fgenesList.index(queryaccName)]+1)+'\t'+str(start[fgenesList.index(queryaccName)-1]-1)
						strandDict[accs]='-'
			#if strandOrder=='-': #[2223371:2223829]
				#positionDict[accs]=str(end[fgenesList.index(queryaccName)])+'\t'+str(start[fgenesList.index(queryaccName)-1]-1)

print('Not found = '+ str(len(count)))

#print(acc_sp) #'WP_089251833.1#7': 'Rhodococcus_kyotonensis'
#print(acc_ftp) #'WP_089251833.1#7': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/188/125/GCF_900188125.1_IMG-taxon_2675903021_annotated_assembly'
#print(acc_chrm) #'WP_089251833.1#7': 'NZ_FZOW01000023.1'
#print(acc_gcf) #'WP_089251833.1#7': 'GCF_900188125.1'
#print(positionDict) #'WP_089251833.1#7': '27226\t35524'

def revComp(sp, seqs):
	if sp=='+':
		return seqs
	if sp=='-':
		if 'U' not in seqs:
			complementSeq=seqs[::-1]
			revcomplementSeq=complementSeq.translate(str.maketrans('ATGC','TACG'))              # reverse complementery sequence for ACGT
			return revcomplementSeq
		elif 'T' not in seqs:
			complementSeq=seqDict[seqid][::-1]
			revcomplementSeq=complementSeq.translate(str.maketrans('AUGC','UACG'))              # reverse complementery sequence for ACGU
			return revcomplementSeq
#import sys
#sys.exit()
seqDict={}
spDict={}
for acc in acc_ftp:
	if acc_ftp[acc]!='not_found':
		if os.path.isfile(acc_gcf[acc]+'.gfna.gz'):
			Line=gzip.open(acc_gcf[acc]+'.gfna.gz', 'rb').read().decode('utf-8').split('>') #Download and read gff.gz
			if Line:
				for item in Line:
					#print(item.split(' ')[0])
					if item.split(' ')[0]==acc_chrm[acc]:
						#print(item.split('\n')[0])
						seqId=item.split(' ')[0]
						#print(seqId)
						#spInfo='_'+item+'_'+accnr_list_dict[item].split('\t')[0]
						#print(spInfo)
						sequence=('').join(item.split('\n')[1:])
						seqDict[seqId]=sequence
						#spDict[seqId]=spInfo
					else:
						pass
			else:
				pass
		if not os.path.isfile(acc_gcf[acc]+'.gfna.gz'):
			ftpLine=acc_ftp[acc]
			print(ftpLine)
			ftpsplitDir = ftpLine.split('/')[3:]
			ftp_path = '/'.join(map(str,ftpsplitDir))
			ftp = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
			ftp.cwd('/'+ftp_path)
			files = ftp.nlst()
			for elements in files:
				if ftpsplitDir[-1]+'_genomic.fna.gz' in elements:
					#print(elements) #Check if GFF.gz is there
					#try:
					ftp.retrbinary('RETR ' + elements, open(acc_gcf[acc]+'.gfna.gz', 'wb').write)
					print('Download complete: '+elements)
					Line=gzip.open(acc_gcf[acc]+'.gfna.gz', 'rb').read().decode('utf-8').split('>') #Download and read gff.gz
					if Line:
						for item in Line:
							#print(item.split(' ')[0])
							if item.split(' ')[0]==acc_chrm[acc]:
								#print(item.split('\n')[0])
								seqId=item.split(' ')[0]
								#print(seqId)
								#spInfo='_'+item+'_'+accnr_list_dict[item].split('\t')[0]
								#print(spInfo)
								sequence=('').join(item.split('\n')[1:])
								seqDict[seqId]=sequence
								#spDict[seqId]=spInfo
							else:
								pass
					else:
						pass
				else:
					pass
	else:
		pass



if args.keep:
	pass
else:
	subprocess.Popen("rm *.db", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	subprocess.Popen("rm G*F*.gz", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


with open(args.out_prefix+'_upstreamIntergenic.txt', 'w') as ftab:
	for acc in acc_chrm:  #{'WP_042359506.1': 'NZ_BAXL01000009.1', 'WP_088005586.1': 'NZ_NHTN01000003.1
		if acc_chrm[acc]!='not_found':
			if acc in positionDict:
				fromPos=str(positionDict[acc].split('\t')[0])
				toPos=str(positionDict[acc].split('\t')[1])
				strandPos=strandDict[acc]
				if int(fromPos)>0 and int(toPos)>0:
					if acc_chrm[acc] in seqDict:
						if int(fromPos)<len(seqDict[acc_chrm[acc]])+1 and int(toPos)<len(seqDict[acc_chrm[acc]])+1 and int(fromPos)<int(toPos):
							print('>'+str(args.out_prefix)+'_'+acc_chrm[acc]+'_'+acc+'_'+acc_sp[acc]+'_'+acc_gcf[acc]\
							+';length='+str(len(seqDict[acc_chrm[acc]][int(fromPos)-1:int(toPos)]))+'/'+str(len(seqDict[acc_chrm[acc]]))+\
							'['+fromPos+':'+toPos+']'+'\n'+revComp(strandPos, seqDict[acc_chrm[acc]][int(fromPos)-1:int(toPos)]),file=ftab)
						else:
							if int(fromPos)<len(seqDict[acc_chrm[acc]])+1 and int(fromPos)<int(toPos):
								if not int(toPos)<len(seqDict[acc_chrm[acc]])+1:
									print('>'+str(args.out_prefix)+'_'+acc_chrm[acc]+'_'+acc+'_'+acc_sp[acc]+'_'+acc_gcf[acc]\
									+';length(adj)='+str(len(seqDict[acc_chrm[acc]][int(fromPos)-1:]))+'/'+str(len(seqDict[acc_chrm[acc]]))+\
									'['+fromPos+':'+toPos+']'+'\n'+revComp(strandPos, seqDict[acc_chrm[acc]][int(fromPos)-1:]),file=ftab)

					else:
						pass
				else:
					print('Wrong input, please check again')
			else:
				pass



if (glob.glob(args.out_prefix+'_upstreamIntergenic.txt')):
	fasIn=open(args.out_prefix+'_upstreamIntergenic.txt', 'r').read().rstrip().split('>')
	if args.iposition:
		with open(args.out_prefix+'_upstreamIntergenic.fasta', 'w') as fileOut:
			for items in fasIn:
				if items!='':
					if len(items.split('\n')[1])>=int(args.iposition):
						fastaID='>'+items.split('\n')[0].split(';')[0]+';length='+args.iposition
						sequence=textwrap.fill(items.split('\n')[1][(int(args.iposition)+1*-1):],80)
						print(fastaID,sequence,sep='\n', file=fileOut)
					else:
						fastaID='>'+items.split('\n')[0].split(';')[0]+';length='+str(len(items.split('\n')[1]))
						sequence=textwrap.fill(items.split('\n')[1],80)
						print(fastaID,sequence,sep='\n', file=fileOut)

	if not args.iposition:
		with open(args.out_prefix+'_upstreamIntergenic.fasta', 'w') as fileOut:
			for items in fasIn:
				if items!='':
					fastaID='>'+items.split('\n')[0]
					sequence=textwrap.fill(items.split('\n')[1],80)
					print(fastaID,sequence,sep='\n', file=fileOut)

if (glob.glob(args.out_prefix+'_upstreamIntergenic.fasta')):
	subprocess.Popen("rm %s"%(args.out_prefix+'_upstreamIntergenic.txt'), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

print('\n\nCompleted!\n')
sys.exit()
