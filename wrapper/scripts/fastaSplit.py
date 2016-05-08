#!/usr/bin/env python

############################################################
### MODULES ################################################

import sys,random,re,math
from optparse import OptionParser
from bisect import *

############################################### /MODULES ###
############################################################




############################################################
### FUNCTIONS ##############################################

def read_fasta(filename):
	fd={}
	fl=[]
	sid=None
	if filename!='stdin':
		ff=open(filename,'r')
	else:
		ff=sys.stdin
	for line in ff:
		l=line[:-1]
		if l!='':
			if l[0]=='>':
				sid=l[1:]
				fd[sid]=[]
				fl.append(sid)
			else:
				fd[sid].append(l)	
		else:
			continue
	if filename!='stdin':
		ff.close()
	for sid in fd:
		fd[sid]=''.join(fd[sid])
	return fd,fl

	
		
def wrap(seq,cpl):
	l=len(seq)
	for x in range(0,l,cpl):
		print seq[x:x+cpl]
		
def padd(tot,x):
	return str(x).zfill(len(str(tot-1)))

############################################# /FUNCTIONS ###
############################################################




############################################################
### ARGUMENTS,OPTIONS ######################################

parser = OptionParser(usage="\n%prog [options]", version="%prog 0.1")
parser.add_option("-i", metavar="FILE", type="string", dest="input_file", default = 'stdin', help="input FASTA filename (default = 'stdin')")
parser.add_option("-s", metavar="STR", type="int", dest="splits", default = 1, help="splits (n). Default = 1")
parser.add_option("-o", metavar="FILE", type="string", dest="output_file", default = 'stdin', help="output FASTA tag (default = 'fastafile')")

(opt, args) = parser.parse_args()
if opt.input_file==None:
	parser.print_help()
	sys.exit(-1)
        
##################################### /ARGUMENTS,OPTIONS ###
############################################################




############################################################
### CONSTANTS ##############################################

############################################# /CONSTANTS ###
############################################################




############################################################
### MAIN ###################################################


fd, fl = read_fasta(opt.input_file)
sort_l=[]
for s in fl:
	insort(sort_l,(len(fd[s]),fd[s],s))


files=[]
totlen=[]
for x in xrange(opt.splits):
	totlen.append(0)
	files.append(open(opt.output_file+'_'+padd(opt.splits,x)+'.fa','w'))

for seq in sort_l[::-1]:
	i=totlen.index(min(totlen))
	totlen[i]+=seq[0]
	print>>files[i], '>'+seq[2]
	for x in range(0,len(seq[1]),60):
		print>>files[i], seq[1][x:x+60]
		

for x in xrange(opt.splits):
	files[x].close()
	
################################################## /MAIN ###
############################################################
