#!/usr/bin/python 

################################################################################
### MODULES ####################################################################

import sys, re, math, random, os
from commands import getoutput
from optparse import OptionParser

tmpdir=os.getenv('TMPDIR')
if tmpdir == None :
	tmpdir="/tmp"


################################################################### /MODULES ###
################################################################################




################################################################################
### FUNCTIONS ##################################################################

def fasta_seq(sid, seq, cpl):
	print '>' + sid
	for x in range(0, len(seq), cpl):
		print seq[x:x+cpl]

def read_fasta(ff):
	fd = {}
	fl = []
	sid = None
	for l in ff:
		if l != '':
			if l[0] == '>':
				sid = l[1:]
				fd[sid] = []
				fl.append(sid)
			else:
				fd[sid].append(l.upper())	
		else:
			continue
	for sid in fd:
		fd[sid] = ''.join(fd[sid])
	return fd, fl

def read_alignment(alignment):
	fd = {}
	fl = []
	sid = None
	for l in alignment.split('\n'):
		if l != '':
			if l[0] == '>':
				sid = l[1:]
				fd[sid] = []
				fl.append(sid)
			else:
				fd[sid].append(l.upper().replace('X','-'))	
		else:
			continue
	for sid in fd:
		fd[sid] = ''.join(fd[sid])
	return fd, fl


def gap_clean(seq, side):
	if side == '5' or side == 'B':
		for m in re.finditer('^-+', seq):
			ms = m.span()
			seq = '*'*(ms[1] - ms[0]) + seq[ms[1]:]
	if side == '3' or side == 'B':
		for m in re.finditer('-+$', seq):
			ms = m.span()
			seq = seq[:ms[0]] + '*'*(ms[1] - ms[0])
	return seq


def get_consensus(temp_file_name,to_align):
	alpha = ['A', 'C', 'G', 'T', 'N', '-']
	od, ol = read_fasta(to_align.split('\n'))
	nseqs = len(ol)

	if nseqs == 1: 
		if ol[0][-1]!= 'N':
			sys.stderr.write('only one read '+ol[0]+'\n')
		return od[ol[0]]

	#elif nseqs == 2 and (ol[0][-1]== 'N' or ol[1][-1]== 'N'): # both equally consensual
	#	if ol[0][-1]== 'N':
	#		if ol[1][-1]== 'N':
	#			return od[random.choice(ol)]
	#		else:
	#			return od[ol[0]]
	#	else:
	#		if ol[1][-1]== 'N':
	#			return od[ol[1]]
	
	elif nseqs >= 2: # align, build mm0, score, pick, if more than one: random choice
		ta = ''
		param_open = '1.2'
		count5 = 0
		count3 = 0
		countN = 0
		for sid in ol:
			ta += '>'+sid+'\n'+od[sid]
			if od[sid] == '':
				ta += 'X'
			ta += '\n'
		for head in ol:
			if head[-1] == 'N':
				countN += 1
			elif head[-1] == '5':
				param_open = '2'
				count5 += 1
			elif head[-1] == '3':
				param_open = '2'
				count3 += 1
			else :
				param_open ='2'
				count3 +=1
				count5 +=1
		if countN == 0 and (count5 == 0 or count3 == 0) :
			sys.stderr.write('only one read '+ol[0]+'\n')
		temp_file =  tmpdir + "/" + temp_file_name + ".tmp"
		ta_file=open(temp_file,'w')
		print>>ta_file, ta
		ta_file.close()
		#alignment = getoutput('printf \''+ta+'\'| mafft.bat --retree 1 --quiet --op '+param_open+' -')
		alignment = getoutput('cat ' + temp_file +' | mafft --retree 1 --quiet --op '+param_open+' -')
		os.remove(temp_file)
		fd, fl = read_alignment(alignment)
		al = len(fd[fl[0]])
		model = {}
		for c in alpha:
			model[c] = [0]*al
		tts = [float(0)]*al
		max_sco = [(None,None)]
	
		for sid in fl:
			fd[sid] = gap_clean(fd[sid], sid[-1])
			for x in xrange(al):
				c = fd[sid][x]
				if c in model:
					model[c][x] += 1
					tts[x] += 1
	
	    	# Commented so that always consensus is made
		#for sid in fl:
			#sco = '-Inf' ### no need
			#if sid[-1] == 'N':
			#	sco = 0
			#	for x in xrange(al):
			#		obs = model[fd[sid][x]][x]/tts[x]
			#		exp = float(1)/len(alpha)
			#		sco += math.log(obs/exp, 2)
			#	if sco > max_sco[0][1]:
			#		max_sco = [(sid, sco)]
			#	elif sco == max_sco[0][1]:
			#		max_sco.append((sid, sco))
			
			#print sid+'\t'+fd[sid]+'\t'+str(sco) ### no need
		#print ### no need
			
		if max_sco != [(None,None)]:
			consensus = fd[random.choice(max_sco)[0]]
		else:
			consensus = ''
			for x in xrange(al):
				maxc = 'A'
				for c in alpha:
					if model[c][x] > model[maxc][x]:
						maxc = c
				consensus = consensus + maxc
		return consensus.replace('-', '')

	else: # nseqs == 0 should raise error
		sys.stderr.write('no sequence \n')
		return ''
	


################################################################# /FUNCTIONS ###
################################################################################




################################################################################
### ARGUMENTS,OPTIONS ##########################################################

parser = OptionParser(usage = "\n%prog [options]", version="%prog 0.1")

parser.add_option(
	"-i",
	metavar = "FILE",
	type = "string",
	dest = "input_file",
	default = 'stdin',
	help = "alignment FASTA filename (default = 'stdin')"
	)

parser.add_option(
        "-f",
        metavar = "FILE",
        type = "string",
        dest = "temp_file",
        default = 'out',
        help = "alignment temporal prefix FASTA filename (default = 'out')"
        )
(opt, args) = parser.parse_args()

######################################################### /ARGUMENTS,OPTIONS ###
################################################################################




################################################################################
### CONSTANTS ##################################################################

################################################################# /CONSTANTS ###
################################################################################




################################################################################
### MAIN #######################################################################
			
if __name__ == "__main__":
		
	if opt.input_file != 'stdin':
		ff = open(opt.input_file, 'r')
	else:
		ff = sys.stdin

	bid = None
	bseq = ''
	to_align = ''
	
	for line in ff:
		l = line.strip()
		if l:
			if l[0] == '>':
				if l[1] != '>':
					if to_align:
						bseq += get_consensus(opt.temp_file + bid, to_align)
						to_align = ''
					if bid:
						fasta_seq(bid,bseq,60)
					bid = l[1:]
					#sys.stderr.write(bid)
					bseq = ''
				else: 
					if to_align:
						bseq += get_consensus(opt.temp_file + bid, to_align)
						to_align = ''
		
			elif l[0] == '<':
				if to_align:
					to_align += '\n'
				to_align += l.replace('<', '>') + '\n'
			
			else:
				if to_align:
					to_align += l
				else:
					bseq += l
	if to_align:
		bseq += get_consensus(opt.temp_file + bid, to_align)
	if bseq:
		fasta_seq(bid,bseq,60)
			 		
	if opt.input_file != 'stdin':
		ff.close()

###################################################################### /MAIN ###
################################################################################
