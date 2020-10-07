#!/usr/bin/env python

import numpy as np, sys, os, random
import formula_processing as fp

'''
ARGS: 1) file path(s) to input mgf(s). If multiple, must be separated by
commas. 2) Metabolite class names of each mgf 3) binwidth 4) output file path
'''

def onehot(lin, bigl):
	'''
	##lin is a list with a sublist of MS/MS peaks
	##bigl is the list of ranges for encoding
	##RETURNS: lin with the sublist of peaks
	##one-hot encoded.
	'''
	peaks = lin[2]
	ohpeaks = []
	for ran in range(len(bigl)-1):
		st = float(bigl[ran])
		end = float(bigl[ran + 1])
		if (any(peak >= st and peak <= end for peak in peaks)):
			ohpeaks.append(1)
		else: ohpeaks.append(0)

	lin[2] = ohpeaks
	return lin

def parse_args(mgf, nam):
	'''
	Both mgf and nam are strings
	This function makes the next one more efficient
	'''
	mgfs = mgf.split(',')
	nams = nam.split(',')
	if len(mgfs) == len(nams):
		return dict(zip(nams, mgfs))
	else:
		sys.exit('Error: # of mgfs does not equal # of names!')
		
def parse_mgf(mgfnamd):
	'''
	mgf is a dict, output of parse_args
	Returns: dict of entries of all mgfs,
	Keys are the class names
	'''
	print(mgfnamd)
	outd = {}
	
	for nam, path in mgfnamd.items():
		mgf = open(path, 'r')
		mgfl = mgf.readline()
		outd[nam] = []
		ct = 0
		while mgfl:
			if mgfl.startswith('BEGIN IONS'):
				mgfl = mgf.readline()
				name = str(ct) + '_' + str(mgfl.split('=')[1])
				outd[nam].append([name, nam])
				ct += 1
			if mgfl.startswith('ION'):
				peaks = []
				mgfl = mgf.readline()
				while mgfl[1].isdigit():
					##we are just logging mz, not abund
					peaks.append(float(mgfl.strip().split()[0]))
					mgfl = mgf.readline()
				outd[nam][-1].append(peaks)
			mgfl = mgf.readline()
		
		mgf.close()
		print('Number of entries: ', str(ct))
	return outd
		
def make_output(mgfs, binw, outp):
	'''
	mgfs is a dictionary (output of parse_mgf)
	binw is a float
	outp is a string (output path)
	Returns: nothing
	Outputs: training and test files for iRF
	'''
	##making list of bins
	biglist=np.arange(50,2500,float(binw))

	##getting balanced dictionary (same number of instances for each class)
	lenl = [len(i) for i in mgfs.values()]
	minins = min(lenl)
	for nam, ins in mgfs.items():
		if len(ins) > minins:
			mgfs[nam] = random.sample(ins, minins)

	##reserving instances for the test arff
	rlst = []
	trlst = []
	for insl in mgfs.values():
		if len(insl) >= 10:
			sampnum = int(float(0.1) * len(insl))
		else:
			sampnum = 1
		tmp = random.sample(insl, sampnum)

		##now getting remainder
		remainder = insl
		for x in tmp:
			remainder.remove(x)

		##getting instances not in tmp (these are the training ones)
		rlst.append(tmp)
		trlst.append(remainder)
	
	flatr = [item for sublist in rlst for item in sublist]
	trainr = [item for sublist in trlst for item in sublist]
	print('done with making lists!')

	##initializing output files
	trainx = open(outp + 'xtrain.tab', 'w')
	testx = open(outp + 'xtest.tab', 'w')
	trainy = open(outp + 'ytrain.tab', 'w')
	testy = open(outp + 'ytest.tab', 'w')
	##output files for keeping track of ids of instances
	trainid = open(outp + 'idtrain.tab', 'w')
	testid = open(outp + 'idtest.tab', 'w')

	##writing out feature names (header) for x-files
	news = ''
	for i in range(len(biglist)-1):
		news += f'{biglist[i]}-{biglist[i+1]}\t'

	news = news[:-1]
	trainx.write(f'{news}\n')
	testx.write(f'{news}\n')

	##doing one-hot now and writing out
	##first for test arff instances
	for ind in range(len(flatr)):
		flatr[ind] = onehot(flatr[ind], biglist)
		featstrr = '\t'.join(str(x) for x in flatr[ind][2])
		testx.write(f'{featstrr}\n')
		testy.write(f'{flatr[ind][1]}\n')
		testid.write(f'{flatr[ind][1]}\t{flatr[ind][0]}')
	testx.close()
	testid.close()
	print(len(biglist)-1)
	print(len(flatr[1][2]))
	testy.close()
		
	##now training
	for ind in range(len(trainr)):
		trainr[ind] = onehot(trainr[ind], biglist)
		featstrr = '\t'.join(str(x) for x in trainr[ind][2])
		trainx.write(f'{featstrr}\n')
		trainy.write(f'{trainr[ind][1]}\n')
		trainid.write(f'{trainr[ind][1]}\t{trainr[ind][0]}')
	trainx.close()
	trainy.close()
	trainid.close()
	
	print('Done making inputs!')
	
def main(inmgfs, mgfnams, binw, outp):
	mgfnames = parse_args(inmgfs, mgfnams)
	mgfd = parse_mgf(mgfnames)
	make_output(mgfd, binw, outp)
	print('All done!')
	
if __name__ == '__main__':
	print('poo!')
	print(sys.argv)
	if len(sys.argv) != 5:
		sys.exit('ARGS: 1) file path(s) to input mgf(s). If multiple,\n'\
				 'must be separated by commas. 2) Metabolite class names of \n'\
				 'each mgf 3) binwidth 4) output file path')
	
	pths = sys.argv[1]
	names = sys.argv[2]
	binw = float(sys.argv[3])
	outp = sys.argv[4]
	
	main(pths, names, binw, outp)