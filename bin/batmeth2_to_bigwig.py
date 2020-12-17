#!/usr/bin/env python
import sys, math, multiprocessing, subprocess, os, gzip


# Usage: python3 allc_to_bigwig_pe.py [-keep] [-sort] [-L=labels] [-p=num_proc] [-c=base_mod]  <chrm_sizes>  <allC_file> [allC_file]*

# NOTE: allc file contains the methylation information for all chromosomes

# Steps:
# 1. allC to bedGraph
# 2. sort bedGraph if necessary
# 3. bedGraph to BigWig
# 4. remove temporary files

NUMPROC=1

def processInputs( allCFileAr, chrmFileStr, keepTmp, labelsAr, outID, baseMod, numProc, isSort ,printStrand):
	print( 'Keep temp files: {:s}\nSort bedGraph: {:s}\nContext/Modification: {:s}'.format( str( keepTmp), str (isSort), ('CG, CHG, CHH' if baseMod == None else baseMod)))
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(allCFileAr[i], chrmFileStr, labelsAr[i], outID, baseMod,  keepTmp, isSort,printStrand) ) for i in range(len(allCFileAr)) ]
	suc = [ p.get() for p in results ]

	print( 'Done' )


def processFile( allCFileStr, chrmFileStr, label, outID, baseMod, keepTmp, isSort, printStrand):
	tmpAllc = allCFileStr.replace( '.tsv','' ).replace( 'allc_','' ).replace('.gz','').replace('.gzip', '').replace('methratio.txt', 'mr')
	if outID == None and label == None:
		outID = tmpAllc
	elif outID == None:
		outID = label
	elif label == None:
		outID += '_' + tmpAllc
	else:
		outID += '_' + label
	
	outID += '_bt2'
	if printStrand:
		outID += '_strand'

	print( 'Reading allC file {:s}'.format( allCFileStr ) )
	bedGraphStr =  outID # + '.bedGraph'

	# default base modifications
	if baseMod == None:
		bedGraphAr = [bedGraphStr + '.' + x for x in ['cg.bedGraph','chg.bedGraph','chh.bedGraph'] ]
	# custom modifications
	else:
		bedGraphAr = [bedGraphStr + '.' + x for x in [baseMod+'.bedGraph'] ]
	#bedGraphStr = bedGraphStr + '.bedGraph'
	OutCfileStr = bedGraphStr + '.c.bedGraph'

	readAllC( allCFileStr, bedGraphAr, OutCfileStr ,printStrand)

	if isSort:
		print( 'Sorting bedGraph files' )
		for b in bedGraphAr:
			sortBedFile( b )
		sortBedFile( OutCfileStr )

	print( 'Converting {:s} files to BigWig'.format(bedGraphStr ) )
	# bedGraph to bigWig
	for b in bedGraphAr:
		processBedGraph( b, chrmFileStr )
	processBedGraph(OutCfileStr, chrmFileStr)
	# remove temporary
	if keepTmp == False:
		print( 'Removing temporary files...' )
		for b in bedGraphAr:
			os.remove( b )
		os.remove(OutCfileStr)
	print( 'BigWig finished for {:s}.bw.*'.format( outID ) )

def readAllC( allCFileStr, outFileAr, OutCfileStr, printStrand):

	if allCFileStr.endswith('.gz') or allCFileStr.endswith('.gzip'):
		allCFile = gzip.open(allCFileStr, 'rt')
	else:
		allCFile = open( allCFileStr, 'r' )
	outFileAr = [open( outFileStr, 'w' ) for outFileStr in outFileAr]
	OutCfile = open( OutCfileStr, 'w' )

	mTypes = [ 'CG', 'CHG', 'CHH' ]

	isContext = (len(outFileAr) > 1)

	for line in allCFile:
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
                #chrM    1       -       CHH     2       1795    0.001114
		if line.startswith('#') or len(lineAr) < 7:  #or lineAr[6].isdigit() == False:
			continue
		elif int(lineAr[4]) >= 0:
			chrm = lineAr[0]
			pos = int( lineAr[1] ) - 1
			mInd = 0
			if isContext:
				methType = decodeMethType( lineAr[3] )
				try:
					mInd = mTypes.index( methType )
				except ValueError:
					continue
			value = float( lineAr[4] ) / float( lineAr[5] )
			#if value == 0:  #skip 0
			#	continue
			isMeth = 0
			if methType == 'CG':
				if value >= 0.7:
					isMeth = 1
			elif methType == 'CHG':
				if value >= 0.5:
					isMeth = 1
			elif methType == 'CHH':
				if value >= 0.3:
					isMeth = 1
			#isMeth = int(lineAr[6])
			# adjust for negative strand
			if printStrand and lineAr[2] == '-':
				value = value * -1

			# (0) chrm (1) start (2) end (3) value
			outStr = '{:s}\t{:d}\t{:d}\t{:.6f}{:d}\n'.format( chrm, pos, pos+1, value,isMeth )
			outFile = outFileAr[mInd]
			outFile.write( outStr )
			OutCfile.write( outStr )
		# end if
	# end for
	allCFile.close()
	[ outFile.close() for outFile in outFileAr ]
	OutCfile.close()

def decodeMethType( mStr ):

	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	elif mStr == 'CNN':
		return 'CNN'
	else:
		return 'CHH'

def sortBedFile( bedFileStr ):
	command = 'bedSort {:s} {:s}'.format( bedFileStr, bedFileStr )
	subprocess.call( command, shell=True )

def processBedGraph( bedGraphStr, chrmFileStr ):

	bigWigStr = bedGraphStr.replace( '.bedGraph', '.bw' )
	#print( bigWigStr )
	# bedGraphToBigWig in.bedGraph chrom.sizes out.bw
	command = 'bedGraphToBigWig {:s} {:s} {:s}'.format( bedGraphStr, chrmFileStr, bigWigStr )
	subprocess.call( command, shell=True)


def parseInputs( argv ):
	# Usage: python3 allc_to_bigwig_pe.py [-keep] [-sort] [-L=labels] [-p=num_proc] [-c=base_mod] <chrm_sizes>  <allC_file> [allC_file]*

	keepTmp = False
	labelsAr = []
	numProc = NUMPROC
	isSort = False
	outID = None
	baseMod = None
	startInd = 0
	printStrand = False

	for i in range( min(5, len(argv)-2) ):
		if argv[i] == '-keep':
			keepTmp = True
			startInd += 1
		elif argv[i] == '-sort':
			isSort = True
			startInd += 1
		elif argv[i].startswith( '-L=' ):
			labelsAr = argv[i][3:].split( ',' )
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			baseMod = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-strand' ):
			printStrand = True
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be an integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	# set use all if not already set
	chrmFileStr = argv[startInd]
	allCFileAr = []
	for j in range(startInd+1, len(argv) ):
		allCFileAr += [ argv[j] ]

	if len(labelsAr) == 0:
		labelsAr = [None] * len(allCFileAr)
	elif len(labelsAr) != len(allCFileAr):
		print( "ERROR: number of labels doesn't match number of allC files" )
		exit()

	processInputs( allCFileAr, chrmFileStr, keepTmp, labelsAr, outID, baseMod, numProc, isSort ,printStrand)

def printHelp():
	print ("Usage: python3 allc_to_bigwig_pe.py [-keep] [-sort] [-L=labels] [-o=out_id] [-p=num_proc] <chrm_sizes>  <allC_file> [allC_file]*")
	print( 'Converts allC files to context-specific BigWig files' )
	print( 'Note: bedGraphToBigWig and bedSort programs must be in the path' )
	print( 'Required:' )
	print( 'chrm_file\ttab-delimited file with chromosome names and lengths,\n\t\ti.e. fasta index file' )
	print( 'allc_file\tallc file with all chrms and contexts\n\t\tcan be gzip compressed' )
	print( 'Optional:' )
	print( '-keep\t\tkeep intermediate files' )
	print( '-sort\t\tcalls bedSort; add this option if bigwig conversion fails' )
	print( '-L=labels\tcomma-separated list of labels to use for the allC files;\n\t\tdefaults to using information from the allc file name' )
	print( '-c=base_mod\tbase modification for file extension; single modification for\n\t\tCall input files, i.e. 5hmC; overrides default G,CHG,CHH output' )
	print( '-o=out_id\toptional identifier to be added to the output file names' )
	print( '-p=num_proc\tnumber of processors to use [default 1]' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
