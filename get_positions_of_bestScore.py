### ----------------------------
### TFBS Interdistances program
### ----------------------------

'''
This program allows to calculate interdistances between transcription factor binding sites.
You need a matrix with frequency values and fasta sequences.
This program was written by Adrien Bessy and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import numpy as np
from Bio import SeqIO
import time
import sys
#from tqdm import *  
from operator import truediv
import argparse
from optparse import OptionParser

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type = str, default= "ARF2")
parser.add_argument("--scores", "-scores", type = str, default= "bestScores")
parser.add_argument("--pseudoCount", "-pc",type = float, default = 0.001)
args = parser.parse_args()

scores = args.scores
factorTranscription = args.factor
pseudoCount = args.pseudoCount

# python get_positions_of_bestScore.py -fac "ARF2" -pc 0.001 -scores "bestScores"

if factorTranscription == "ARF2" :
	FastaFile = "../sequences/ARF2_bound_sequences.fas" 
	MatrixFile = "../matrices/ARF2_matrix1_MEME_1500FirstSeq.txt" 
	matrixType = "freq" 
	dependencyFile = "" 
	
if factorTranscription == "ARF5" :
	FastaFile = "../sequences/ARF5_bound_sequences.fas" 
	MatrixFile = "../matrices/ARF5_allSeq_3prime_freq_pasteTo'OMalley.txt" 
	matrixType = "freq" 
	dependencyFile = "" 
	
if factorTranscription == "TGTC_on_ARF2seq" :
	FastaFile = "../sequences/ARF2_bound_sequences.fas" 
	MatrixFile = "../matrices/TGTC.txt" 
	matrixType = "freq" 
	dependencyFile = "" 
	
if factorTranscription == "TGTC_on_ARF5seq" :
	FastaFile = "../sequences/ARF5_bound_sequences.fas" 
	MatrixFile = "../matrices/TGTC.txt" 
	matrixType = "freq" 
	dependencyFile = "" 

''' separation between numbers can be spaces, tabulation, comas...
                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
                                                                        '''

#FastaFile = "bound_sequences.fas"
#FastaFileN = "MP_neg_1.fas"
#FastaFile = "oneSequence.fas"

def all(MatrixFile,FastaFile) :
	# These 3 lines allows to retrieve the matrix from the file
	F = open(MatrixFile,"r")
	matrix = F.read().replace("\r","\n") + "\n"
	F.close()

	# These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list
	import re
	num = re.compile(r"([+-]?\d+[.,]\d+)")
	Mdata = num.findall(matrix)

	## These lines allows to transform the frequency values into scores values
	matF = []
	lenMotif=0
	for i in range(0,len(Mdata)):
		if i%4==0:
			lenMotif=lenMotif+1
			fmax = float(max(Mdata[i],Mdata[i+1],Mdata[i+2],Mdata[i+3])) + pseudoCount
			for j in range (0,4):
				matF.append(np.log(float(float(Mdata[i+j]) + pseudoCount) /fmax))

	# This line allows to retrieve all the sequences from the fasta file
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))

	print "  There are %s sequence(s) to analyze"%(len(sequences))

	# The following line allows to produce the reversed matrix
	'''if we take the example given before : A T G C
				Position 1:      0.4444  0.155  0.654   0.645
				Position 2:      0.1645  0.1565 0.21614 0.16456
	Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
	So we can calculate with this reverse matrix, the score of the complementary strand.
	'''
	matRev = list(reversed(matF))
	
	positions =[]
	bestScoreBySeq = []
	# We look at all the fasta sequences:
	#for s in tqdm(sequences):
	if scores == "bestScores" :
		for s in sequences:
			# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
			GoodScorePositionsStrandPos = []

			# This line allows to retrieve the DNA sequence
			seq = sequences[s].seq
			score_seq = []
			# We look at each sub-sequences of the whole sequence. Each sub-sequence has the same length that the matrix length.
			for c in range(len(seq)-lenMotif+1):
				#print(c)
				#print("lenMotif+1 : ",lenMotif+1)
				strandPos = seq[c:c+lenMotif].upper()
				test = 0
				for nu in strandPos :
					if nu not in ["A","C","G","T"]:
						test=1
				if test == 1:
					score = "NA"
				else :
					n=0
					#These lines allows to calculate a score for one sub-sequence
					scoreStrandPos=0
					scoreStrandNeg=0
					while n < lenMotif:
						if strandPos[n] == 'A':
							scoreStrandPos = scoreStrandPos + matF[n*4]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4]
						elif strandPos[n] == 'C':
							scoreStrandPos = scoreStrandPos + matF[n*4+1]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+1]
						elif strandPos[n] == 'G':
							scoreStrandPos = scoreStrandPos + matF[n*4+2]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+2]
						elif strandPos[n] == 'T':
							scoreStrandPos = scoreStrandPos + matF[n*4+3]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+3]		
						n += 1
					if factorTranscription == "ARF2":
						score_seq.append([scoreStrandPos,c+1])
						score_seq.append([scoreStrandNeg,c+9])
					if factorTranscription == "ARF5":
						score_seq.append([scoreStrandPos,c+3])
						score_seq.append([scoreStrandNeg,c+7])
					#if scoreStrandPos == 0:
						#positions.append(c+1)
					#if scoreStrandNeg == 0:
						#positions.append(c+1)       
			bestScoreBySeq.append(max(score_seq, key=lambda x: x[0]))
			#print("bestScoreBySeq : ",bestScoreBySeq)
		positions = [item[1] for item in bestScoreBySeq]
	if scores == "all_scores" :
		for s in sequences:
			# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
			GoodScorePositionsStrandPos = []

			# This line allows to retrieve the DNA sequence
			seq = sequences[s].seq
			score_seq = []
			# We look at each sub-sequences of the whole sequence. Each sub-sequence has the same length that the matrix length.
			for c in range(len(seq)-lenMotif+1):
				#print(c)
				#print("lenMotif+1 : ",lenMotif+1)
				strandPos = seq[c:c+lenMotif].upper()
				test = 0
				for nu in strandPos :
					if nu not in ["A","C","G","T"]:
						test=1
				if test == 1:
					score = "NA"
				else :
					n=0
					#These lines allows to calculate a score for one sub-sequence
					scoreStrandPos=0
					scoreStrandNeg=0
					while n < lenMotif:
						if strandPos[n] == 'A':
							scoreStrandPos = scoreStrandPos + matF[n*4]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4]
						elif strandPos[n] == 'C':
							scoreStrandPos = scoreStrandPos + matF[n*4+1]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+1]
						elif strandPos[n] == 'G':
							scoreStrandPos = scoreStrandPos + matF[n*4+2]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+2]
						elif strandPos[n] == 'T':
							scoreStrandPos = scoreStrandPos + matF[n*4+3]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+3]		
						n += 1
					if scoreStrandPos == 0:
						positions.append(c)
					#if scoreStrandNeg == 0:
						#positions.append(c+3)       
			#bestScoreBySeq.append(max(score_seq, key=lambda x: x[0]))
			#print("bestScoreBySeq : ",bestScoreBySeq)
		#positions = [item[1] for item in bestScoreBySeq]
	return(positions)

positions = all(MatrixFile,FastaFile)
#positionsN = all(MatrixFile,FastaFileN)
	

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import matplotlib.patches as mpatches
from matplotlib import pylab

fig = plt.figure(1,figsize= (18,10))
ax = plt.subplot(111)

from scipy import stats

labels1, values1 = zip(*Counter(positions).items() )
#labels2, values2 = zip(*Counter(positionsN).items() )

#values3 = [a_i - b_i for a_i, b_i in zip(values1, values2)]
fig.suptitle(sys.argv, fontsize = 11, fontweight='bold')

ax.bar(labels1, values1, color='cornflowerblue')
ax.axis([1, 201, 0, max(values1)])
ax.set_xlabel("Positions of the first T of the motif TGTC that has the best score on each sequence", size = 16)
ax.set_ylabel("Occurences", size = 16)

plt.show()


