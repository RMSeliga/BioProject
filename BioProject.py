
#Reader Function: Reads the sequence and the associated probabilities.
#Could work one of two ways
#1: Creates two arrays of the same length as the sequence, each index has the corresponding nucleotide and probability
#2: Creates one array of tuples with length equal to the sequence, each tuple contains the nucleotide + probability (could probably also do with sets for better runtime)
import numpy as np
from string import *
sequence = []
sequenceN = ""
def Reader(file1, file2, s):
	s=""
	f = open(file1)
	for line in f:
		for nucleotide in line:
			sequence.append([nucleotide,None])
			s+=nucleotide

	# print len(s)
	# print s
	
	i=0;
	f2 = open(file2)
	for line in f2:
		for p in line.split():
			sequence[i][1]=float(p)
			i+=1
	return s

sequenceN = Reader("sequence.txt", "probabilities.txt", sequenceN)
print len(sequenceN)

def Diagonal(n1,n2,pt):
    if(n1 == n2):
        return pt['MATCH']
    else:
        return pt['MISMATCH']

def Pointers(di,ho,ve):

    pointer = max(di,ho,ve) #based on python default maximum(return the first element).

    if(di == pointer):
        return 'd'
    elif(ho == pointer):
        return 'h'
    else:
         return 'v' 

def NeedlemanWunsch(query, dictionary, match, mismatch, gap, start):
	
	penalty = {'MATCH': match, 'MISMATCH': mismatch, 'GAP': gap}
	n = len(query)+1
	m = len(dictionary)+1
	D = np.zeros((m,n),dtype = str)#[["" for x in range(len(T))]for y in range(len(S))]
	M = np.zeros((m,n),dtype = float)#[[0 for x in range(len(T))]for y in range(len(S))]
	#print len(M)
	#print len(M[0])
	# m[0][0] = 0
	query = "-" + query
	dictionary = "-" + dictionary
	for i in range(m):
		M[i][0] = gap*i * sequence[i][1]
		D[i][0] = 'v'
	for j in range(n):
		M[0][j] = gap*j * sequence[j][1]
		D[0][j] = 'h'
	M[0][0] = 0
	for i in range(1,m):
		for j in range(1,n):
			diag = M[i-1][j-1] + Diagonal(query[j-1],dictionary[i-1], penalty)
			horizontal = M[i][j-1] + gap
			vertical = M[i-1][j] + gap
			M[i][j] = max(diag,horizontal,vertical) * sequence[i][1]
			D[i][j] = Pointers(diag, horizontal, vertical)
	# print np.matrix(M)
	# print np.matrix(D)
	finalSeq=""
	finalDictSeq=""
	x = 0
	y = 0
	
	#print m

	# for p in range(m)[::-1]:
	# 	for q in range(n)[::-1]:
	# 		if M[p][q] > 0 and M[p][q] > M[p-1][q-1] and M[p][q] > M[p-1][q] and M[p][q] > M[p][q-1]:
	# 			# print M[p][q]
	# 			x = p
	# 			y = q
	# 			print "score: " + str(M[p][q])
	# 			break
	# 	if (x!=0):
	# 		break
				
	# print p 
	# print q
	


	# x = m-2
	# y = n -2
	x=m-2
	y=n-2
	while True:
		#backtrack 
		#print str(x) + " " + str(y)
		if D[x][y] == 'h':
			finalSeq += query[y]
			finalDictSeq += "-"#dictionary[x]
			y -= 1
		elif D[x][y] == 'v':
			finalSeq += "-"#query[y]
			finalDictSeq += dictionary[x]#"-"
			x -= 1
		elif D[x][y] == 'd':
			finalSeq += query[y]
			finalDictSeq += dictionary[x]
			x-=1
			y-=1
		if (x==0 or y==0):
			break
		#if x == 0 or y == 0: break
	# print query
	# print finalSeq[::-1]
	# print finalDictSeq[::-1]
	# # print len(finalSeq[::-1])
	# # print len(dictionary)
	# print "Score: " + str(M[m-1][n-1])
	return [M[m-1][n-1], finalSeq[::-1], finalDictSeq[::-1], start]
			#check if cell[x-+1][y+1], [x-1][y], [x][y-1]


	# for x in range(m):

	# 	# for y in range(n):
	# 		# if D[x][y]=='h':
	# 		# 	finalseq += "-"
	# 		# 	finalDictSeq += sequenceN[i:(i+len(target)*2)][x]
	# 		# else if D[x][y] == 'v':
	# 		# 	finalseq += target[x]
	# 		# 	finalDictSeq += "-"
	# 		# else if D[x][y] == 'd':
	# 		# 	finalseq += target[x]
	# 		# 	finalDictSeq += sequenceN[i:(i+len(target)*2)][x]
	# 		# if M[x][y+1] > M[x][y]:
	# 	if M[x][x]<M[x-1][x-1]:
	# 		print target[len(target)-x:len(target)]








#Finder/Matcher
#This should use the blast algorithm with a small addition
#Basically it does the standard best fit search to find the most likely origin point for the query sequence, only now the confidence value also incorporates the probabilites
#I.E. at each point, when comparing nucleotides from the query to the dictionary if the probability of the dictionary nucleotide being correct is low, then the confidence value from that index should also be low
#E.G a matching that has lots of matches but very low probabilities for those matches will be less likely to be selected than a matching with fewer matches, but very high probabilites at those matches
#This should also be true for the opposite case where there is a mismatch
target = 'AGACGGTTTCTCTCTTTGCAGCATCACCCAGGCTGGGAGT'
def Matcher(target, wordSize):
	counter = 0
	matches = []
	chunks = [target[i:i+wordSize] for i in range(0, len(target))]
	chunks = chunks[0:-(wordSize-1)]
	#for seqChunk in [sequenceN[index:index+wordSize] for i in range(0, len(sequenceN))]:
	for i in range(0, len(sequenceN)):
		seqChunk = sequenceN[i:i+wordSize]
		for c in chunks:
			if c == seqChunk:
				start = i-chunks.index(c)
				end = start+len(target)
				matches.append(NeedlemanWunsch(target,sequenceN[start:end],1,-1,-1, start))
				
				
	maximum = 0
	mSeq = ""
	mDict = ""
	location = 0
	for j in matches:
		if j[0]>maximum:
			
			maximum=j[0]
			mSeq = j[1]
			mDict = j[2]
			location = j[3]
	print maximum
	print mSeq
	print mDict
	print location
	

	#print counter
	


Matcher(target, 11)

#M * N from comparing each chunk of query to each chunk of dictionary
#N^2 from needleman wunsch
#O(M*N^3)