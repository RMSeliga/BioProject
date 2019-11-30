#Reader Function: Reads the sequence and the associated probabilities.
#Could work one of two ways
#1: Creates two arrays of the same length as the sequence, each index has the corresponding nucleotide and probability
#2: Creates one array of tuples with length equal to the sequence, each tuple contains the nucleotide + probability (could probably also do with sets for better runtime)
sequence = []
def Reader(file1, file2):
	f = open(file1)
	for line in f:
		for nucleotide in line:
			sequence.append([nucleotide,None])
	print len(sequence)
	
	i=0;
	f2 = open(file2)
	for line in f2:
		for p in line.split():
			sequence[i][1]=p
			i+=1

Reader("sequence.txt", "probabilities.txt")
#Finder/Matcher
#This should use the blast algorithm with a small addition
#Basically it does the standard best fit search to find the most likely origin point for the query sequence, only now the confidence value also incorporates the probabilites
#I.E. at each point, when comparing nucleotides from the query to the dictionary if the probability of the dictionary nucleotide being correct is low, then the confidence value from that index should also be low
#E.G a matching that has lots of matches but very low probabilities for those matches will be less likely to be selected than a matching with fewer matches, but very high probabilites at those matches
#This should also be true for the opposite case where there is a mismatch

 