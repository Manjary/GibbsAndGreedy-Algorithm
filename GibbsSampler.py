

import random
import operator

def MultipleSeedsGibbsSampling(dna, numSeeds, k, N):
	results = GibbsSampler(dna, k, N)
	best_score = score(results)
	best_motifs = list(results)
	for i in range(1, numSeeds):
		results = GibbsSampler(dna, k, N)
		if score(results) < best_score:
			best_score = score(results)
			best_motifs = list(results)
	return best_motifs
	
def GibbsSampler(dna, k, N):
	t = len(dna)
	motifs = []
	for strand in dna:
		i = random.randrange(len(strand)-k+1)
		substr = strand[i:i+k]
		motifs.append(substr)
		
	best_motifs = list(motifs)
	best_motifs_score = score(best_motifs)
	for j in range(1,N):
		i = random.randrange(t)
		subset_motifs = motifs[0:i]+motifs[i+1:t]
		replacement_motif = SingleReplacementMotif(subset_motifs, dna[i])
		motifs[i] = replacement_motif
		
		if score(motifs) < best_motifs_score:
			best_motifs = list(motifs)
			best_motifs_score = score(best_motifs)
	return best_motifs

	
def SingleReplacementMotif(motifs, dna_i):


	k = len(motifs[0])
	profile = BuildProfile(motifs)
	
	kmer_densities = [0 for x in range(len(dna_i)-k+1)]
	for i in range(len(dna_i)-k+1):
		prob = 1
		for j in range(k):
			if dna_i[i+j] == 'A':
				prob *= profile[j][0]
			elif dna_i[i+j] == 'C':
				prob *= profile[j][1]
			elif dna_i[i+j] == 'G':
				prob *= profile[j][2]
			elif dna_i[i+j] == 'T':
				prob *= profile[j][3]
		kmer_densities[i] = prob

	# normalize probabilities
	normalization_tot = sum(kmer_densities)
	for i in range(len(dna_i)-k+1):
		kmer_densities[i] = kmer_densities[i]/normalization_tot

	# construct prefix sum for lookup		
	kmer_densities = list(accumulate(kmer_densities))
	
	# randomly select a k-mer
	rand_val = random.random()
	for i in range(len(dna_i)-k+1):
		if rand_val < kmer_densities[i]:
			break
	replacement_kmer = dna_i[i:i+k]

	return replacement_kmer

def BuildProfile(motif):
	k = len(motif[0])
	profile = [[0 for y in range(4)] for x in range(k)]
	for count in range(k):
		A, C, G, T = 0, 0, 0, 0
		A += 1
		C += 1
		G += 1
		T += 1
		for string in motif:
			if string[count] == 'A':
				A += 1
			elif string[count] == 'C':
				C += 1
			elif string[count] == 'G':
				G += 1
			elif string[count] == 'T':
				T += 1
		# Insert frequencies if base A
		profile[count][0] = float(A) / (A + C + G + T)
		# Insert frequencies if base C
		profile[count][1] = float(C) / (A + C + G + T)
		# Insert frequencies if base G
		profile[count][2] = float(G) / (A + C + G + T)
		# Insert frequencies if base T
		profile[count][3] = float(T) / (A + C + G + T)
	return profile

def BuildMotifs(profile, dna, k):

	motif = []
	for string in dna:
		best_sub_str = ''
		for i in range(len(string)+1-k):
			substr = string[i:i+k]
			prob = 1
			best_prob = -1
			for j in range(k):
				if substr[j] == 'A':
					prob *= profile[j][0]
				elif substr[j] == 'C':
					prob *= profile[j][1]
				elif substr[j] == 'G':
					prob *= profile[j][2]
				elif substr[j] == 'T':
					prob *= profile[j][3]
			if prob > best_prob:
				best_prob = prob
				best_sub_str = substr
		motif.append(best_sub_str)
	return motif

def score(motifs):
	k = len(motifs[0])
	pattern = []
	for i in range(k):
		A=0
		C=0
		G=0
		T=0
		for string in motifs:
			if string[i]=='A':
				A+=1
			elif string[i]=='C':
				C+=1
			elif string[i]=='G':
				G+=1
			elif string[i]=='T':
				T+=1
		
		if A >= C and A >= G and A >= T:
			pattern.append('A')
		elif C >= G and C >= T:
			pattern.append('C')
		elif G >= T:
			pattern.append('G')
		else:
			pattern.append('T')

	pattern = "".join(pattern)
			
	score = 0
	for string in motifs:
		score += HammingDistance(string, pattern)
	return score


def d(kmer, dna):
	k = len(kmer)
	motif = []
	totDist = 0
	for phrase in dna:
		localDist = len(phrase)+len(kmer)
		word = ""
		for i in range(len(phrase)-k+1):
			subPattern = phrase[i:i+k]
			if localDist > HammingDistance(kmer, subPattern):
				localDist = HammingDistance(kmer, subPattern)
				word = subPattern
		motif.append(word)
		totDist += localDist
	return totDist


def HammingDistance(str1, str2):
	diffs = 0
	for ch1, ch2 in zip(str1, str2):
		if ch1 != ch2:
			diffs += 1
		return diffs
		
		
def accumulate(iterable, func=operator.add):
	it = iter(iterable)
	try:
		total = next(it)
	except StopIteration:
		return
	yield total
	for element in it:
		total = func(total, element)
		yield total
