
import sys
import os
import inspect
import re
import GibbsSampler

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.insert(0, path)



def DemoMotifFinder(numSeeds, k, N):
	filename = inspect.getframeinfo(inspect.currentframe()).filename
	path = os.path.dirname(os.path.abspath(filename))
	mm9_file = open(path + "/dm01r.txt", 'r')
	
	Dna = []
	for row in mm9_file:
		if row[0] != '>':
			row = row.strip()
			Dna.append(row)

	BestMotif = GibbsSampler.MultipleSeedsGibbsSampling(Dna, numSeeds, k, N)
	mm9_file.close()
	return BestMotif


numSeeds = 20
k = 11
N = 1000
BestMotif = DemoMotifFinder(numSeeds, k, N)
best_score = GibbsSampler.score(BestMotif)
print("Gibbs Sampler Motifs")
print("best_score:   ", best_score)

mm9SolLoc = path + "/dm01rSol.txt"
mm9_solutions = open(mm9SolLoc, 'r')
solutions_motif = []
for line in mm9_solutions:
	if line[0] != '>':
		line = re.sub('[^A-Z]', '', line)
		solutions_motif.append(line)
print("Real Motifs")
print("Real Score:   ", GibbsSampler.score(solutions_motif))
print()
print("Algorithms Best Pick", "     Match/Wrong      ", "Real Motif")
for i in range(len(solutions_motif)):
	if solutions_motif[i] == BestMotif[i]:
		print(BestMotif[i], "  Match     ", solutions_motif[i])
	else:
		print(BestMotif[i], "  Wrong     ", solutions_motif[i])

	
