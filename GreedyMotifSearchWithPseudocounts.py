import numpy

def greedy_motif_search(dna, k, t):
    best_motifs = []

    for i in range(t):
        best_motifs.append(dna[i][0:k])

    for i in range(len(dna[0]) - k + 1):
        motifs = []
        motifs.append(dna[0][i:i+k])

        for j in range(1, t):
            profile = build_profile(motifs, k)

            motifs.append(profile_most_probable(profile, dna[j], k))

        if score(motifs) < score(best_motifs):
            best_motifs = motifs.copy()

    return best_motifs


def build_profile(motifs, k):
    profile = numpy.zeros((4, k))

    for r in profile:
        for c in range(len(r)):
            r[c] = 1

    for col in range(k):
        for r in range(len(motifs)):
            if motifs[r][col] == 'A':
                profile[0][col] += 1
            if motifs[r][col] == 'C':
                profile[1][col] += 1
            if motifs[r][col] == 'G':
                profile[2][col] += 1
            if motifs[r][col] == 'T':
                profile[3][col] += 1

    return profile / (len(motifs) + 4)



def profile_most_probable(profile, dna, k):
    max_prob = -1
    probable = []

    # loop through each kmer in dna
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prob = 1
        j = 0

        for n in kmer:
            if n == 'A':
                prob *= profile[0][j]
            if n == 'C':
                prob *= profile[1][j]
            if n == 'G':
                prob *= profile[2][j]
            if n == 'T':
                prob *= profile[3][j]
            j += 1
        if prob > max_prob:
            max_prob = prob
            probable.append([kmer, max_prob])
            most_probable = kmer

    return most_probable



def read_dna(dna, k, t):
    dnalist = []
    n = len(dna)/t
    print(n)
    print(len(dna))
    for i in range(0, len(dna)):
        if dna[i][0] == '>':
            print(dna[i][0])
            continue
        dnalist.append(dna[i])
    return dnalist

def score(motif):
    score = 0
    consensus = max(set(motif), key=motif.count)
    for idx in motif:
        score += hamming_distance(idx, consensus)
    return score


def hamming_distance(str1, str2):
    hd = 0
    for idx in range(len(str1)):
        if str1[idx] != str2[idx]:
            hd += 1
    return hd


if __name__ == '__main__':
    # Set values for k and t
    # The file name should be a text file containing a DNA sequence
    file = 'dm01r.fasta'
    with open(file) as f:
        dna = f.read().splitlines()

    x = read_dna(dna, 5, 4)
    for i in greedy_motif_search(x, 5, 4):
        print(i)

