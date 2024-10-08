import random


def Count(Motifs):
    count = {}
    for nucleotide in "ATCG":
        count[nucleotide] = [0] * len(Motifs[0])

    for motif in Motifs:
        for index, nucleotide in enumerate(motif):
            count[nucleotide][index] += 1
    return count


print(Count(["ACATTG", 'CATGCA', 'GAGTCA', 'TACAGT']))


def Profile(Motifs):
    t = len(Motifs)
    profile = {}
    counts = Count(Motifs)
    for nucleotide in counts:
        for _ in counts[nucleotide]:
            profile[nucleotide] = [x / t for x in counts[nucleotide]]
    return profile


def Consensus(Motifs):
    consensus = ""
    k = len(Motifs[0])
    count = Count(Motifs)
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    for motif in Motifs:
        for m, c in zip(motif, consensus):
            if m != c:
                score += 1
    return score


# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    probability = 1
    for i, nucleotide in enumerate(Text):
        probability *= Profile[nucleotide][i]
    return probability


def ProfileMostProbableKmer(text, k, profile):
    kmers_probabilities = {}
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        probability = Pr(pattern, profile)
        kmers_probabilities[pattern] = probability
    return max(kmers_probabilities, key=kmers_probabilities.get)


def Motifs(Profile, Dna):
    motifs = []
    k = len(Profile["A"])
    for text in Dna:
        motifs.append(ProfileMostProbableKmer(text, k, Profile))
    return motifs


def GreedyMotifSearch(Dna, k, t):
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(best_motifs):
            best_motifs = Motifs
    return best_motifs


# Week four

def CountWithPseudocounts(Motifs):
    count = {}
    for nucleotide in "ATCG":
        count[nucleotide] = [1] * len(Motifs[0])

    for motif in Motifs:
        for index, nucleotide in enumerate(motif):
            count[nucleotide][index] += 1
    return count


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    profile = {}
    counts = CountWithPseudocounts(Motifs)
    for nucleotide in counts:
        for _ in counts[nucleotide]:
            profile[nucleotide] = [x / (t + 4) for x in counts[nucleotide]]
    return profile


def ConsensusWithPseudocounts(Motifs):
    consensus = ""
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def ScoreWithPseudocounts(Motifs):
    score = 0
    consensus = ConsensusWithPseudocounts(Motifs)
    for motif in Motifs:
        for m, c in zip(motif, consensus):
            if m != c:
                score += 1
    return score


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if ScoreWithPseudocounts(Motifs) < ScoreWithPseudocounts(best_motifs):
            best_motifs = Motifs
    return best_motifs


def RandomMotifs(Dna, k, t):
    import random
    motifs = []
    for text in Dna:
        i = random.randint(0, len(text) - k)
        motifs.append(text[i:i + k])
    return motifs


def RandomizedMotifSearch(Dna, k, t):
    # M = RandomMotifs(Dna, k, t)
    M = ["GTC", "CCC", "ATA", "GCT"]
    BestMotifs = M
    print("iteration")
    Profile = ProfileWithPseudocounts(M)
    M = Motifs(Profile, Dna)
    print(M)
    if ScoreWithPseudocounts(M) < ScoreWithPseudocounts(BestMotifs):
        BestMotifs = M
    else:
        return BestMotifs


dna = ["ATGAGGTC", "GCCCTAGA", "AAATAGAT", "TTGTGCTA"]
print(RandomizedMotifSearch(dna, 3, 1))


def Normalize(Probabilities):
    total = sum(Probabilities.values())
    for key in Probabilities:
        Probabilities[key] /= total
    return Probabilities


prs = {"A": 0.15, "B": 0.6, "C": 0.225, "D": 0.225, "E": 0.3}
print(Normalize(prs).values())
print(sum([0.1, 0.4, 0.15, 0.15, 0.2]))


def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    for kmer, probability in Probabilities.items():
        p -= probability
        if p <= 0:
            return kmer


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


def GibbsSampler(Dna, k, t, N):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    for j in range(N):
        i = random.randint(0, t - 1)
        Profile = ProfileWithPseudocounts(M[:i] + M[i + 1:])
        Motifs = M[:i] + [ProfileGeneratedString(Dna[i], Profile, k)] + M[i + 1:]
        if ScoreWithPseudocounts(Motifs) < ScoreWithPseudocounts(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# print(GibbsSampler(["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
#                     "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
#                     "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"], 8, 5, 100))

# text = "AAACCCAAACCC"
# profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
# print(ProfileGeneratedString(text, profile, 2))
# probs = {"A": 0.1, "C": 0.1, "G": 0.3, "T": 0.5}
# result = {"A": 0, "C": 0, "G": 0, "T": 0}
# for i in range(10000):
#     result[WeightedDie(probs)] += 1
# print(result)

# str = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
# # print(GreedyMotifSearch(str, 3, 5))
# # print(GreedyMotifSearchWithPseudocounts(str, 3, 5))
# print(GreedyMotifSearchWithPseudocounts(str, 3, 5))
# print(RandomizedMotifSearch(str, 3, 5))
#
# pr = {"A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9], "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1], "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
#       "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}
