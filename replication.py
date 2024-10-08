def PatternCount(Text, Pattern):
    count = 0
    p_size = len(Pattern)
    for i in range(len(Text) - p_size + 1):
        if Text[i:i + p_size] == Pattern:
            count += 1
    return count


def FrequencyMap(Text, k):
    freq = {}
    text_size = len(Text)
    for i in range(text_size - k + 1):
        pattern = Text[i:i + k]
        freq[pattern] = 0
    for item in freq:
        freq[item] = PatternCount(Text, item)
    return freq


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def FrequentWordsWithCount(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append({key: freq[key]})
    return words


def PatternMatching(Pattern, Genome):
    positions = []
    p_size = len(Pattern)
    for i in range(len(Genome) - p_size + 1):
        if Genome[i:i + p_size] == Pattern:
            positions.append(i)
    return positions


def Reverse(Pattern):
    return Pattern[::-1]


def Complement(Pattern):
    pairs = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    complementary = ""
    for i in Pattern:
        complementary += pairs[i]
    return complementary


def ReverseComplement(Pattern):
    return Reverse(Complement(Pattern))


### Week two

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]

    # look at the first half of Genome to compute first array value
    array[0] = Genome[0:n // 2].count(symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i - 1]
        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i - 1] == symbol:
            array[i] = array[i] - 1
        if ExtendedGenome[i + (n // 2) - 1] == symbol:
            array[i] = array[i] + 1
    return array


# Analyze E. coli genome
# with open('e_coli_genome.txt', 'r') as file:
#     content = file.read()
#     print(SymbolArray(content, "C")[0::10])


def SkewArray(Genome):
    skew = [0]
    weights = {"C": -1, "G": 1, "A": 0, "T": 0}
    for nucleotide in Genome:
        skew.append(skew[-1] + weights[nucleotide])
    return skew


def MinSkew(Genome):
    skew = SkewArray(Genome)
    minimum = min(skew)
    min_positions = [i for i, v in enumerate(skew) if v == minimum]
    return min_positions


def HammingDistance(p, q):
    distance = 0
    for a, b in zip(p, q):
        if a != b:
            distance += 1
    return distance


def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    p_size = len(Pattern)
    for i in range(len(Text) - p_size + 1):
        if HammingDistance(Text[i:i + p_size], Pattern) <= d:
            positions.append(i)
    return positions


def ApproximatePatternCount(Text, Pattern, d):
    count = 0
    p_size = len(Pattern)
    for i in range(len(Text) - p_size + 1):
        if HammingDistance(Text[i:i + p_size], Pattern) <= d:
            count += 1
    return count
