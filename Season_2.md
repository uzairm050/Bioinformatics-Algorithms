# 🧬 Season 2 — Adding Steel Rods as Skeleton

> *"Season 1 taught me to find patterns. Season 2 taught me to understand the structure beneath them."*

Season 2 moves beyond simple pattern matching into **array-based analysis** and **distance metrics**. These tools are essential for understanding genome-wide structure and biological variation — the kind of problems that come up constantly in real NGS analysis.

---

## 2.1 — Symbol Array in DNA Sequence

### 🔬 Biological Purpose
A symbol array tracks **how many times a specific nucleotide has appeared** up to each position in a DNA sequence. This running count helps reveal regions of the genome that are enriched or depleted for a particular nucleotide — which has implications for replication and gene expression.

### 💻 Code

```python
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome  # circular genome trick
    count = PatternCount(Genome, symbol)
    array[0] = count
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] -= 1
        if ExtendedGenome[i + n - 1] == symbol:
            array[i] += 1
    return array

# Example
genome = "AGATCGGCTGA"
print(SymbolArray(genome, 'G'))
```

### 🧠 Understanding Developed
The most important thing here was the **circular genome trick** — concatenating the genome with itself (`Genome + Genome`). Bacterial chromosomes are circular, so when you reach the end you wrap back to the beginning. This was my first exposure to how biological reality (circular chromosomes) forces you to modify your algorithms. I also started seeing how maintaining a **running count** (rather than recounting from scratch each time) is a core efficiency strategy.

---

## 2.2 — Extended Symbol Array (Faster Version)

### 🔬 Biological Purpose
Same purpose as above — but optimized. The basic `SymbolArray` recalculates counts repeatedly. The extended version uses a **sliding update** approach, making it significantly faster for large genomes.

### 💻 Code

```python
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome
    array[0] = PatternCount(Genome[:n//2], symbol)
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] -= 1
        if ExtendedGenome[i + n//2 - 1] == symbol:
            array[i] += 1
    return array
```

### 🧠 Understanding Developed
Seeing two implementations of the same algorithm side by side was hugely valuable. The faster version uses the idea of **incremental updates** — instead of counting from scratch at every position, you adjust the previous count by +1 or -1 based on what entered and left the window. This is the same concept used in **sliding window algorithms** across all of computer science. It also introduced me to thinking about **time complexity** — the faster version is O(n) vs O(n²) for the naive approach.

---

## 2.3 — Skew Array in DNA Sequence

### 🔬 Biological Purpose
The **skew** measures the running difference between the count of G and C nucleotides along a genome. Due to mutation patterns during DNA replication, G tends to accumulate on one strand. The point where skew is **minimum** often corresponds to the **origin of replication** — one of the most important locations in any genome.

### 💻 Code

```python
# Method 1 — Direct calculation
def SkewArray(Genome):
    skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            skew.append(skew[i] + 1)
        elif Genome[i] == 'C':
            skew.append(skew[i] - 1)
        else:
            skew.append(skew[i])
    return skew

# Method 2 — Using dictionary mapping
def SkewArray_v2(Genome):
    score = {'A': 0, 'T': 0, 'G': 1, 'C': -1}
    skew = [0]
    for nucleotide in Genome:
        skew.append(skew[-1] + score[nucleotide])
    return skew

# Example
print(SkewArray("CATGGGCATCGGCCATACGCC"))
```

### 🧠 Understanding Developed
The Skew Array was one of the most biologically meaningful algorithms in this whole series. The idea that a simple G-C count difference can **locate the origin of replication** in a bacterial genome blew my mind — this is real, published science. Method 2 was more elegant: using a dictionary to map nucleotides to scores is cleaner and more Pythonic. I now use this pattern constantly when encoding biological symbols as numerical values.

---

## 2.4 — Minimum Skew

### 🔬 Biological Purpose
The **minimum skew position** is where the genome transitions from the lagging strand to the leading strand — in other words, it marks the **origin of replication (oriC)**. This function finds that position.

### 💻 Code

```python
def MinimumSkew(Genome):
    skew = SkewArray(Genome)
    min_skew = min(skew)
    positions = []
    for i in range(len(skew)):
        if skew[i] == min_skew:
            positions.append(i)
    return positions

# Example
print(MinimumSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGAAATGAAAGGAGTCCATGCAATCAAGCCATGCCGAAATGCGGCCCCGCATCTGGCATGGCGGCGACCGCATGCAAATGCAAACCACCGAAATGCATGATGAAATGCGCAATGCATGATGCATGATGATGCATGATGATGCATGATGATGAAATGCATGATGATGATGATGATGCATGATGATGATGAAATGCG"))
```

### 🧠 Understanding Developed
This function combined everything from Season 2 so far. The pattern of **computing an array first, then finding the minimum** is a recurring design pattern in bioinformatics. I also learned that the minimum skew can occur at **multiple positions** — biological systems don't always have a single clean answer, and algorithms must handle that ambiguity.

---

## 2.5 — Hamming Distance Between DNA Sequences

### 🔬 Biological Purpose
**Hamming Distance** measures how different two same-length DNA sequences are — specifically, how many positions have different nucleotides. This directly models **point mutations** between sequences and is fundamental to phylogenetics, variant calling, and approximate motif finding.

### 💻 Code

```python
def HammingDistance(p, q):
    distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance

# Example
print(HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC"))  # Output: 3
```

### 🧠 Understanding Developed
This is beautifully simple but biologically profound. The Hamming Distance is how we quantify **evolutionary divergence** between sequences at the most basic level. It also introduced me to the concept of a **metric** — a formal measure of distance — which appears throughout bioinformatics in clustering, alignment scoring, and phylogenetic tree building. In variant calling, every SNP is essentially a Hamming distance of 1.

---

## 2.6 — Approximate Pattern Matching

### 🔬 Biological Purpose
In real biological data, patterns are **never perfectly conserved**. A transcription factor binding site may tolerate 1-2 mutations. Approximate pattern matching finds all positions where a pattern occurs with **at most d mismatches** — far more biologically realistic than exact matching.

### 💻 Code

```python
def ApproximatePatternMatching(Pattern, Text, d):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

# Example
print(ApproximatePatternMatching("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCGGCGGACTATGGACTGAACGTAAAGGCGTATGTAAACAACCGCCCCGATTCCAG", 3))
```

### 🧠 Understanding Developed
This was the most impactful algorithm in Season 2 for me. Combining the sliding window from Season 1 with Hamming Distance from earlier in Season 2 created something genuinely useful — a basic **fuzzy search** for biological sequences. This is conceptually what tools like BLAST do at a much larger scale. It also reinforced that in biology, **tolerance for variation is the rule, not the exception**.

---

## ✅ Season 2 Summary

| Algorithm | Core Concept | Biological Application |
|-----------|-------------|----------------------|
| Symbol Array | Running count | Nucleotide enrichment analysis |
| Extended Symbol Array | Sliding window update (O(n)) | Efficient genome-wide scanning |
| Skew Array | G-C running difference | Replication strand bias |
| Minimum Skew | Argmin of array | Origin of replication (oriC) |
| Hamming Distance | Position-wise mismatch count | Mutation quantification, variant calling |
| Approximate Pattern Matching | Fuzzy sliding window search | Motif finding with mismatches |

**Key Takeaway:** Real biology is messy. Sequences mutate, genomes are circular, and patterns are approximate. Season 2 taught me to build algorithms that reflect biological reality rather than assuming perfect, clean data.
