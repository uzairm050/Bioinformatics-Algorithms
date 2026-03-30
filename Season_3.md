# 🧬 Season 3 — Adding Concrete Slurry to the Structure

> *"Season 3 is where individual functions stop being exercises and start becoming a real algorithm."*

Season 3 introduces **motif finding** — one of the central unsolved problems in computational biology. The approach here is deliberately modular: small, focused functions are built one by one, and then assembled into **Greedy Motif Search**, a complete algorithm for finding shared patterns across multiple DNA sequences.

---

## 3.1 — Basic NumPy Array Understanding

### 🔬 Biological Purpose
NumPy is the numerical backbone of almost all bioinformatics in Python — from storing count matrices to running statistical tests. Before building a profile matrix, understanding how NumPy arrays work is essential.

### 💻 Code

```python
import numpy as np

# Creating arrays
arr = np.array([1, 2, 3, 4, 5])

# 2D array (like a matrix)
matrix = np.array([[1, 2, 3],
                   [4, 5, 6]])

# Basic operations
print(arr.shape)        # (5,)
print(matrix.shape)     # (2, 3)
print(np.sum(matrix, axis=0))   # Column sums: [5, 7, 9]
print(np.sum(matrix, axis=1))   # Row sums: [6, 15]
```

### 🧠 Understanding Developed
NumPy's concept of **axis** was initially confusing — `axis=0` means "collapse rows" (work down columns), and `axis=1` means "collapse columns" (work across rows). Once this clicked, it completely changed how I think about matrix operations. A profile matrix in bioinformatics is just a 2D array where **rows = nucleotides** and **columns = positions** — and NumPy makes operating on it trivial.

---

## 3.2 — Counting Nucleotide Frequencies in DNA

### 🔬 Biological Purpose
Given a collection of DNA strings (motifs), we need to count how often each nucleotide (A, C, G, T) appears at each position. This **count matrix** is the raw material for building a profile matrix and finding consensus sequences.

### 💻 Code

```python
def CountMatrix(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = [0] * k
    for motif in Motifs:
        for i, nucleotide in enumerate(motif):
            count[nucleotide][i] += 1
    return count

# Example
motifs = ["AACGTA",
          "CCCGTT",
          "CACCTT",
          "GGATTA",
          "TTTTCA"]
print(CountMatrix(motifs))
```

### 🧠 Understanding Developed
Using `enumerate()` to get both the index and the character simultaneously was a Python technique I started using everywhere after this. The count matrix captures **positional nucleotide frequency** — not just how many G's are in a sequence, but how many G's are at **position 3** across all motifs. This positional thinking is fundamental to understanding how transcription factors recognize binding sites.

---

## 3.3 — Profile Matrix from DNA

### 🔬 Biological Purpose
A **profile matrix** converts raw counts into **probabilities**. Each cell tells us: "given that we're at position i, what is the probability of seeing nucleotide X?" This is how we model a **position weight matrix (PWM)** — the standard representation of transcription factor binding sites in databases like JASPAR.

### 💻 Code

```python
def ProfileMatrix(Motifs):
    count = CountMatrix(Motifs)
    profile = {}
    t = len(Motifs)
    for symbol in "ACGT":
        profile[symbol] = [x / t for x in count[symbol]]
    return profile

# Example
motifs = ["AACGTA",
          "CCCGTT",
          "CACCTT",
          "GGATTA",
          "TTTTCA"]
print(ProfileMatrix(motifs))
```

### 🧠 Understanding Developed
Dividing by the total number of motifs (`t`) converts counts to probabilities — this is the **normalization** step that appears in virtually every statistical model in bioinformatics. The profile matrix is essentially a simple **probabilistic model** of a sequence motif. I recognized this as the same principle behind Naive Bayes classifiers in machine learning — assign probabilities to each feature at each position.

---

## 3.4 — Consensus Sequence from DNA

### 🔬 Biological Purpose
The **consensus sequence** is the "most representative" motif — at each position, we pick the nucleotide with the highest probability. For a transcription factor binding site, this is the optimal binding sequence. It's used as a compact summary of a motif.

### 💻 Code

```python
def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountMatrix(Motifs)
    consensus = ""
    for i in range(k):
        max_count = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][i] > max_count:
                max_count = count[symbol][i]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus

# Example
motifs = ["AACGTA",
          "CCCGTT",
          "CACCTT",
          "GGATTA",
          "TTTTCA"]
print(Consensus(motifs))  # Output: CACCTA
```

### 🧠 Understanding Developed
The consensus sequence taught me about **information loss through summarization** — by picking just one nucleotide per position, we lose all the uncertainty information. This is why PWMs (which keep the full probability distribution) are preferred in real tools. But the consensus is still useful as a quick summary and is what appears in papers when researchers describe a binding motif.

---

## 3.5 — Scoring DNA Motifs

### 🔬 Biological Purpose
A motif score measures how **far a set of motifs is from a perfect consensus**. Lower score = more conserved motifs = stronger biological signal. This is used to evaluate and compare candidate motif sets during motif search.

### 💻 Code

```python
def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    for motif in Motifs:
        for i in range(len(motif)):
            if motif[i] != consensus[i]:
                score += 1
    return score

# Example
motifs = ["AACGTA",
          "CCCGTT",
          "CACCTT",
          "GGATTA",
          "TTTTCA"]
print(Score(motifs))
```

### 🧠 Understanding Developed
Score is essentially the **total Hamming distance** from all motifs to the consensus. I started seeing how everything in Season 1 and 2 connects — Hamming Distance from Season 2 is now being used inside a scoring function for motif finding. This cumulative building of functions is exactly how real bioinformatics software is structured. Also: minimizing this score is an **optimization problem**, which directly connects to machine learning concepts.

---

## 3.6 — Profile Most Probable K-mer

### 🔬 Biological Purpose
Given a profile matrix (representing a known or candidate motif), this function finds the k-mer in a new DNA sequence that is **most likely to have been generated by that profile**. This is how we extend a motif model to new sequences.

### 💻 Code

```python
def ProfileMostProbableKmer(Text, k, Profile):
    max_prob = -1
    best_kmer = Text[:k]
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i+k]
        prob = 1.0
        for j, nucleotide in enumerate(kmer):
            prob *= Profile[nucleotide][j]
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer

# Example
Text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCTTGGGGGTTTATTATCAAAGATTGATGATGATAATTTCCTGCATGAAATGGCTTCAAGCTTTATAAGATCAAGGCAATTTGCCAATTTTTCATAAATAAATAATGAGATCAATCTGAAATGAGCGTGAAGCAAGAGGGAGATGGGCATGGCAAACCACCTCAAATGGCATCCTCAACAAAATGATTTAAAATATAAATTTGAATTTTATAAATGT"
k = 12
Profile = {'A': [0.2, 0.2, 0.3, 0.2, 0.3],
           'C': [0.4, 0.3, 0.1, 0.5, 0.1],
           'G': [0.3, 0.3, 0.5, 0.2, 0.4],
           'T': [0.1, 0.2, 0.1, 0.1, 0.2]}
print(ProfileMostProbableKmer(Text, 5, Profile))
```

### 🧠 Understanding Developed
Multiplying probabilities position by position is essentially computing a **likelihood** under a probabilistic model. I recognized this as the same operation used in **Hidden Markov Models** and other sequence models in bioinformatics. One important limitation I noticed: if any probability is zero, the entire product becomes zero — which is a real problem. This sets up the need for pseudocounts in Season 4 perfectly.

---

## 3.7 — Greedy Motif Search (The Grand Finale)

### 🔬 Biological Purpose
**Greedy Motif Search** is a complete algorithm for finding a motif of length k that is **shared across t DNA sequences** (e.g., the promoter regions of co-regulated genes). It iterates through all possible starting k-mers in the first sequence, builds a profile from them, and extends to each remaining sequence greedily.

### 💻 Code

```python
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = [sequence[:k] for sequence in Dna]
    
    for i in range(len(Dna[0]) - k + 1):
        Motifs = [Dna[0][i:i+k]]
        
        for j in range(1, t):
            profile = ProfileMatrix(Motifs)
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))
        
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    
    return BestMotifs

# Example
Dna = [
    "GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG"
]
k = 3
t = 5
print(GreedyMotifSearch(Dna, k, t))
```

### 🧠 Understanding Developed
This was the most satisfying moment of the entire repository — watching every function from Season 1 through Season 3 come together into one meaningful algorithm. `CountMatrix` → `ProfileMatrix` → `ProfileMostProbableKmer` → `Score` — all assembled into `GreedyMotifSearch`. The greedy strategy (always pick the locally best option) is fast and intuitive, but I understood its key weakness: it can get **stuck in local optima** and miss the globally best motif set. This limitation is precisely why Season 4 introduces randomization.

---

## 3.8 — Finding Patterns with Mismatches (Bonus)

### 🔬 Biological Purpose
This bonus algorithm finds the most frequent k-mers in a sequence while **tolerating up to d mismatches** — the approximate version of the frequent k-mer problem. This is more biologically realistic since binding sites and functional sequences often have minor variations.

### 💻 Code

```python
def Neighbors(Pattern, d):
    if d == 0:
        return {Pattern}
    if len(Pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    neighborhood = set()
    SuffixNeighbors = Neighbors(Pattern[1:], d)
    for text in SuffixNeighbors:
        if HammingDistance(Pattern[1:], text) < d:
            for nucleotide in 'ACGT':
                neighborhood.add(nucleotide + text)
        else:
            neighborhood.add(Pattern[0] + text)
    return neighborhood

def FrequentWordsWithMismatches(Text, k, d):
    patterns = []
    freq = {}
    for i in range(len(Text) - k + 1):
        neighborhood = Neighbors(Text[i:i+k], d)
        for neighbor in neighborhood:
            if neighbor not in freq:
                freq[neighbor] = 1
            else:
                freq[neighbor] += 1
    max_count = max(freq.values())
    for pattern in freq:
        if freq[pattern] == max_count:
            patterns.append(pattern)
    return patterns
```

### 🧠 Understanding Developed
The `Neighbors` function uses **recursion** — my first real encounter with recursive thinking in bioinformatics. It generates all k-mers within Hamming distance d of a given pattern by building up from the suffix. This is computationally clever and connects directly to how **sequence neighborhoods** are used in variant effect prediction and mutational scanning experiments.

---

## ✅ Season 3 Summary

| Algorithm | Core Concept | Biological Application |
|-----------|-------------|----------------------|
| Count Matrix | Positional frequency counting | Motif characterization |
| Profile Matrix | Normalization to probabilities | Position Weight Matrix (PWM) |
| Consensus Sequence | Argmax per position | Binding site summarization |
| Score | Total deviation from consensus | Motif quality evaluation |
| Profile Most Probable K-mer | Likelihood computation | Motif model extension |
| Greedy Motif Search | Greedy algorithm assembly | Regulatory motif discovery |
| Patterns with Mismatches | Recursive neighborhood | Approximate motif finding |

**Key Takeaway:** Complex algorithms are just well-organized collections of simple functions. Season 3 showed me how modular design — building small, testable functions and composing them — is both good software engineering and good science.
