# 🧬 Season 1 — Digging Down Deep

> *"Before you can run any complex pipeline, you need to understand what a pattern in a sequence even means."*

Season 1 covers the most fundamental operations in bioinformatics — treating DNA as a string and extracting meaningful biological information from it using core Python logic. These algorithms might look simple, but they are the backbone of tools like BLAST, genome assemblers, and replication origin finders.

---

## 1.1 — Pattern Count

### 🔬 Biological Purpose
In molecular biology, certain short DNA sequences (called **k-mers**) appear more frequently near the **origin of replication** — the region where DNA copying begins. Counting how often a pattern appears helps us locate these biologically significant regions.

### 💻 Code

```python
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

# Example
dna = "AGATCGAGATCG"
pattern = "AGATC"
print(PatternCount(dna, pattern))  # Output: 2
```

### 🧠 Understanding Developed
This was my first real encounter with **sliding window** logic — moving a fixed-size window across a sequence one step at a time. The key insight is `range(len(Text) - len(Pattern) + 1)` — without the `+1` correction, you miss the last valid position. Simple, but it forces you to think carefully about **off-by-one errors**, which show up everywhere in bioinformatics.

---

## 1.2 — Frequency Map of DNA Sequence

### 🔬 Biological Purpose
Rather than searching for one specific pattern, a frequency map gives us a **complete picture** of all k-mers present in a sequence and how often each one appears. This is essential for identifying overrepresented sequences that may have functional significance.

### 💻 Code

```python
def FrequencyMap(Text, k):
    freq = {}
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i+k]
        if kmer in freq:
            freq[kmer] += 1
        else:
            freq[kmer] = 1
    return freq

# Example
dna = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
print(FrequencyMap(dna, 4))
```

### 🧠 Understanding Developed
This introduced me to **Python dictionaries** as a tool for counting — a pattern I now use constantly. Instead of searching for one pattern at a time, we scan once and build a lookup table. This is fundamentally more efficient and is the same idea behind **hash tables** used in genome indexing tools like BWA and Bowtie.

---

## 1.3 — Most Frequent K-mers

### 🔬 Biological Purpose
Building on the frequency map, this function finds the **most overrepresented k-mers** in a sequence. In bacterial genomes, the most frequent 9-mer near the replication origin is often the **DnaA box** — a sequence that proteins bind to initiate DNA replication.

### 💻 Code

```python
def FrequentWords(Text, k):
    freq = FrequencyMap(Text, k)
    max_count = max(freq.values())
    frequent_kmers = []
    for kmer in freq:
        if freq[kmer] == max_count:
            frequent_kmers.append(kmer)
    return frequent_kmers

# Example
dna = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
print(FrequentWords(dna, 4))  # Output: ['CATG', 'GCAT']
```

### 🧠 Understanding Developed
This taught me to **build on top of previous functions** rather than rewriting everything. The `FrequencyMap` does the heavy lifting; `FrequentWords` just extracts the maximum. I also learned that there can be **multiple equally frequent k-mers** — biological signals are rarely unique, and code must handle that gracefully.

---

## 1.4 — Complementary DNA String

### 🔬 Biological Purpose
DNA is **double-stranded**. Every strand has a complementary strand running in the opposite direction (5'→3' becomes 3'→5'). When searching for a pattern like a DnaA box, we must search **both strands** — so generating the reverse complement is essential.

### 💻 Code

```python
def ReverseComplement(Pattern):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp = ''
    for nucleotide in reversed(Pattern):
        reverse_comp += complement[nucleotide]
    return reverse_comp

# Example
print(ReverseComplement("ATGCCA"))  # Output: TGGCAT
```

### 🧠 Understanding Developed
This function taught me two things: first, the biological rule of **Watson-Crick base pairing** (A↔T, C↔G); second, how to use a **dictionary as a lookup table** for transformations. The `reversed()` step captures the antiparallel nature of DNA. This concept is used in virtually every sequence alignment tool — you always search both strands.

---

## 1.5 — Pattern Matching

### 🔬 Biological Purpose
Once we know our pattern of interest (e.g., a DnaA box), we need to find **every position** in the genome where it occurs. This is the basis of sequence search tools and motif scanning.

### 💻 Code

```python
def PatternMatching(Pattern, Genome):
    positions = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

# Example
genome = "AATATATATAT"
pattern = "ATAT"
print(PatternMatching(pattern, genome))  # Output: [0, 2, 4, 6]
```

### 🧠 Understanding Developed
This looks almost identical to `PatternCount` — but the **output is different**: positions instead of a count. This reinforced an important lesson: small changes in what you return can completely change how useful a function is. Returning positions lets you **visualize where** a motif occurs across a genome, which is far more informative than just knowing it occurs 4 times.

---

## ✅ Season 1 Summary

| Algorithm | Core Concept | Biological Application |
|-----------|-------------|----------------------|
| Pattern Count | Sliding window | Replication origin detection |
| Frequency Map | Dictionary/hash table | K-mer profiling |
| Most Frequent K-mers | Argmax over dictionary | DnaA box finding |
| Reverse Complement | Lookup table + reversal | Double-strand searching |
| Pattern Matching | Sliding window + position logging | Motif location mapping |

**Key Takeaway:** DNA is a string. String operations are biology. Everything in Season 1 is a building block used — often invisibly — inside every modern bioinformatics tool.
