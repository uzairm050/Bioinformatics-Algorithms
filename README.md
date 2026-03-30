# 🧬 Bioinformatics Algorithms — A Learning Journey

---

## 📌 About This Repository

This repository documents my personal, hands-on journey through core bioinformatics algorithms implemented in Python. Structured across **four progressive seasons**, it starts from the most basic DNA pattern operations and builds all the way up to probabilistic motif-finding algorithms used in real computational biology research.

My background is in Biotechnology (BSc, University of Peshawar), and I am currently pursuing my MS in Bioinformatics at NUST. Working through these algorithms helped me bridge the gap between my biological knowledge and computational thinking — understanding not just *what* tools like BLAST or MEME do, but *how* and *why* they work.

> **This is not just a code dump — each season contains detailed explanations, biological context, and honest reflections on what I learned and struggled with.**

---

## 🗺️ Repository Structure

```
📁 Bioinformatics-Algorithms/
│
├── 📄 README.md               ← You are here
├── 📄 Season_1.md             ← Pattern counting, k-mers, reverse complement
├── 📄 Season_2.md             ← Skew arrays, Hamming distance, approximate matching
├── 📄 Season_3.md             ← Profile matrices, motif scoring, Greedy Motif Search
└── 📄 Season_4.md             ← Pseudocounts, Randomized Motif Search
```

---

## 🧬 Season Overview

### [🔗 Season 1 — Digging Down Deep](./Season_1.md)
> *Foundational DNA string operations*

The starting point. Treating DNA as a string and building the core operations that every bioinformatics pipeline relies on under the hood.

| Topic | Key Concept |
|-------|------------|
| Pattern Count | Sliding window |
| Frequency Map | Dictionary/hash table |
| Most Frequent K-mers | K-mer overrepresentation |
| Reverse Complement | Watson-Crick base pairing |
| Pattern Matching | Position logging |

---

### [🔗 Season 2 — Adding Steel Rods as Skeleton](./Season_2.md)
> *Array-based analysis and distance metrics*

Moving from simple pattern matching to understanding the structure and variation within DNA sequences. Introduces the biological concept of the replication origin.

| Topic | Key Concept |
|-------|------------|
| Symbol Array | Running nucleotide count |
| Extended Symbol Array | O(n) sliding window update |
| Skew Array | G-C imbalance tracking |
| Minimum Skew | Origin of replication (oriC) |
| Hamming Distance | Mutation quantification |
| Approximate Pattern Matching | Fuzzy biological search |

---

### [🔗 Season 3 — Adding Concrete Slurry to the Structure](./Season_3.md)
> *Modular algorithm design and motif finding*

The most complex season. Individual functions are built step by step and assembled into **Greedy Motif Search** — a complete algorithm for finding regulatory motifs across multiple DNA sequences.

| Topic | Key Concept |
|-------|------------|
| Count Matrix | Positional frequency |
| Profile Matrix | Position Weight Matrix (PWM) |
| Consensus Sequence | Motif summarization |
| Score | Motif quality metric |
| Profile Most Probable K-mer | Likelihood computation |
| Greedy Motif Search | Complete motif finder |
| Patterns with Mismatches | Recursive neighborhood |

---

### [🔗 Season 4 — Colouring the Building and Final Touch-ups](./Season_4.md)
> *Probabilistic refinement and stochastic optimization*

Addressing the weaknesses of greedy algorithms. Introduces Laplace smoothing and randomization — concepts that connect directly to modern machine learning.

| Topic | Key Concept |
|-------|------------|
| Greedy Search + Pseudocounts | Laplace smoothing |
| Randomized Motif Search | Stochastic optimization |
| Repeated Random Search | Random restarts |

---

## 🧠 Key Takeaways Across All Seasons

**1. DNA is data.** The moment you start treating a genome as a string and applying algorithmic thinking to it, a huge space of tools and techniques becomes available.

**2. Build modular.** Every complex algorithm in this repo is just a composition of smaller, testable functions. This is good engineering and good science.

**3. Biology is approximate.** Exact matching is rarely biologically realistic. Hamming distance, pseudocounts, and approximate matching all reflect the messy reality of biological sequences.

**4. Bioinformatics and ML share the same math.** Laplace smoothing, probabilistic models, stochastic optimization, random restarts — these appear in both fields because they're solving the same class of problems.

**5. Understand the tool before using the tool.** Working through these algorithms from scratch gave me an intuition for what BLAST, MEME, and BioPython are actually doing — which makes me a better scientist.

---

## 🔭 What's Next

Building on this foundation, I am currently working on:
- 🔬 RNA-seq differential expression analysis (HISAT2, DESeq2)
- 🧫 Single-cell RNA-seq analysis (Scanpy, Seurat)
- 🤖 Machine learning applications in genomics
- 📊 Variant calling pipelines (GATK)

---

## 📚 Resources That Helped

- [Bioinformatics Algorithms (Compeau & Pevzner)](https://www.bioinformaticsalgorithms.org/) — the textbook behind most of these algorithms
- [Rosalind](http://rosalind.info) — great platform for practicing bioinformatics problems
- [BioPython Documentation](https://biopython.org/docs/latest/api/)

---


