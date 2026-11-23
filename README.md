# HIVIntact

**Find the genetically intact HIV-1 proviruses in a set of near full-length genomes.**

Most proviruses that persist on antiretroviral therapy carry lethal defects. The rare intact ones make up the replication-competent reservoir, the target of HIV cure research. HIVIntact aligns each sequence to a subtype reference, tests it for every major class of defect, and sorts it into intact or defective with the reasons attached.

It is a command-line tool with FASTA in and FASTA + JSON out, built to drop into high-throughput sequencing pipelines.

> Wright IA, Bale MJ, Shao W, Hu W-S, Coffin JM, Van Zyl GU, Kearney MF.
> **HIVIntact: a python-based tool for HIV-1 genome intactness inference.**
> *Retrovirology* 18, 16 (2021). [doi:10.1186/s12977-021-00561-5](https://doi.org/10.1186/s12977-021-00561-5)

Validated on the Proviral Sequence Database: HIVIntact called 623 genomes intact where the database's own annotation gives 624.

## What it checks

| Defect | Reported as | Default |
|---|---|---|
| Missing, misplaced, truncated or frameshifted *gag*, *pol*, *env* | `WrongORFNumber`, `MisplacedORF`, `DeletionInOrf`, `FrameshiftInOrf` | on |
| Packaging signal (Ψ) deletion | `PackagingSignalDeletion` | on |
| Rev response element deletion | `RevResponseElementDeletion` | on |
| Major splice donor site mutation | `MajorSpliceDonorSiteMutated` | on |
| APOBEC3G/F hypermutation (Fisher's exact test, after Rose & Korber 2000) | `APOBECHypermutationDetected` | on |
| Defects in *vif*, *vpr*, *tat*, *rev*, *vpu*, *nef* | same ORF codes | reported; counted with `--include-small-orfs` |

Reverse-complemented input is detected and flipped automatically. Every check has an on/off switch; see `proviral intact --help`.

## Install

Requires Python 3 and [MAFFT](https://mafft.cbrc.jp/alignment/software/) on your `PATH` (`apt install mafft`, `brew install mafft`, or `conda install -c bioconda mafft`).

```
git clone https://github.com/ramics/HIVIntact
cd HIVIntact
python3 -m venv env && source env/bin/activate
pip install .
```

## Run

```
proviral intact --subtype B sequences.fasta
```

Subtypes: `A1`, `A2`, `B`, `C`, `D`, `F1`, `F2`, `G`, `H`.

Four files land in the current directory:

| File | Contents |
|---|---|
| `intact.fasta` | Sequences that passed every enabled check |
| `nonintact.fasta` | Sequences with at least one defect |
| `errors.json` | The defects found in each sequence |
| `orfs.json` | ORF locations for every sequence, in HXB2 coordinates |

## Pairs with

[NanoHIV](https://github.com/ramics/NanoHIV) builds the near full-length consensus genomes from Oxford Nanopore reads; HIVIntact tells you which of them are intact.
