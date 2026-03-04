# UBQ Gene Family Identification Pipeline (PF00240)

This repository contains a reproducible pipeline to identify UBQ family members from a proteome using:
1) Pfam HMM search (PF00240 ubiquitin domain)
2) merged-repeat counting to detect polyubiquitin
3) BLASTp validation against UniProt Swiss-Prot (reviewed) for single-domain candidates

## What is included
- `scripts/pfam_pipeline.sh`: utilities to parse Pfam domtblout and extract FASTA
- `run_ubq.sh`: end-to-end UBQ identification
- `config/UBQ.yml`: rules used for UBQ family
- `example/`: example outputs from a real run (final IDs/FASTA)

## What is NOT included (too large)
- Pfam-A.hmm
- UniProt Swiss-Prot fasta and BLAST database files

You should download these resources locally and provide paths when running the pipeline.

## Requirements
- HMMER (hmmscan)
- BLAST+ (makeblastdb, blastp)
- seqkit
- awk, sort, comm (coreutils)

## Quick start
Prepare:
- proteome FASTA: e.g. `data/DB.protein.fa`
- Pfam HMM: `hmm/Pfam-A.hmm`
- Swiss-Prot BLAST DB prefix: `db/uniprot/sprot` (created by `makeblastdb`)

Run:
```bash
bash run_ubq.sh data/DB.protein.fa hmm/Pfam-A.hmm db/uniprot/sprot
