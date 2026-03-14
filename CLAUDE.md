# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A bash tool that extracts medically relevant SNP (Single Nucleotide Polymorphism) variants from Ancestry DNA raw data files and formats them into reports suitable for local LLM analysis.

## Usage

```bash
# Full report
./extract_snps.sh in/AncestryDNA.txt

# Compact report (~2k tokens, for small context windows)
./extract_snps.sh --compact in/AncestryDNA.txt
```

Output files are written to `out/` with datetime-stamped filenames (e.g. `snp_report_20260313_195100.txt` and the corresponding `_prompt_` file for LLM ingestion).

## Architecture

**Single-script tool** (`extract_snps.sh`): everything runs in one bash script using associative arrays.

- **SNP database** (`variant-list.txt`): Pipe-delimited file with format `rsID|Gene|Category|Description|Risk_Allele|Normal_Allele`. The script reads this at runtime. Add new variants or categories by editing this file — categories are derived dynamically from whatever appears in the data.
- **Lookup**: For each rsID in the database, greps the raw Ancestry file (tab-delimited: `rsID\tchrom\tpos\tallele1\tallele2`) and classifies genotype as HOMOZYGOUS RISK, HETEROZYGOUS, or NORMAL.
- **Derived data**: APOE type (combined rs429358 + rs7412) and MTHFR compound heterozygote status are computed and included as labeled derived sections.
- **Output philosophy**: Reports contain only factual data (genotypes, risk classifications, descriptions). No diagnostic conclusions, medical advice, or recommendations — the downstream LLM handles interpretation.

## Directory Layout

- `in/` — Raw Ancestry DNA data files (private, gitignored)
- `out/` — Generated reports (gitignored)

## TODO.md

Tracks planned features and open questions. Check it before starting new work to stay aligned with project direction.
