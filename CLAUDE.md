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

Output files are written to the current directory: `snp_report.txt` / `snp_report_compact.txt` (variant report) and `snp_report_prompt.txt` / `snp_report_compact_prompt.txt` (LLM-ready prompt wrapping the report).

## Architecture

**Single-script tool** (`extract_snps.sh`): everything runs in one bash script using associative arrays.

- **SNP database**: Currently embedded as a pipe-delimited heredoc (`SNP_DB`) inside the script. Each entry: `rsID|Gene|Category|Description|Risk_Allele|Normal_Allele`. A planned refactor (see TODO.md) will externalize this to `variant-list.txt`.
- **variant-list.txt**: Extended variant list with additional categories (FH/LDL, Lp(a), statin safety/response, CAD risk, lipid metabolism, blood pressure, inflammation). Same pipe-delimited format. Not yet consumed by the script.
- **Lookup**: For each rsID in the database, greps the raw Ancestry file (tab-delimited: `rsID\tchrom\tpos\tallele1\tallele2`) and classifies genotype as homozygous risk, heterozygous, or normal.
- **Special interpretation**: APOE type (combined rs429358 + rs7412), MTHFR compound heterozygote detection, G6PD summary.
- **Categories**: FAVISM, MTHFR, ALZHEIMERS, CARDIOVASCULAR, CLOTTING, CANCER, PHARMACOGENOMICS, IMMUNE, METABOLISM, NUTRITION, DETOX, MENTAL_HEALTH — with plans to make these dynamic based on input data.

## Directory Layout

- `in/` — Raw Ancestry DNA data files (private, gitignored)
- `out/` — Generated reports (gitignored)

## Planned Changes (TODO.md)

Key refactors: externalize SNP_DB to variant-list.txt, dynamically derive categories from input data, output to `out/` with datetime stamps, expand variant list (including Fragile X), and minimize diagnostic text in output since an LLM will interpret results downstream.
