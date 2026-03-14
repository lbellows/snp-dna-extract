# snp-dna-extract

Extract medically relevant SNP variants from Ancestry DNA raw data and format them for analysis by a local LLM.

## Usage

```bash
# Full report — all variants with descriptions
./extract_snps.sh in/AncestryDNA.txt

# Compact report — found variants only, one line each (~2k tokens)
./extract_snps.sh --compact in/AncestryDNA.txt
```

Reports are written to `out/` with datetime-stamped filenames:
- `snp_report_YYYYMMDD_HHMMSS.txt` — variant data
- `snp_report_prompt_YYYYMMDD_HHMMSS.txt` — LLM-ready prompt wrapping the report

Feed the prompt file to your local model:
```bash
cat out/snp_report_prompt_*.txt | ollama run llama3:70b
```

## Variant Database

Variants are defined in `variant-list.txt` using pipe-delimited format:
```
rsID|Gene|Category|Description|Risk_Allele|Normal_Allele
```

Categories are derived dynamically from the data — add new categories simply by using a new category name in the file. Lines starting with `#` are comments.

## Privacy

**Your DNA data is sensitive and uniquely identifying.** Treat raw genotype files and generated reports with the same care as medical records.

- **Run inference locally.** Use a local model (Ollama, llama.cpp, LM Studio, etc.) rather than uploading reports to cloud AI services. Your genetic data cannot be changed if it is leaked.
- **Keep raw data files in `in/`** — this directory is gitignored.
- **Do not commit** raw DNA files, generated reports, or any output containing genotype data to version control.
- **Be cautious sharing reports** — even partial SNP data can reveal disease risk, carrier status, and identity.

## Disclaimer

This tool is for **educational and informational purposes only**. It is not a medical device and does not provide medical diagnoses. Genetic variants are one factor among many. Always consult a qualified genetic counselor or physician for medical decisions.
