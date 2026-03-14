#!/usr/bin/env bash
# =============================================================================
# extract_snps.sh — Extract medically relevant SNPs from Ancestry DNA raw data
# and format them for local AI analysis.
#
# Usage:
#   ./extract_snps.sh /path/to/AncestryDNA.txt            # full report
#   ./extract_snps.sh --compact /path/to/AncestryDNA.txt   # compact (~2k tokens)
#
# Options:
#   --compact, -c   Condensed report for small context windows (4-8k).
#                   Only found variants, one line each.
#
# Output (written to out/ with datetime stamp):
#   snp_report_YYYYMMDD_HHMMSS.txt           — Variant data
#   snp_report_prompt_YYYYMMDD_HHMMSS.txt    — LLM prompt wrapping the report
# =============================================================================

set -euo pipefail

COMPACT=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --compact|-c) COMPACT=1; shift ;;
        -*) echo "Unknown option: $1"; echo "Usage: $0 [--compact] <ancestry_raw_data.txt>"; exit 1 ;;
        *) break ;;
    esac
done

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 [--compact] <ancestry_raw_data.txt>"
    exit 1
fi

RAW_FILE="$1"
if [[ ! -f "$RAW_FILE" ]]; then
    echo "Error: File not found: $RAW_FILE"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VARIANT_FILE="$SCRIPT_DIR/variant-list.txt"
OUT_DIR="$SCRIPT_DIR/out"
mkdir -p "$OUT_DIR"

DATETIME=$(date '+%Y%m%d_%H%M%S')
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

if [[ "$COMPACT" == "1" ]]; then
    REPORT="$OUT_DIR/snp_report_compact_${DATETIME}.txt"
    PROMPT_FILE="$OUT_DIR/snp_report_compact_prompt_${DATETIME}.txt"
else
    REPORT="$OUT_DIR/snp_report_${DATETIME}.txt"
    PROMPT_FILE="$OUT_DIR/snp_report_prompt_${DATETIME}.txt"
fi

if [[ ! -f "$VARIANT_FILE" ]]; then
    echo "Error: Variant database not found: $VARIANT_FILE"
    exit 1
fi

TOTAL_SNPS=$(grep -vc "^#" "$RAW_FILE" | tr -d ' ' || echo "0")
echo "Loaded ~$TOTAL_SNPS lines from $RAW_FILE"

# ── Load variant database from file ─────────────────────────────────────────
declare -A SNP_GENE SNP_CAT SNP_DESC SNP_RISK SNP_NORMAL
RS_LIST=()
CATEGORIES=()
declare -A SEEN_CATS

while IFS='|' read -r rs gene cat desc risk normal; do
    [[ -z "$rs" ]] && continue
    RS_LIST+=("$rs")
    SNP_GENE["$rs"]="$gene"
    SNP_CAT["$rs"]="$cat"
    SNP_DESC["$rs"]="$desc"
    SNP_RISK["$rs"]="$risk"
    SNP_NORMAL["$rs"]="$normal"
    if [[ -z "${SEEN_CATS[$cat]+x}" ]]; then
        CATEGORIES+=("$cat")
        SEEN_CATS["$cat"]=1
    fi
done < <(grep -v '^\s*#' "$VARIANT_FILE" | grep -v '^\s*$')

echo "Searching for ${#RS_LIST[@]} variants across ${#CATEGORIES[@]} categories..."

# ── Look up each SNP ─────────────────────────────────────────────────────────
declare -A FOUND_GENO FOUND_CHR FOUND_POS FOUND_FLAG
FOUND_COUNT=0

for rs in "${RS_LIST[@]}"; do
    line=$(grep -P "^${rs}\t" "$RAW_FILE" 2>/dev/null | head -1 | tr -d '\r' || true)
    if [[ -n "$line" ]]; then
        chrom=$(echo "$line" | awk -F'\t' '{print $2}')
        pos=$(echo "$line" | awk -F'\t' '{print $3}')
        a1=$(echo "$line" | awk -F'\t' '{print $4}')
        a2=$(echo "$line" | awk -F'\t' '{print $5}')
        FOUND_GENO["$rs"]="${a1}${a2}"
        FOUND_CHR["$rs"]="$chrom"
        FOUND_POS["$rs"]="$pos"
        FOUND_FLAG["$rs"]="1"
        ((FOUND_COUNT++)) || true
    else
        FOUND_FLAG["$rs"]="0"
    fi
done

echo "Found $FOUND_COUNT of ${#RS_LIST[@]} queried variants."
echo ""

# ── Risk classification ──────────────────────────────────────────────────────
classify_risk() {
    local rs="$1"
    local geno="${FOUND_GENO[$rs]}"
    local risk="${SNP_RISK[$rs]}"
    local normal="${SNP_NORMAL[$rs]}"

    if [[ "$risk" == "—" || -z "$risk" ]]; then
        echo "PRESENT"
        return
    fi

    local a1="${geno:0:1}"
    local a2="${geno:1:1}"

    if [[ "$a1" == "$risk" && "$a2" == "$risk" ]]; then
        echo "HOMOZYGOUS RISK ($geno)"
    elif [[ "$a1" == "$risk" || "$a2" == "$risk" ]]; then
        echo "HETEROZYGOUS ($geno)"
    elif [[ "$a1" == "$normal" && "$a2" == "$normal" ]]; then
        echo "NORMAL ($geno)"
    else
        echo "GENOTYPE: $geno (risk allele: $risk)"
    fi
}

# ── APOE type derivation ────────────────────────────────────────────────────
determine_apoe() {
    local g1="${FOUND_GENO[rs429358]}"
    local g2="${FOUND_GENO[rs7412]}"
    if [[ "$g1" == "TT" && "$g2" == "TT" ]]; then echo "e2/e2"
    elif [[ "$g1" == "TT" && ( "$g2" == "CT" || "$g2" == "TC" ) ]]; then echo "e2/e3"
    elif [[ "$g1" == "TT" && "$g2" == "CC" ]]; then echo "e3/e3"
    elif [[ ( "$g1" == "CT" || "$g1" == "TC" ) && "$g2" == "CC" ]]; then echo "e3/e4"
    elif [[ "$g1" == "CC" && "$g2" == "CC" ]]; then echo "e4/e4"
    elif [[ ( "$g1" == "CT" || "$g1" == "TC" ) && ( "$g2" == "CT" || "$g2" == "TC" ) ]]; then echo "e2/e4"
    else echo "UNDETERMINED (rs429358=$g1, rs7412=$g2)"
    fi
}

# ── Category display label ──────────────────────────────────────────────────
cat_label() {
    echo "$1" | tr '_' ' '
}

# =============================================================================
#  WRITE REPORT
# =============================================================================

if [[ "$COMPACT" == "1" ]]; then
# ── COMPACT REPORT ───────────────────────────────────────────────────────────
{
    echo "SNP VARIANT REPORT (COMPACT)"
    echo "Generated: $TIMESTAMP  |  Source: $(basename "$RAW_FILE")"
    echo "Variants found: $FOUND_COUNT / ${#RS_LIST[@]}"
    echo ""

    for cat in "${CATEGORIES[@]}"; do
        label=$(cat_label "$cat")
        section=""
        has_entries=0

        for rs in "${RS_LIST[@]}"; do
            [[ "${SNP_CAT[$rs]}" != "$cat" ]] && continue
            [[ "${FOUND_FLAG[$rs]}" != "1" ]] && continue
            has_entries=1
            risk_text=$(classify_risk "$rs")
            section+="  $rs | ${SNP_GENE[$rs]} | ${FOUND_GENO[$rs]} | $risk_text"$'\n'
        done

        if [[ "$has_entries" == "1" ]]; then
            echo "[$label]"
            echo -n "$section"
            echo ""
        fi
    done

    # APOE derived type
    if [[ "${FOUND_FLAG[rs429358]:-0}" == "1" && "${FOUND_FLAG[rs7412]:-0}" == "1" ]]; then
        echo "[DERIVED: APOE TYPE]"
        echo "  rs429358=${FOUND_GENO[rs429358]} + rs7412=${FOUND_GENO[rs7412]} = $(determine_apoe)"
        echo ""
    fi

    # MTHFR compound status
    if [[ "${FOUND_FLAG[rs1801133]:-0}" == "1" && "${FOUND_FLAG[rs1801131]:-0}" == "1" ]]; then
        g677="${FOUND_GENO[rs1801133]}"
        g1298="${FOUND_GENO[rs1801131]}"
        r677=0; r1298=0
        [[ "${g677:0:1}" == "T" || "${g677:1:1}" == "T" ]] && r677=1
        [[ "${g1298:0:1}" == "G" || "${g1298:1:1}" == "G" ]] && r1298=1
        echo "[DERIVED: MTHFR COMPOUND STATUS]"
        echo "  C677T (rs1801133): $g677  |  A1298C (rs1801131): $g1298"
        if [[ "$g677" == "TT" ]]; then
            echo "  Classification: Homozygous C677T"
        elif [[ "$r677" == "1" && "$r1298" == "1" ]]; then
            echo "  Classification: Compound heterozygote (C677T + A1298C)"
        elif [[ "$r677" == "1" ]]; then
            echo "  Classification: Heterozygous C677T"
        elif [[ "$r1298" == "1" ]]; then
            echo "  Classification: Heterozygous A1298C"
        else
            echo "  Classification: No risk alleles detected"
        fi
        echo ""
    fi

} > "$REPORT"

else
# ── FULL REPORT ──────────────────────────────────────────────────────────────
{
    echo "============================================================================="
    echo "  SNP VARIANT REPORT"
    echo "  Generated: $TIMESTAMP"
    echo "  Source:    $(basename "$RAW_FILE")"
    echo "  Variants found: $FOUND_COUNT / ${#RS_LIST[@]}"
    echo "============================================================================="
    echo ""

    for cat in "${CATEGORIES[@]}"; do
        label=$(cat_label "$cat")

        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "  $label"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""

        for rs in "${RS_LIST[@]}"; do
            [[ "${SNP_CAT[$rs]}" != "$cat" ]] && continue

            if [[ "${FOUND_FLAG[$rs]}" == "1" ]]; then
                risk_text=$(classify_risk "$rs")
                echo "  $rs  |  ${SNP_GENE[$rs]}  |  Chr${FOUND_CHR[$rs]}:${FOUND_POS[$rs]}"
                echo "  Genotype: ${FOUND_GENO[$rs]}  |  Status: $risk_text"
                echo "  Risk allele: ${SNP_RISK[$rs]}  |  Normal allele: ${SNP_NORMAL[$rs]}"
                echo "  ${SNP_DESC[$rs]}"
                echo ""
            else
                echo "  $rs  |  ${SNP_GENE[$rs]}  |  NOT ON CHIP"
                echo "  ${SNP_DESC[$rs]}"
                echo ""
            fi
        done
    done

    # APOE derived type
    if [[ "${FOUND_FLAG[rs429358]:-0}" == "1" && "${FOUND_FLAG[rs7412]:-0}" == "1" ]]; then
        echo "============================================================================="
        echo "  DERIVED: APOE TYPE"
        echo "============================================================================="
        echo ""
        echo "  rs429358=${FOUND_GENO[rs429358]}  rs7412=${FOUND_GENO[rs7412]}  =  $(determine_apoe)"
        echo ""
    fi

    # MTHFR compound status
    if [[ "${FOUND_FLAG[rs1801133]:-0}" == "1" && "${FOUND_FLAG[rs1801131]:-0}" == "1" ]]; then
        g677="${FOUND_GENO[rs1801133]}"
        g1298="${FOUND_GENO[rs1801131]}"
        r677=0; r1298=0
        [[ "${g677:0:1}" == "T" || "${g677:1:1}" == "T" ]] && r677=1
        [[ "${g1298:0:1}" == "G" || "${g1298:1:1}" == "G" ]] && r1298=1
        echo "============================================================================="
        echo "  DERIVED: MTHFR COMPOUND STATUS"
        echo "============================================================================="
        echo ""
        echo "  C677T (rs1801133): $g677  |  A1298C (rs1801131): $g1298"
        if [[ "$g677" == "TT" ]]; then
            echo "  Classification: Homozygous C677T"
        elif [[ "$r677" == "1" && "$r1298" == "1" ]]; then
            echo "  Classification: Compound heterozygote (C677T + A1298C)"
        elif [[ "$r677" == "1" ]]; then
            echo "  Classification: Heterozygous C677T"
        elif [[ "$r1298" == "1" ]]; then
            echo "  Classification: Heterozygous A1298C"
        else
            echo "  Classification: No risk alleles detected"
        fi
        echo ""
    fi

} > "$REPORT"

fi

echo "============================================"
echo "Report written to: $REPORT"
echo "============================================"
echo ""

# ── Generate LLM prompt ──────────────────────────────────────────────────────
{
    cat << 'PROMPTEOF'
I have my genetic SNP data from an Ancestry DNA test. Below is an extracted
report of medically relevant variants with genotypes and risk classifications.

Please analyze these results and provide a comprehensive health interpretation.

Key for reading the data:
- "NOT ON CHIP" = not tested by this platform, NOT a negative result
- "HOMOZYGOUS RISK" = both copies carry the risk allele
- "HETEROZYGOUS" = one risk allele + one normal allele
- "NORMAL" = no risk alleles at that position
- This is SNP chip data, not whole genome sequencing

--- BEGIN SNP REPORT ---

PROMPTEOF
    cat "$REPORT"
    echo ""
    echo "--- END SNP REPORT ---"
} > "$PROMPT_FILE"

echo "LLM prompt written to: $PROMPT_FILE"
echo ""
echo "  Feed the prompt to your local AI:"
echo "    cat $PROMPT_FILE | ollama run llama3:70b"
if [[ "$COMPACT" == "0" ]]; then
echo ""
echo "  TIP: For smaller context windows, re-run with:"
echo "    $0 --compact $RAW_FILE"
fi
echo ""
