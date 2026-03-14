#!/usr/bin/env bash
# =============================================================================
# extract_snps.sh — Extract medically relevant SNPs from Ancestry DNA raw data
# and format them for local AI analysis.
#
# Usage:
#   chmod +x extract_snps.sh
#   ./extract_snps.sh /path/to/AncestryDNA.txt            # full report
#   ./extract_snps.sh --compact /path/to/AncestryDNA.txt   # compact (~2k tokens)
#
# Options:
#   --compact, -c   Condensed report for small context windows (4-8k).
#                   Only found variants, one line each, plus key summaries.
#
# Output:
#   snp_report.txt / snp_report_compact.txt           — Variant report
#   snp_report_prompt.txt / snp_report_compact_prompt.txt — LLM prompt
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

if [[ "$COMPACT" == "1" ]]; then
    REPORT="snp_report_compact.txt"
    PROMPT_FILE="snp_report_compact_prompt.txt"
else
    REPORT="snp_report.txt"
    PROMPT_FILE="snp_report_prompt.txt"
fi
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

TOTAL_SNPS=$(grep -vc "^#" "$RAW_FILE" | tr -d ' ' || echo "0")
echo "Loaded ~$TOTAL_SNPS lines from $RAW_FILE"

# ── SNP lookup table (pipe-delimited) ────────────────────────────────────────
read -r -d '' SNP_DB << 'SNPEOF' || true
rs1050828|G6PD|FAVISM|G6PD A- variant (Val68Met) — most common cause of G6PD deficiency worldwide. Triggers hemolytic anemia on exposure to fava beans, certain drugs (primaquine, dapsone, sulfonamides), and infections. X-linked: one copy in males = deficient.|T|C
rs1050829|G6PD|FAVISM|G6PD A variant (Asn126Asp) — often found alongside rs1050828. Alone it causes mild reduction; combined with rs1050828 creates the A- deficiency phenotype.|T|C
rs5030868|G6PD|FAVISM|G6PD Mediterranean variant (Ser188Phe) — causes severe G6PD deficiency. Common in Mediterranean, Middle Eastern, and North African populations. Triggers severe hemolysis with fava beans and oxidant drugs.|T|C
rs137852328|G6PD|FAVISM|G6PD Canton variant — common in East and Southeast Asian populations. Causes moderate-to-severe deficiency.|A|G
rs76723693|G6PD|FAVISM|G6PD Mahidol variant — most common G6PD variant in Southeast Asia (especially Myanmar/Thailand). Causes mild-to-moderate deficiency.|A|G
rs2230037|G6PD|FAVISM|G6PD synonymous tag SNP — does not directly cause deficiency but tags G6PD haplotypes. Useful for inferring G6PD status when causal variants are not on the chip.|A|G
rs1801133|MTHFR|MTHFR|MTHFR C677T — the most studied MTHFR variant. TT genotype reduces enzyme activity to ~30% of normal, CT to ~65%. Leads to elevated homocysteine and reduced conversion of folic acid to active methylfolate. TT carriers often benefit from methylfolate (5-MTHF) instead of folic acid.|T|C
rs1801131|MTHFR|MTHFR|MTHFR A1298C — second most common MTHFR variant. CC genotype causes mild reduction in enzyme activity. Most impactful when compound heterozygous (one copy of C677T + one copy of A1298C together).|G|T
rs2274976|MTHFR|MTHFR|MTHFR R594Q — rare functional variant that reduces enzyme activity. Less studied than C677T/A1298C but can compound with them.|T|C
rs17367504|MTHFR region|MTHFR|Near MTHFR — GWAS hit for blood pressure regulation. The A allele is associated with lower blood pressure.|A|G
rs429358|APOE|ALZHEIMERS|APOE epsilon-4 determinant — this SNP plus rs7412 together determine your APOE type (e2/e3/e4). The C allele here contributes to the e4 allele, the strongest common genetic risk factor for late-onset Alzheimer's. e4/e4 carriers have ~8-12x increased risk.|C|T
rs7412|APOE|ALZHEIMERS|APOE epsilon-2 determinant — the T allele here contributes to the e2 allele, which is protective against Alzheimer's and associated with longevity. Must be interpreted together with rs429358.|T|C
rs1333049|9p21|CARDIOVASCULAR|9p21 locus — the strongest and most replicated common genetic risk factor for coronary artery disease. CC genotype carries ~1.6x risk vs GG. Not modifiable, but knowledge supports aggressive lifestyle management.|C|G
rs10757274|9p21|CARDIOVASCULAR|9p21 region — additional tag SNP for the chromosome 9p21 heart disease locus. GG genotype associated with ~1.3x increased CAD risk.|G|A
rs4977574|9p21|CARDIOVASCULAR|9p21 region — another independent signal at the 9p21 locus for myocardial infarction and coronary artery disease risk.|G|A
rs1801282|PPARG|CARDIOVASCULAR|PPARG Pro12Ala — the G (Ala) allele is protective against type 2 diabetes. CC (Pro/Pro) is the higher-risk genotype. Affects insulin sensitivity and fat cell differentiation.|G|C
rs7903146|TCF7L2|CARDIOVASCULAR|TCF7L2 — the single strongest common genetic risk factor for type 2 diabetes. TT genotype carries ~1.7x risk. Affects incretin signaling and insulin secretion.|T|C
rs12255372|TCF7L2|CARDIOVASCULAR|TCF7L2 secondary signal — TT genotype increases type 2 diabetes risk by ~1.4x. In linkage disequilibrium with rs7903146.|T|G
rs5186|AGTR1|CARDIOVASCULAR|Angiotensin II receptor A1166C — CC genotype associated with hypertension, aortic stiffness, and increased response to angiotensin. May influence choice of antihypertensive medication.|C|A
rs1799983|NOS3|CARDIOVASCULAR|eNOS Glu298Asp — TT genotype reduces nitric oxide production, associated with endothelial dysfunction, hypertension, and coronary spasm. Exercise and dietary nitrate can help compensate.|T|G
rs662|PON1|CARDIOVASCULAR|Paraoxonase Q192R — affects ability to break down oxidized lipids and organophosphate pesticides. AA genotype has lower enzyme activity. Relevant to cardiovascular protection and pesticide sensitivity.|A|G
rs6025|F5|CLOTTING|Factor V Leiden — the most common inherited thrombophilia in people of European descent. A single copy (AG) increases venous clot risk ~5x; two copies (AA) increase risk ~50x. Critical for decisions about oral contraceptives, HRT, surgery, and long flights.|A|G
rs1799963|F2|CLOTTING|Prothrombin G20210A — second most common inherited clotting disorder. AG genotype increases venous thrombosis risk ~3x. Compounds with Factor V Leiden. Important for contraceptive and surgical planning.|A|G
rs8176719|ABO|CLOTTING|ABO blood type determinant — helps determine O vs non-O blood type. Non-O blood types (A, B, AB) have ~1.5x higher risk of venous thrombosis and higher von Willebrand factor levels.|—|—
rs80357906|BRCA2|CANCER|BRCA2 pathogenic variant — if present, significantly increases lifetime risk of breast (~45-70%) and ovarian (~10-20%) cancer, plus elevated prostate and pancreatic cancer risk. Warrants genetic counseling.|—|—
rs28897696|BRCA2|CANCER|BRCA2 pathogenic variant — another known deleterious BRCA2 mutation. Same implications as above.|—|—
rs80357713|BRCA1|CANCER|BRCA1 pathogenic variant — if present, lifetime breast cancer risk ~55-70%, ovarian ~40-45%. Warrants immediate genetic counseling and likely increased screening.|—|—
rs16942|BRCA1|CANCER|BRCA1 common polymorphism (K1183R) — very common, generally benign. Used as a tag SNP for BRCA1 haplotype analysis. Does NOT indicate BRCA1 mutation on its own.|T|C
rs1042522|TP53|CANCER|p53 Pro72Arg — the Arg (G) allele is more efficient at inducing apoptosis, while Pro (C) is better at DNA repair/cell cycle arrest. Context-dependent cancer associations; neither is clearly bad.|C|G
rs61764370|KRAS|CANCER|KRAS-variant — the T allele disrupts a let-7 microRNA binding site, associated with increased risk of non-small cell lung cancer (especially in never-smokers) and ovarian cancer.|T|G
rs1800566|NQO1|CANCER|NQO1 Pro187Ser — TT genotype has near-zero NQO1 enzyme activity, reducing ability to detoxify quinones (benzene metabolites, some chemo drugs). Increases risk of therapy-related leukemia.|T|C
rs4244285|CYP2C19|PHARMACOGENOMICS|CYP2C19*2 — the most common loss-of-function allele. Carriers are poor metabolizers of clopidogrel (Plavix — reduced efficacy, higher cardiac event risk), some SSRIs, and PPIs. FDA black box warning on clopidogrel for poor metabolizers.|A|G
rs4986893|CYP2C19|PHARMACOGENOMICS|CYP2C19*3 — another loss-of-function allele, more common in East Asian populations. Same drug implications as *2.|A|G
rs12248560|CYP2C19|PHARMACOGENOMICS|CYP2C19*17 — gain-of-function allele creating ultra-rapid metabolism. May need higher doses of PPIs and some SSRIs; increased clopidogrel activation (actually beneficial for Plavix).|T|C
rs1065852|CYP2D6|PHARMACOGENOMICS|CYP2D6*10 — reduced-function allele. Affects metabolism of codeine (less morphine conversion), tamoxifen (reduced activation), many antidepressants, beta-blockers, and antipsychotics.|A|G
rs3892097|CYP2D6|PHARMACOGENOMICS|CYP2D6*4 — the most common non-functional CYP2D6 allele in Europeans (~20% carrier rate). Affects ~25% of all prescribed drugs. Poor metabolizers get stronger effects from codeine alternatives, SSRIs.|A|G
rs1799853|CYP2C9|PHARMACOGENOMICS|CYP2C9*2 — reduced-function allele. Carriers need lower warfarin doses (bleeding risk at standard doses). Also affects NSAID metabolism and some oral hypoglycemics.|T|C
rs1057910|CYP2C9|PHARMACOGENOMICS|CYP2C9*3 — more severely reduced function than *2. Strongest impact on warfarin dosing. *2/*3 or *3/*3 carriers may need 50-75% dose reductions.|C|A
rs9923231|VKORC1|PHARMACOGENOMICS|VKORC1 -1639G>A — the single biggest genetic determinant of warfarin dose. AA genotype needs ~50% less warfarin than GG. Combined with CYP2C9 status, explains ~40% of warfarin dose variability.|A|G
rs1045642|ABCB1|PHARMACOGENOMICS|MDR1/P-glycoprotein C3435T — affects drug absorption and brain penetration. TT genotype has lower P-gp expression, meaning higher bioavailability of many drugs including digoxin, HIV protease inhibitors, and some statins.|T|C
rs776746|CYP3A5|PHARMACOGENOMICS|CYP3A5*3 — the most common non-expressing allele. Non-expressers (most Europeans) metabolize tacrolimus slowly, needing lower immunosuppressant doses post-transplant.|T|C
rs28399504|CYP2C19|PHARMACOGENOMICS|CYP2C19*4 — rare loss-of-function variant. Same clinical impact as *2 and *3 when present.|A|G
rs2187668|HLA-DQ2.5|IMMUNE|HLA-DQ2.5 — the strongest genetic risk factor for celiac disease. ~95% of celiac patients carry DQ2.5 or DQ8. However, ~30% of the general population carries it too, so it is necessary but not sufficient. Negative result nearly rules out celiac.|T|C
rs7454108|HLA-DQ8|IMMUNE|HLA-DQ8 tag — second major celiac risk haplotype, also a risk factor for type 1 diabetes. Relevant if DQ2.5 is negative but celiac is still suspected.|C|T
rs2395029|HLA-B*5701|IMMUNE|HLA-B*5701 tag — pharmacogenomic screening marker. Carriers must NEVER take abacavir (HIV drug) due to severe hypersensitivity reaction risk. Also associated with autoimmune conditions.|G|T
rs3135388|HLA-DRB1*1501|IMMUNE|HLA-DRB1*1501 tag — the strongest genetic risk factor for multiple sclerosis. AA genotype associated with ~3x increased risk. Also relevant to vitamin D metabolism in MS context.|A|G
rs2476601|PTPN22|IMMUNE|PTPN22 R620W — a master variant for autoimmunity. The A allele increases risk for type 1 diabetes, rheumatoid arthritis, lupus, Hashimoto thyroiditis, and Graves disease. One of the strongest non-HLA autoimmune variants.|A|G
rs12722489|IL2RA|IMMUNE|IL-2 receptor alpha — risk variant for multiple sclerosis and type 1 diabetes. Affects regulatory T cell function and immune tolerance.|A|G
rs11209026|IL23R|IMMUNE|IL-23 receptor R381Q — the A allele is PROTECTIVE against Crohn disease, psoriasis, and ankylosing spondylitis. One of the clearest examples of a protective rare variant.|A|G
rs2104286|IL2RA|IMMUNE|IL-2 receptor alpha (second signal) — autoimmune disease susceptibility, particularly MS and T1D. Affects immune regulation pathways.|A|G
rs4988235|LCT/MCM6|METABOLISM|Lactase persistence — the A allele allows continued lactose digestion in adulthood (dominant). GG = lactose intolerant (the ancestral state). One of the strongest signals of recent human evolution.|A|G
rs182549|LCT/MCM6|METABOLISM|Lactase persistence European tag — T allele = lactase persistent. CC = likely lactose intolerant. Highly correlated with rs4988235 in European populations.|T|C
rs1801394|MTRR|METABOLISM|MTRR A66G — methionine synthase reductase. GG genotype reduces B12 recycling efficiency, increasing homocysteine. Compounds with MTHFR variants. Consider monitoring B12 and homocysteine levels.|G|A
rs1805087|MTR|METABOLISM|MTR A2756G — methionine synthase. GG genotype may increase homocysteine by reducing methionine synthase activity. Part of the folate/B12/methylation cycle.|G|A
rs602662|FUT2|METABOLISM|FUT2 secretor status — affects whether you secrete blood type antigens in saliva and gut mucus. Non-secretors have different gut microbiome composition and lower B12 absorption but higher norovirus resistance.|A|G
rs601338|FUT2|METABOLISM|FUT2 nonsecretor (W143X) — AA genotype = nonsecretor. ~20% of Europeans. Resistant to most norovirus strains but may need more B12. Affects gut microbiome diversity.|A|G
rs12934922|BCMO1|NUTRITION|Beta-carotene monooxygenase — TT genotype reduces conversion of beta-carotene to vitamin A by ~50%. Vegetarians/vegans with this variant may need preformed vitamin A (retinol) sources.|T|A
rs7501331|BCMO1|NUTRITION|Beta-carotene conversion (second variant) — TT genotype further reduces carotenoid conversion. Combined with rs12934922 can reduce conversion by up to 70%.|T|C
rs2282679|GC|NUTRITION|Vitamin D binding protein (GC/DBP) — CC genotype associated with significantly lower circulating 25(OH)D levels. May need higher vitamin D supplementation to reach adequate levels, especially at northern latitudes.|C|A
rs10741657|CYP2R1|NUTRITION|Vitamin D 25-hydroxylase — AA genotype reduces the liver enzyme that activates vitamin D. Combined with GC variants, can substantially impact vitamin D status.|A|G
rs1544410|VDR|NUTRITION|Vitamin D receptor BsmI — AA genotype may reduce vitamin D signaling efficiency at the cellular level even with adequate blood levels. Relevant to bone density and immune function.|A|G
rs7946|PEMT|NUTRITION|PEMT V175M — TT genotype (especially in premenopausal women) increases dietary choline requirement substantially. Choline is critical for liver function, brain development, and methylation. Found in eggs, liver, and soybeans.|T|C
rs174547|FADS1|NUTRITION|Fatty acid desaturase 1 — CC genotype is more efficient at converting plant omega-3 (ALA) to EPA/DHA. TT genotype may benefit more from direct fish oil/DHA supplementation rather than relying on flaxseed/walnuts.|C|T
rs1800562|HFE|METABOLISM|HFE C282Y — the primary hereditary hemochromatosis mutation. AA = homozygous, high risk of iron overload. AG = carrier, usually asymptomatic but should monitor ferritin. Most common in Northern European descent.|A|G
rs1799945|HFE|METABOLISM|HFE H63D — milder hemochromatosis variant. GG homozygotes have slight iron overload risk. Compound heterozygotes (C282Y + H63D) have moderate risk. Monitor ferritin periodically.|G|C
rs4680|COMT|DETOX|COMT Val158Met — one of the most studied neurogenetic variants. GG (Val/Val) = fast dopamine breakdown (warrior — stress resilient but lower baseline dopamine). AA (Met/Met) = slow breakdown (worrier — higher dopamine, better focus but more stress sensitive). AG = intermediate.|A|G
rs1695|GSTP1|DETOX|GSTP1 Ile105Val — GG (Val/Val) has reduced glutathione conjugation activity, impairing phase II detoxification of carcinogens, heavy metals, and environmental toxins. May benefit from cruciferous vegetables and NAC.|G|A
rs1056836|CYP1B1|DETOX|CYP1B1 Leu432Val — CC (Val/Val) has higher enzyme activity, producing more 4-OH estrogen metabolites (potentially genotoxic). Relevant to estrogen-related cancer risk. DIM and I3C from cruciferous vegetables may help.|C|G
rs762551|CYP1A2|DETOX|CYP1A2*1F — AA genotype = fast caffeine metabolizer (coffee may be cardioprotective). AC/CC = slow metabolizer (caffeine lingers longer, higher heart attack risk with heavy coffee consumption). Affects caffeine, melatonin, and some drug metabolism.|C|A
rs6265|BDNF|MENTAL_HEALTH|BDNF Val66Met — TT (Met/Met) significantly reduces activity-dependent BDNF secretion, affecting neuroplasticity, memory, and recovery from depression. Exercise is especially beneficial for Met carriers as it boosts BDNF through alternative pathways.|T|C
rs4570625|TPH2|MENTAL_HEALTH|Tryptophan hydroxylase 2 — TT genotype reduces serotonin synthesis in the brain. Associated with increased anxiety, emotional reactivity, and amygdala activation. May respond well to tryptophan-rich foods and stress management.|T|G
rs53576|OXTR|MENTAL_HEALTH|Oxytocin receptor — AA genotype associated with lower empathy scores, reduced social sensitivity, but potentially greater resilience to social stress. GG associated with higher empathy but also greater sensitivity to childhood adversity.|A|G
rs1800497|DRD2/ANKK1|MENTAL_HEALTH|Taq1A (DRD2) — AA (A1/A1) genotype has ~30% fewer D2 dopamine receptors. Associated with increased risk of addictive behaviors, reward-seeking, and reduced impulse control. Exercise and mindfulness can help upregulate dopamine signaling.|A|G
rs6311|HTR2A|MENTAL_HEALTH|Serotonin 2A receptor -1438A/G — TT genotype associated with altered serotonin signaling, may affect antidepressant response (especially SSRIs). Some studies link it to higher response to certain antidepressants.|T|C
rs25531|SLC6A4|MENTAL_HEALTH|5-HTTLPR related — this SNP modifies the serotonin transporter long/short allele expression. G allele effectively converts L allele to function like S allele (lower serotonin reuptake). Relevant to SSRI response and stress sensitivity.|G|A
SNPEOF

# ── Parse SNP DB ─────────────────────────────────────────────────────────────
declare -A SNP_GENE SNP_CAT SNP_DESC SNP_RISK SNP_NORMAL
RS_LIST=()

while IFS='|' read -r rs gene cat desc risk normal; do
    [[ -z "$rs" ]] && continue
    RS_LIST+=("$rs")
    SNP_GENE["$rs"]="$gene"
    SNP_CAT["$rs"]="$cat"
    SNP_DESC["$rs"]="$desc"
    SNP_RISK["$rs"]="$risk"
    SNP_NORMAL["$rs"]="$normal"
done <<< "$SNP_DB"

echo "Searching for ${#RS_LIST[@]} variants..."

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

# ── Risk interpretation ──────────────────────────────────────────────────────
interpret_risk() {
    local rs="$1"
    local geno="${FOUND_GENO[$rs]}"
    local risk="${SNP_RISK[$rs]}"
    local normal="${SNP_NORMAL[$rs]}"

    if [[ "$risk" == "—" || -z "$risk" ]]; then
        echo "PRESENT (interpret with context)"
        return
    fi

    local a1="${geno:0:1}"
    local a2="${geno:1:1}"

    if [[ "$a1" == "$risk" && "$a2" == "$risk" ]]; then
        echo "** HOMOZYGOUS RISK ($geno) ** — both copies are risk allele"
    elif [[ "$a1" == "$risk" || "$a2" == "$risk" ]]; then
        echo "HETEROZYGOUS ($geno) — one risk, one normal"
    elif [[ "$a1" == "$normal" && "$a2" == "$normal" ]]; then
        echo "Normal ($geno) — no risk alleles"
    else
        echo "Genotype: $geno (review against risk allele: $risk)"
    fi
}

# ── APOE helper ──────────────────────────────────────────────────────────────
determine_apoe() {
    local g1="${FOUND_GENO[rs429358]}"
    local g2="${FOUND_GENO[rs7412]}"
    if [[ "$g1" == "TT" && "$g2" == "TT" ]]; then echo "e2/e2 — lowest Alzheimer risk"
    elif [[ "$g1" == "TT" && ( "$g2" == "CT" || "$g2" == "TC" ) ]]; then echo "e2/e3 — below-average risk"
    elif [[ "$g1" == "TT" && "$g2" == "CC" ]]; then echo "e3/e3 — average risk (most common)"
    elif [[ ( "$g1" == "CT" || "$g1" == "TC" ) && "$g2" == "CC" ]]; then echo "e3/e4 — increased risk (~3x)"
    elif [[ "$g1" == "CC" && "$g2" == "CC" ]]; then echo "e4/e4 — high risk (~8-12x)"
    elif [[ ( "$g1" == "CT" || "$g1" == "TC" ) && ( "$g2" == "CT" || "$g2" == "TC" ) ]]; then echo "e2/e4 — mixed"
    else echo "rs429358=$g1 rs7412=$g2 (check manually)"
    fi
}

# ── Category arrays ──────────────────────────────────────────────────────────
CATEGORIES=("FAVISM" "MTHFR" "ALZHEIMERS" "CARDIOVASCULAR" "CLOTTING" "CANCER"
            "PHARMACOGENOMICS" "IMMUNE" "METABOLISM" "NUTRITION" "DETOX" "MENTAL_HEALTH")

CAT_LABELS_FULL=(
    "FAVISM / G6PD DEFICIENCY"
    "MTHFR / FOLATE METABOLISM"
    "ALZHEIMER'S / APOE"
    "CARDIOVASCULAR / DIABETES"
    "CLOTTING / THROMBOPHILIA"
    "CANCER RISK VARIANTS"
    "PHARMACOGENOMICS (DRUG RESPONSE)"
    "IMMUNE / AUTOIMMUNE"
    "METABOLISM"
    "NUTRITION / VITAMINS"
    "DETOXIFICATION"
    "MENTAL HEALTH / NEUROTRANSMITTERS"
)

CAT_LABELS_SHORT=(
    "G6PD/FAVISM" "MTHFR/FOLATE" "APOE/ALZHEIMERS" "CARDIOVASCULAR"
    "CLOTTING" "CANCER" "DRUG METABOLISM" "IMMUNE/AUTOIMMUNE"
    "METABOLISM" "NUTRITION" "DETOX" "MENTAL HEALTH"
)

# =============================================================================
#  WRITE REPORT
# =============================================================================

if [[ "$COMPACT" == "1" ]]; then
# ── COMPACT REPORT ───────────────────────────────────────────────────────────
{
    echo "SNP VARIANT REPORT (COMPACT)"
    echo "Generated: $TIMESTAMP  |  Source: $(basename "$RAW_FILE")"
    echo "Variants found: $FOUND_COUNT / ${#RS_LIST[@]}"
    echo "DISCLAIMER: Educational only, not medical advice."
    echo ""

    for i in "${!CATEGORIES[@]}"; do
        cat="${CATEGORIES[$i]}"
        label="${CAT_LABELS_SHORT[$i]}"
        section=""
        has_entries=0

        for rs in "${RS_LIST[@]}"; do
            [[ "${SNP_CAT[$rs]}" != "$cat" ]] && continue
            [[ "${FOUND_FLAG[$rs]}" != "1" ]] && continue
            has_entries=1
            risk_text=$(interpret_risk "$rs")
            section+="  ${SNP_GENE[$rs]} ($rs): ${FOUND_GENO[$rs]} — $risk_text"$'\n'
        done

        if [[ "$has_entries" == "1" ]]; then
            echo "[$label]"
            echo -n "$section"
            echo ""
        fi
    done

    # APOE
    if [[ "${FOUND_FLAG[rs429358]}" == "1" && "${FOUND_FLAG[rs7412]}" == "1" ]]; then
        echo "[APOE TYPE]"
        echo "  rs429358=${FOUND_GENO[rs429358]} + rs7412=${FOUND_GENO[rs7412]} → $(determine_apoe)"
        echo ""
    fi

    # MTHFR
    if [[ "${FOUND_FLAG[rs1801133]}" == "1" && "${FOUND_FLAG[rs1801131]}" == "1" ]]; then
        g677="${FOUND_GENO[rs1801133]}"
        g1298="${FOUND_GENO[rs1801131]}"
        echo "[MTHFR SUMMARY]"
        echo "  C677T (rs1801133): $g677  |  A1298C (rs1801131): $g1298"
        r677=0; r1298=0
        [[ "${g677:0:1}" == "T" || "${g677:1:1}" == "T" ]] && r677=1
        [[ "${g1298:0:1}" == "G" || "${g1298:1:1}" == "G" ]] && r1298=1
        if [[ "$g677" == "TT" ]]; then
            echo "  → Homozygous C677T: ~30% enzyme activity. Methylfolate recommended."
        elif [[ "$r677" == "1" && "$r1298" == "1" ]]; then
            echo "  → Compound heterozygote: reduced activity. Consider methylfolate."
        elif [[ "$r677" == "1" ]]; then
            echo "  → Heterozygous C677T: ~65% activity. Mild impact."
        elif [[ "$r1298" == "1" ]]; then
            echo "  → Heterozygous A1298C only: mild reduction."
        else
            echo "  → Normal MTHFR function."
        fi
        echo ""
    fi

    # G6PD
    echo "[G6PD SUMMARY]"
    g6pd_found=0
    for rs in rs1050828 rs1050829 rs5030868 rs137852328 rs76723693; do
        if [[ "${FOUND_FLAG[$rs]:-0}" == "1" ]]; then
            risk="${SNP_RISK[$rs]}"
            geno="${FOUND_GENO[$rs]}"
            a1="${geno:0:1}"; a2="${geno:1:1}"
            if [[ "$a1" == "$risk" || "$a2" == "$risk" ]]; then
                echo "  ${SNP_GENE[$rs]} ($rs): $geno — VARIANT DETECTED"
                g6pd_found=1
            fi
        fi
    done
    if [[ "$g6pd_found" == "0" ]]; then
        echo "  No G6PD risk alleles detected in tested variants."
    else
        echo "  Avoid: fava beans, naphthalene, primaquine, dapsone, sulfonamides."
    fi
    echo ""

} > "$REPORT"

else
# ── FULL REPORT ──────────────────────────────────────────────────────────────
{
    echo "============================================================================="
    echo "  SNP VARIANT REPORT"
    echo "  Generated: $TIMESTAMP"
    echo "  Source:    $RAW_FILE"
    echo "  Variants found: $FOUND_COUNT / ${#RS_LIST[@]}"
    echo "============================================================================="
    echo ""
    echo "IMPORTANT DISCLAIMER: This report is for INFORMATIONAL and EDUCATIONAL"
    echo "purposes only. It is NOT a medical diagnosis. Genetic variants are only"
    echo "one factor among many. Always consult a genetic counselor or physician"
    echo "for medical decisions."
    echo ""

    for i in "${!CATEGORIES[@]}"; do
        cat="${CATEGORIES[$i]}"
        label="${CAT_LABELS_FULL[$i]}"

        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "  $label"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""

        for rs in "${RS_LIST[@]}"; do
            [[ "${SNP_CAT[$rs]}" != "$cat" ]] && continue

            if [[ "${FOUND_FLAG[$rs]}" == "1" ]]; then
                risk_text=$(interpret_risk "$rs")
                echo "  ┌─ $rs  ·  ${SNP_GENE[$rs]}  ·  Chr${FOUND_CHR[$rs]}:${FOUND_POS[$rs]}"
                echo "  │"
                echo "  │  Genotype: ${FOUND_GENO[$rs]}"
                echo "  │  Status:   $risk_text"
                echo "  │  Risk allele: ${SNP_RISK[$rs]}  ·  Normal allele: ${SNP_NORMAL[$rs]}"
                echo "  │"
                echo "  │  ${SNP_DESC[$rs]}" | fold -s -w 76 | sed '1!s/^/  │  /'
                echo "  └──────────────────────────────────────────────────────────────────────"
            else
                echo "  ┌─ $rs  ·  ${SNP_GENE[$rs]}  ·  NOT ON CHIP (not tested)"
                echo "  │"
                echo "  │  ${SNP_DESC[$rs]}" | fold -s -w 76 | sed '1!s/^/  │  /'
                echo "  └──────────────────────────────────────────────────────────────────────"
            fi
            echo ""
        done
    done

    echo ""
    echo "============================================================================="
    echo "  APOE GENOTYPE QUICK REFERENCE"
    echo "============================================================================="
    echo ""
    echo "  Your APOE type is determined by rs429358 + rs7412 together:"
    echo "    e2/e2 = rs429358 TT + rs7412 TT  → lowest Alzheimer risk"
    echo "    e2/e3 = rs429358 TT + rs7412 CT  → below-average risk"
    echo "    e3/e3 = rs429358 TT + rs7412 CC  → most common, average risk"
    echo "    e3/e4 = rs429358 CT + rs7412 CC  → ~3x increased risk"
    echo "    e4/e4 = rs429358 CC + rs7412 CC  → ~8-12x increased risk"
    echo "    e2/e4 = rs429358 CT + rs7412 CT  → mixed effects"
    if [[ "${FOUND_FLAG[rs429358]}" == "1" && "${FOUND_FLAG[rs7412]}" == "1" ]]; then
        echo ""
        echo "  >>> YOUR DATA: rs429358=${FOUND_GENO[rs429358]}  rs7412=${FOUND_GENO[rs7412]}"
        echo "  >>> $(determine_apoe)"
    fi
    echo ""
    echo "============================================================================="
    echo "  G6PD / FAVISM NOTES"
    echo "============================================================================="
    echo ""
    echo "  G6PD deficiency is X-linked:"
    echo "    Males (XY): One variant allele = deficient"
    echo "    Females (XX): Usually need two variant alleles, but can be"
    echo "                  symptomatic with one due to X-inactivation patterns"
    echo ""
    echo "  If ANY G6PD variant allele is found, discuss with your doctor."
    echo "  Common triggers to avoid: fava beans, mothballs (naphthalene),"
    echo "  primaquine, dapsone, sulfonamides, high-dose vitamin C,"
    echo "  certain antibiotics. Always check medication labels."
    echo ""
    echo "============================================================================="
    echo "  MTHFR NOTES"
    echo "============================================================================="
    echo ""
    echo "  MTHFR affects your ability to convert folic acid to methylfolate."
    echo "  If you have one or more risk alleles:"
    echo "    - Consider methylfolate (5-MTHF) instead of folic acid"
    echo "    - Eat folate-rich foods (leafy greens, legumes)"
    echo "    - Ask your doctor to check homocysteine levels"
    echo "    - B12 (methylcobalamin form) supports the same pathway"
    echo ""
    echo "  Compound heterozygote (C677T + A1298C) can be as impactful as"
    echo "  C677T homozygous. Check both rs1801133 and rs1801131 above."
    echo ""

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
report of medically relevant variants with genotypes and risk interpretations.
Please analyze these results and provide:

1. A plain-English summary of the most significant findings
2. My APOE genotype and what it means for Alzheimer's risk
3. Whether I carry any G6PD (favism) variants and what I should avoid
4. My MTHFR status and whether I should take methylfolate vs folic acid
5. Any pharmacogenomic findings that affect how I metabolize common drugs
6. Key nutritional considerations (vitamin D, B12, folate, iron, choline)
7. Any clotting disorder risk (Factor V Leiden, prothrombin)
8. Celiac / autoimmune risk markers
9. Anything else notable or actionable

IMPORTANT:
- "NOT ON CHIP" = not tested, NOT negative
- "HOMOZYGOUS RISK" = both copies carry risk allele (highest impact)
- "HETEROZYGOUS" = one risk + one normal (moderate impact)
- "Normal" = no risk alleles at that position
- This is SNP chip data, not whole genome sequencing

--- BEGIN SNP REPORT ---

PROMPTEOF
    cat "$REPORT"
    echo ""
    echo "--- END SNP REPORT ---"
} > "$PROMPT_FILE"

echo "LLM prompt written to: $PROMPT_FILE"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  NEXT STEPS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  1. Review $REPORT for your variant data"
echo ""
echo "  2. Feed the prompt to your local AI:"
echo "     cat $PROMPT_FILE | ollama run llama3:70b"
echo "     — or paste into your local web UI"
echo ""
if [[ "$COMPACT" == "0" ]]; then
echo "  TIP: If the report is too large for your model's context window, re-run:"
echo "     $0 --compact $RAW_FILE"
echo ""
fi
echo "  REMINDER: This is educational, not medical advice."
echo ""