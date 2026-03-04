#!/bin/bash
# Pre-run check for micro-C general pipeline
# Verifies inputs, references, tools, and conda environments before running Snakemake

set -euo pipefail

PASS=0
FAIL=0
WARN=0

pass() { echo "  [PASS] $1"; PASS=$((PASS + 1)); }
fail() { echo "  [FAIL] $1"; FAIL=$((FAIL + 1)); }
warn() { echo "  [WARN] $1"; WARN=$((WARN + 1)); }

# ---- 1. fastqList.txt ----
echo "== Checking fastqList.txt =="
if [[ ! -f fastqList.txt ]]; then
    fail "fastqList.txt not found"
else
    pass "fastqList.txt exists"
    # Check format: 3 tab-separated columns
    bad_lines=$(awk -F'\t' 'NF != 3 && NF > 0' fastqList.txt | wc -l | tr -d ' ')
    if [[ "$bad_lines" -gt 0 ]]; then
        fail "fastqList.txt has $bad_lines line(s) without exactly 3 tab-separated columns"
    else
        pass "fastqList.txt format OK (3 columns)"
    fi
    # Check each sample's fastq files exist
    while IFS=$'\t' read -r name r1 r2; do
        [[ -z "$name" ]] && continue
        if [[ ! -f "$r1" ]]; then
            fail "Sample '$name' R1 not found: $r1"
        else
            pass "Sample '$name' R1 exists"
        fi
        if [[ ! -f "$r2" ]]; then
            fail "Sample '$name' R2 not found: $r2"
        else
            pass "Sample '$name' R2 exists"
        fi
    done < fastqList.txt
fi

# ---- 2. Reference files ----
echo ""
echo "== Checking reference files =="

ref_fasta="/tgen_labs/barthel/references/CHM13v2/chm13v2.0.fasta"
ref_genome="/tgen_labs/barthel/references/CHM13v2/chm13v2.0.size.genome"

if [[ ! -f "$ref_fasta" ]]; then
    fail "Reference FASTA not found: $ref_fasta"
else
    pass "Reference FASTA exists"
fi

# BWA index files
for ext in amb ann bwt pac sa; do
    if [[ ! -f "${ref_fasta}.${ext}" ]]; then
        fail "BWA index missing: ${ref_fasta}.${ext}"
    else
        pass "BWA index .${ext} exists"
    fi
done

if [[ ! -f "$ref_genome" ]]; then
    fail "Genome size file not found: $ref_genome"
else
    pass "Genome size file exists"
fi

# ---- 3. External scripts/tools ----
echo ""
echo "== Checking external scripts =="

srcget_qc="/tgen_labs/barthel/software/github/external/Micro-C/get_qc.py"
srcJuicer="/tgen_labs/barthel/software/github/external/Micro-C/juicer_tools_1.22.01.jar"
karyotype="/tgen_labs/barthel/software/miniforge3/envs/micro-C/data/karyotype/karyotype.human.chm13v2.txt"

if [[ ! -f "$srcget_qc" ]]; then
    fail "get_qc.py not found: $srcget_qc"
else
    pass "get_qc.py exists"
fi

if [[ ! -f "$srcJuicer" ]]; then
    fail "juicer_tools jar not found: $srcJuicer"
else
    pass "juicer_tools jar exists"
fi

if [[ ! -f "$karyotype" ]]; then
    fail "Circos karyotype file not found: $karyotype"
else
    pass "Circos karyotype file exists"
fi

# Local scripts
for script in scripts/fixMAforCircos.py scripts/makeCircosConf.sh; do
    if [[ ! -f "$script" ]]; then
        fail "Local script not found: $script"
    else
        pass "$script exists"
    fi
done

if [[ ! -f "config/circos.template.conf" ]]; then
    fail "Circos template config not found: config/circos.template.conf"
else
    pass "Circos template config exists"
fi

# ---- 4. Software in PATH ----
echo ""
echo "== Checking software in PATH =="

for tool in bwa samtools snakemake java python; do
    if command -v "$tool" &>/dev/null; then
        pass "$tool found: $(command -v $tool)"
    else
        fail "$tool not found in PATH"
    fi
done

# ---- 5. Conda environments ----
echo ""
echo "== Checking conda environments =="

if command -v conda &>/dev/null; then
    pass "conda found"
    for env in micro-C telomereC.py3.1; do
        if conda env list 2>/dev/null | grep -qw "$env"; then
            pass "Conda env '$env' exists"
        else
            fail "Conda env '$env' not found"
        fi
    done
else
    warn "conda not found — cannot verify conda environments"
fi

# Tools expected inside conda envs (check if available via conda run)
echo ""
echo "== Checking tools in conda environments =="

# micro-C env: pairtools, cooler, pairix, bgzip, circos
for tool in pairtools cooler pairix bgzip circos; do
    if conda run -n micro-C which "$tool" &>/dev/null; then
        pass "$tool found in micro-C env"
    else
        fail "$tool not found in micro-C env"
    fi
done

# telomereC.py3.1 env: bamCoverage
if conda run -n telomereC.py3.1 which bamCoverage &>/dev/null; then
    pass "bamCoverage found in telomereC.py3.1 env"
else
    fail "bamCoverage not found in telomereC.py3.1 env"
fi

# ---- 6. Snakefile config ----
echo ""
echo "== Checking Snakefile =="

if grep -q '^include: "workflow/pipeline.smk"' Snakefile; then
    pass "General pipeline (pipeline.smk) is active"
else
    warn "General pipeline (pipeline.smk) is NOT active in Snakefile"
fi

if grep -q '^include: "workflow/pipeline2.smk"' Snakefile; then
    warn "Telomere-C pipeline (pipeline2.smk) is also active — is this intended?"
fi

# ---- Summary ----
echo ""
echo "=============================="
echo "  PASS: $PASS  |  FAIL: $FAIL  |  WARN: $WARN"
echo "=============================="

if [[ "$FAIL" -gt 0 ]]; then
    echo "Fix the above FAILs before running the pipeline."
    exit 1
else
    echo "All checks passed. Ready to run."
    exit 0
fi
