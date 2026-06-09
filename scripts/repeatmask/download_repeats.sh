#!/bin/bash
# Download RepeatMasker .out.gz files from UCSC GenArk hubs
# for all GCF assemblies, convert to BED, and report masked fractions.
# GCA-only and DNA Zoo assemblies are listed separately at the end.

set -euo pipefail

OUTDIR="repeatmasker_beds"
mkdir -p "$OUTDIR"

# Function to convert GCF accession to UCSC hub path
# e.g. GCF_024166365.1 -> hubs/GCF/024/166/365/GCF_024166365.1
gcf_to_path() {
    local acc="$1"
    local prefix="${acc:0:3}"          # GCF
    local nums="${acc:4:9}"            # 024166365
    local p1="${nums:0:3}"             # 024
    local p2="${nums:3:3}"             # 166
    local p3="${nums:6:3}"             # 365
    echo "hubs/${prefix}/${p1}/${p2}/${p3}/${acc}"
}

download_and_convert() {
    local species="$1"
    local gcf="$2"

    local hub_path
    hub_path=$(gcf_to_path "$gcf")
    local url="https://hgdownload.soe.ucsc.edu/${hub_path}/${gcf}.repeatMasker.out.gz"
    local outfile="${OUTDIR}/${species}.repeatMasker.out.gz"
    local bedfile="${OUTDIR}/${species}.repeatMasker.bed"

    echo "=== $species ($gcf) ==="

    if [[ -f "$bedfile" ]]; then
        echo "  BED already exists, skipping."
        return
    fi

    echo "  Downloading: $url"
    if ! wget -q -O "$outfile" "$url"; then
        echo "  ERROR: Download failed for $species - file may not exist on UCSC"
        rm -f "$outfile"
        return
    fi

    echo "  Converting to BED..."
    gunzip -c "$outfile" | \
        tail -n +4 | \
        awk 'BEGIN{OFS="\t"} NF>=7 {print $5, $6, $7}' | \
        sort -k1,1 -k2,2n > "$bedfile"

    local masked
    masked=$(awk '{sum += $3-$2} END {print sum}' "$bedfile")
    local intervals
    intervals=$(wc -l < "$bedfile")
    echo "  Masked bases: $masked"
    echo "  Intervals: $intervals"

    # Clean up the .out.gz to save space (BED is what we need)
    rm -f "$outfile"
    echo "  Done -> $bedfile"
    echo ""
}

# All GCF assemblies (using GCF accession from FTP URLs, not the GCA from TSV)
download_and_convert "Loxodonta_africana"    "GCF_030014295.1"
download_and_convert "Panthera_tigris"       "GCF_018350195.1"
download_and_convert "Diceros_bicornis"      "GCF_020826845.1"
download_and_convert "Panthera_leo"          "GCF_018350215.1"
download_and_convert "Panthera_pardus"       "GCF_024362965.1"
download_and_convert "Panthera_uncia"        "GCF_023721935.1"
download_and_convert "Ceratotherium_simum"   "GCF_000283155.1"

echo "=============================="
echo "Done with GCF assemblies."
echo ""
echo "The following assemblies have NO UCSC RepeatMasker track"
echo "and will need RepeatMasker or WindowMasker run manually:"
echo "  Naja_naja              GCA_009733165.1  (no GCF equivalent)"
echo "  Rhinoceros_unicornis   GCA_028646465.1  (no GCF equivalent)"
echo "  Bos_gaurus             GCA_014182915.2  (no GCF equivalent)"
echo "  Tragelaphus_eurycerus  barney_pseudo2.1 (DNA Zoo assembly)"
