#!/bin/bash

set -e

# --- Configuration (!!! USER MUST REVIEW/EDIT THESE PATHS !!!) ---
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." &> /dev/null && pwd)

GWAS_PROC_DIR="${PROJECT_ROOT}/pipeline_output_no_clumping/processed_gwas"
ADSP_DIR="/path/to/your/adsp/plink/files" # !!! EDIT THIS
OUTPUT_DIR="${PROJECT_ROOT}/pipeline_output_no_clumping"
LOG_DIR="${OUTPUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/pipeline.log"

PYTHON_EXEC="python3" # !!! EDIT IF NEEDED
PLINK_EXEC="/path/to/your/plink" # !!! EDIT THIS
LIFTOVER_EXEC="/path/to/your/liftOver" # !!! EDIT THIS
CHAIN_FILE="${PROJECT_ROOT}/resources/hg19ToHg38.over.chain.gz" # !!! EDIT THIS (Ensure file exists)

# --- Function Definitions ---
log_message() {
    mkdir -p "$(dirname "$LOG_FILE")"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

check_success() {
    local status=$?
    if [ ${status} -ne 0 ]; then
        log_message "Error: Previous command failed with status ${status}. Exiting."
        exit ${status}
    fi
}

# --- Setup ---
mkdir -p "${OUTPUT_DIR}"
check_success
mkdir -p "${LOG_DIR}"
check_success
echo "Starting pipeline script (No Clumping Version) at $(date)" > "${LOG_FILE}"
check_success

# --- Step 0: Update ADSP BIM Files ---
log_message "Step 0: Checking and updating ADSP .bim files in ${ADSP_DIR}..."
BIM_BACKUP_DIR="${ADSP_DIR}/original_bims_backup"
UPDATE_BIM_NEEDED=false
FIRST_BIM=$(find "${ADSP_DIR}" -maxdepth 1 -name 'chr1*_QCnew*.bim' -print -quit)

if [ -n "$FIRST_BIM" ] && [ -f "$FIRST_BIM" ] ; then
    if awk 'BEGIN{FS="[ \t]+"} $2 == "." { found=1; exit } END{ exit !found }' "${FIRST_BIM}"; then
        log_message "Found '.' variant IDs. Update will proceed."
        UPDATE_BIM_NEEDED=true
    else
        log_message "No '.' variant IDs found. Assuming okay or updated."
    fi
else
    log_message "Warning: Cannot find sample chr1 .bim file in ${ADSP_DIR}. Skipping update check."
fi

if [ "$UPDATE_BIM_NEEDED" = true ]; then
    if [ ! -d "${BIM_BACKUP_DIR}" ]; then
        log_message "Creating backup in ${BIM_BACKUP_DIR}..."
        mkdir -p "${BIM_BACKUP_DIR}"
        check_success
        cp ${ADSP_DIR}/chr*_QCnew*.bim "${BIM_BACKUP_DIR}/"
        check_success
        log_message "Backup created."
    else
        log_message "Backup directory ${BIM_BACKUP_DIR} exists."
    fi

    log_message "Updating '.' variant IDs to CHR:POS format..."
    ORIG_DIR=$(pwd)
    cd "${ADSP_DIR}" || exit 1

    BIM_UPDATE_FAILED=false
    shopt -s nullglob
    for bim_file in chr*_QCnew*.bim; do
        if [ -r "${bim_file}" ]; then
            log_message "Processing ${bim_file}..."
            tmp_file="${bim_file}.tmp_update"
            awk 'BEGIN{FS="[ \t]+"; OFS="\t"} { if ($2 == ".") $2 = $1 ":" $4; print $0 }' "${bim_file}" > "${tmp_file}"
            if [ $? -eq 0 ]; then
                if [ -s "${tmp_file}" ]; then
                    mv "${tmp_file}" "${bim_file}"
                else
                    log_message "Error: Empty temp file for ${bim_file}. Skipping."
                    rm -f "${tmp_file}"
                    BIM_UPDATE_FAILED=true
                fi
            else
                log_message "Error: awk failed for ${bim_file}. Skipping."
                rm -f "${tmp_file}"
                BIM_UPDATE_FAILED=true
            fi
        else
             log_message "Warning: Cannot read ${bim_file}. Skipping."
             BIM_UPDATE_FAILED=true
       fi
    done
    shopt -u nullglob
    cd "${ORIG_DIR}" || log_message "Warning: cd back to ${ORIG_DIR} failed"

    if [ "$BIM_UPDATE_FAILED" = true ]; then
        log_message "Warning: Errors during .bim update."
    else
        log_message "Finished updating .bim files."
    fi
else
     log_message "Skipping .bim file update."
fi

# --- Step 1: Process GWAS Summary Statistics ---
log_message "Step 1: Running 02_process_gwas.py..."
PROCESS_GWAS_LOG="${LOG_DIR}/process_gwas.log"
mkdir -p "$(dirname "$GWAS_PROC_DIR")"
mkdir -p "$GWAS_PROC_DIR"
> "${PROCESS_GWAS_LOG}"
${PYTHON_EXEC} "${SCRIPT_DIR}/02_process_gwas.py" >> "${PROCESS_GWAS_LOG}" 2>&1
check_success
log_message "Step 1: GWAS preprocessing finished."

AD_SNPS_FILE="${GWAS_PROC_DIR}/ad_snps.txt"
SCZ_SNPS_FILE="${GWAS_PROC_DIR}/scz_snps.txt"

if [ ! -s "$AD_SNPS_FILE" ]; then log_message "Error: AD SNP file '${AD_SNPS_FILE}' missing or empty."; exit 1; fi
if [ ! -s "$SCZ_SNPS_FILE" ]; then log_message "Error: SCZ SNP file '${SCZ_SNPS_FILE}' missing or empty."; exit 1; fi

# --- Step 2a: Merge ADSP PLINK Files ---
log_message "Step 2a: Merging ADSP PLINK files (GRCh38)..."
MERGE_LIST_FILE="${OUTPUT_DIR}/adsp_merge_list.txt"
ADSP_MERGED_PREFIX="${OUTPUT_DIR}/adsp_nhw_merged"
PLINK_MERGE_LOG="${LOG_DIR}/plink_merge.log"
> "${PLINK_MERGE_LOG}"

if [ -f "${ADSP_MERGED_PREFIX}.bed" ] && [ -f "${ADSP_MERGED_PREFIX}.bim" ] && [ -f "${ADSP_MERGED_PREFIX}.fam" ]; then
    log_message "Merged ADSP file ${ADSP_MERGED_PREFIX} exists. Skipping."
else
    find "${ADSP_DIR}" -maxdepth 1 -name 'chr*_QCnew*.bed' -printf '%h/%f\n' | sed 's/\.bed$//' | sort > "${MERGE_LIST_FILE}"
    check_success
    if [ ! -s "${MERGE_LIST_FILE}" ]; then log_message "Error: Cannot find ADSP PLINK files in ${ADSP_DIR}."; exit 1; fi
    NUM_FILES_TO_MERGE=$(wc -l < "${MERGE_LIST_FILE}")
    log_message "Merge list created: ${NUM_FILES_TO_MERGE} files."
    if [ "${NUM_FILES_TO_MERGE}" -lt 22 ] || [ "${NUM_FILES_TO_MERGE}" -gt 23 ]; then
        log_message "Warning: Found ${NUM_FILES_TO_MERGE} files. Check ${MERGE_LIST_FILE}."
    fi

    ${PLINK_EXEC} --merge-list "${MERGE_LIST_FILE}" --allow-no-sex --make-bed --out "${ADSP_MERGED_PREFIX}" >> "${PLINK_MERGE_LOG}" 2>&1
    check_success
    log_message "Step 2a: ADSP PLINK files merged to ${ADSP_MERGED_PREFIX}."
fi

# --- Step 2b & 2c: Clumping Skipped ---
log_message "Step 2b & 2c: SKIPPING LD Clumping."

# --- Step 2d: Prepare BED files for Liftover ---
log_message "Step 2d: Preparing BED files for Liftover..."
AD_SIG_SNPS_HG19_BED="${OUTPUT_DIR}/ad_significant_snps_hg19.bed"
SCZ_SIG_SNPS_HG19_BED="${OUTPUT_DIR}/scz_significant_snps_hg19.bed"
BED_CREATION_LOG="${LOG_DIR}/bed_creation.log"
> "${AD_SIG_SNPS_HG19_BED}"
> "${SCZ_SIG_SNPS_HG19_BED}"
> "${BED_CREATION_LOG}"

log_message "Creating AD BED..."
awk 'NR>1 { split($1, c, ":"); chr=c[1]; pos=c[2]; if (chr ~ /^[0-9]+$/ && pos ~ /^[0-9]+$/) {OFS="\t"; print "chr"chr, pos-1, pos, chr":"pos} else {print "Warning: Skipping invalid AD line: " $0 > "/dev/stderr"} }' "${AD_SNPS_FILE}" > "${AD_SIG_SNPS_HG19_BED}" 2>> "${BED_CREATION_LOG}"
check_success

log_message "Creating SCZ BED..."
awk 'NR>1 { split($1, c, ":"); chr=c[1]; pos=c[2]; if (chr ~ /^[0-9]+$/ && pos ~ /^[0-9]+$/) {OFS="\t"; print "chr"chr, pos-1, pos, chr":"pos} else {print "Warning: Skipping invalid SCZ line: " $0 > "/dev/stderr"} }' "${SCZ_SNPS_FILE}" > "${SCZ_SIG_SNPS_HG19_BED}" 2>> "${BED_CREATION_LOG}"
check_success
log_message "Step 2d: BED files prepared."

# --- Step 2e: Perform Liftover ---
log_message "Step 2e: Performing Liftover..."
AD_SIG_SNPS_HG38_BED="${OUTPUT_DIR}/ad_significant_snps_hg38.bed"
SCZ_SIG_SNPS_HG38_BED="${OUTPUT_DIR}/scz_significant_snps_hg38.bed"
UNMAPPED_AD_SIG_FILE="${LOG_DIR}/unmapped_ad_significant_liftover.txt"
UNMAPPED_SCZ_SIG_FILE="${LOG_DIR}/unmapped_scz_significant_liftover.txt"
> "${UNMAPPED_AD_SIG_FILE}" > "${UNMAPPED_SCZ_SIG_FILE}" > "${AD_SIG_SNPS_HG38_BED}" > "${SCZ_SIG_SNPS_HG38_BED}"

if [ ! -f "${CHAIN_FILE}" ]; then log_message "Error: Liftover chain file not found: ${CHAIN_FILE}"; exit 1; fi
if [ ! -x "${LIFTOVER_EXEC}" ]; then log_message "Error: Liftover executable not found/executable: ${LIFTOVER_EXEC}"; exit 1; fi
if [ ! -s "${AD_SIG_SNPS_HG19_BED}" ]; then log_message "Warning: Input AD BED is empty."; fi
if [ ! -s "${SCZ_SIG_SNPS_HG19_BED}" ]; then log_message "Warning: Input SCZ BED is empty."; fi

${LIFTOVER_EXEC} "${AD_SIG_SNPS_HG19_BED}" "${CHAIN_FILE}" "${AD_SIG_SNPS_HG38_BED}" "${UNMAPPED_AD_SIG_FILE}"
check_success
${LIFTOVER_EXEC} "${SCZ_SIG_SNPS_HG19_BED}" "${CHAIN_FILE}" "${SCZ_SIG_SNPS_HG38_BED}" "${UNMAPPED_SCZ_SIG_FILE}"
check_success
log_message "Step 2e: Liftover finished."

# --- Step 2f: Extract Final SNP Lists ---
log_message "Step 2f: Extracting final GRCh38 SNP lists..."
AD_SIG_SNPS_HG38_LIST="${OUTPUT_DIR}/ad_significant_snps_hg38.txt"
SCZ_SIG_SNPS_HG38_LIST="${OUTPUT_DIR}/scz_significant_snps_hg38.txt"
MATCH_LOG="${LOG_DIR}/liftover_id_matching.log"
> "${MATCH_LOG}" > "${AD_SIG_SNPS_HG38_LIST}" > "${SCZ_SIG_SNPS_HG38_LIST}"

AWK_SCRIPT='BEGIN { FS = "[ \t]+"; OFS = "\n" } NR == FNR { coord_to_id[$1 ":" $4] = $2; next; } { chr_col = $1; sub(/^chr/, "", chr_col); key = chr_col ":" $3; if (key in coord_to_id) { print coord_to_id[key]; } else { print "Warning: Coordinate " key " from " FILENAME " not found in BIM." > "/dev/stderr"; }}'

if [ ! -s "${AD_SIG_SNPS_HG38_BED}" ]; then
    log_message "Warning: Lifted AD BED empty. Skipping ID extraction."
else
    log_message "Matching AD IDs..."
    awk "$AWK_SCRIPT" "${ADSP_MERGED_PREFIX}.bim" "${AD_SIG_SNPS_HG38_BED}" > "${AD_SIG_SNPS_HG38_LIST}" 2>> "${MATCH_LOG}"
    check_success
fi

if [ ! -s "${SCZ_SIG_SNPS_HG38_BED}" ]; then
     log_message "Warning: Lifted SCZ BED empty. Skipping ID extraction."
else
    log_message "Matching SCZ IDs..."
    awk "$AWK_SCRIPT" "${ADSP_MERGED_PREFIX}.bim" "${SCZ_SIG_SNPS_HG38_BED}" > "${SCZ_SIG_SNPS_HG38_LIST}" 2>> "${MATCH_LOG}"
    check_success
fi
log_message "Step 2f: Finished extracting SNP lists."

# --- Step 3a: Create Combined SNP List ---
log_message "Step 3a: Creating combined GRCh38 SNP list..."
COMBINED_SIG_SNPS_HG38_LIST="${OUTPUT_DIR}/combined_significant_snps_hg38.txt"
cat "${AD_SIG_SNPS_HG38_LIST}" "${SCZ_SIG_SNPS_HG38_LIST}" | sort | uniq > "${COMBINED_SIG_SNPS_HG38_LIST}"
check_success
log_message "Step 3a: Combined list created."

# --- Step 3b: Extract Final SNP Sets from ADSP Data ---
log_message "Step 3b: Extracting final SNP sets using PLINK..."
AD_SIG_OUT_PREFIX="${OUTPUT_DIR}/adsp_ad_significant_snps_hg38"
SCZ_SIG_OUT_PREFIX="${OUTPUT_DIR}/adsp_scz_significant_snps_hg38"
COMBINED_SIG_OUT_PREFIX="${OUTPUT_DIR}/adsp_combined_significant_snps_hg38"

PLINK_EXTRACT_AD_LOG="${LOG_DIR}/plink_extract_ad_significant.log"
PLINK_EXTRACT_SCZ_LOG="${LOG_DIR}/plink_extract_scz_significant.log"
PLINK_EXTRACT_COMBINED_LOG="${LOG_DIR}/plink_extract_combined_significant.log"
> "${PLINK_EXTRACT_AD_LOG}" > "${PLINK_EXTRACT_SCZ_LOG}" > "${PLINK_EXTRACT_COMBINED_LOG}"

run_plink_extract() {
    local list_file="$1"; local out_prefix="$2"; local log_file="$3"; local name="$4"
    if [ -s "${list_file}" ]; then
        log_message "Extracting ${name} SNPs..."
        ${PLINK_EXEC} --bfile "${ADSP_MERGED_PREFIX}" --extract "${list_file}" --make-bed --out "${out_prefix}" >> "${log_file}" 2>&1
        check_success
        ${PLINK_EXEC} --bfile "${out_prefix}" --recode A --out "${out_prefix}_raw" >> "${log_file}" 2>&1
        check_success
    else
        log_message "Skipping ${name} SNP extraction: list empty."
        touch "${out_prefix}.bed" "${out_prefix}.bim" "${out_prefix}.fam" "${out_prefix}_raw.raw"
    fi
}

run_plink_extract "${AD_SIG_SNPS_HG38_LIST}" "${AD_SIG_OUT_PREFIX}" "${PLINK_EXTRACT_AD_LOG}" "AD"
run_plink_extract "${SCZ_SIG_SNPS_HG38_LIST}" "${SCZ_SIG_OUT_PREFIX}" "${PLINK_EXTRACT_SCZ_LOG}" "SCZ"
run_plink_extract "${COMBINED_SIG_SNPS_HG38_LIST}" "${COMBINED_SIG_OUT_PREFIX}" "${PLINK_EXTRACT_COMBINED_LOG}" "Combined"

log_message "Step 3b: Final PLINK filesets created:"
log_message "  AD: ${AD_SIG_OUT_PREFIX}.*"
log_message "  SCZ: ${SCZ_SIG_OUT_PREFIX}.*"
log_message "  Combined: ${COMBINED_SIG_OUT_PREFIX}.*"

log_message "--- Pipeline Finished Successfully (Clumping SKIPPED) at $(date) ---"

exit 0
