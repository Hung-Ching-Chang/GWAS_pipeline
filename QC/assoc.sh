#!/bin/bash

# arguments
while getopts 'hi:c:d:o:p:' flag; do
    case $flag in
        h)
            echo "Pipeline for Summary statistics (absolute paths are required)"
            echo "options:"
            echo "-i, the basename of input bfile"
            echo "-c, the filename of covariate"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output bfile"
            echo "-p, the directory of src"
            ;;
        i) BFILE=$OPTARG;;
        c) COVFILE=$OPTARG;;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        p) SRC_DIR=$OPTARG;;
        *) echo "usage: $0 [-i] [-c] [-d] [-o] [-p]"; exit 1;;
    esac
done

if [ ! -f "$BFILE.bed" ] || [ ! -f "$BFILE.bim" ] || [ ! -f "$BFILE.fam" ]; then
    echo "-i missing or designating error"
    exit 1
fi

if [ -z "$WORK_DIR" ]; then
    echo "-d missing"
    exit 1
fi

if [ -z "$BASENAME" ]; then
    echo "-o missing"
    exit 1
fi

if [ ! -f "$SRC_DIR/add_BETA_or_OR.py" ]; then
    echo '-p directory dose not contain add_BETA_or_OR.py'
    exit 1
fi

mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || exit

THRESHOLD="1e-8"

# association
if [ ! -f "$COVFILE" ]; then
    echo "association without covariants"
    plink2 \
        --bfile "${BFILE}" \
        --glm \
        --pfilter 1 \
        --out "${WORK_DIR}/${BASENAME}"
else
    plink2 \
        --bfile "${BFILE}" \
        --glm hide-covar \
        --pfilter 1 \
        --covar "${COVFILE}" \
        --out "${WORK_DIR}/${BASENAME}"
fi


# add BETA or OR
SS_FILE=${WORK_DIR}/${BASENAME}.PHENO1.glm.*
python3 "${SRC_DIR}/add_BETA_or_OR.py" ${SS_FILE}

# extract significant SNPs
awk -v threshold=$THRESHOLD '{if ($12 < threshold) print $0}' ${SS_FILE} > "${WORK_DIR}/${BASENAME}.significant.txt"

# Manhattan plot
manhattan_generator \
    --twopoint ${SS_FILE} \
    --col-chr "#CHROM" \
    --col-name "ID" \
    --col-pos "POS" \
    --col-pvalue "P" \
    --bp \
    --use-pvalues \
    --abline 7.3 \
    --significant-threshold 7.3 \
    --no-annotation \
    --significant-point-size 2 \
    --point-size 1 \
    --graph-title "${BASENAME}" \
    --chr-text-size 10 \
    --exclude-chr 23,24 \
    --output "${WORK_DIR}/${BASENAME}.manhanttan"
