#!/bin/bash

# arguments
while getopts 'hi:r:f:o:p:e:sc' flag; do
    case $flag in
        h)
            echo "Pipeline for GWAS QC (absolute paths are required)"
            echo "options:"
            echo "-i, the basename of input bfile"
            echo "-r, the basename of hapmap bfile"
            echo "-f, the filename of hapmap information"
            echo "-o, the directory of output"
            echo "-p, the directory of src"
            echo '-e, Performing population stratification via IBS values, ex: CHD-CHB-JPT (Select: "CEU", "ASW", "MKK", "MEX", "CHD", "CHB", "JPT", "LWK", "TSI", "GIH", "YRI")'
            echo "-s, check sex"
            echo "-c, apply by-cluster on MDS"
            ;;
        i) IN_FILENAME=$OPTARG;;
        r) HM_FILENAME=$OPTARG;;
        f) HM_INFO=$OPTARG;;
        o) OUT_DIR=$OPTARG;;
        p) SRC_DIR=$OPTARG;;
        e) POPULATION=$OPTARG;;
        s) CHECK_SEX='true';;
        c) BY_CLUSTER='true';;
        *) echo "usage: $0 [-i] [-r] [-f] [-o] [-p] [-e] [-s] [-c]"; exit 1;;
    esac
done

IN_BED=$IN_FILENAME.bed
IN_BIM=$IN_FILENAME.bim
IN_FAM=$IN_FILENAME.fam

HM_BED=$HM_FILENAME.bed
HM_BIM=$HM_FILENAME.bim
HM_FAM=$HM_FILENAME.fam

if [ ! -f "$IN_BED" ] || [ ! -f "$IN_BIM" ] || [ ! -f "$IN_FAM" ]; then
    echo "-i missing or designating error"
    exit 1
fi

if [ ! -f "$HM_BED" ] || [ ! -f "$HM_BIM" ] || [ ! -f "$HM_FAM" ]; then
    echo "-r missing or designating error"
    exit 1
fi

if [ ! -f "$HM_INFO" ]; then
    echo "-f missing or designating error"
    exit 1
fi

if [ -z "$OUT_DIR" ]; then
    echo "-o missing"
    exit 1
fi

if [ ! -f "$SRC_DIR/QC_report.py" ]; then
    echo '-p missing or designating error'
    exit 1
fi

if [ -z "$POPULATION" ]; then
    echo "no population is indicated"
    POPULATION='N'
fi

MAF_QC=0.05
HWE_QC=5e-7
GENO_QC=0.05
MIND_QC=0.05
IBD_QC=0.1875
HETER_SD=3
PRUNE_WINDOW=50
PRUNE_STEP=5
PRUNE_THRESHOLD=0.2

IN_BASENAME=$(basename "${IN_FILENAME}")
HM_BASENAME=$(basename "${HM_FILENAME}")

# preparation
mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit
mkdir QualityControl
mkdir PopulationStratification

# "_" is not allowed in the fam file
sed 's/_/-/g' "$IN_FAM" > QualityControl/"$IN_BASENAME".fam
cp "$IN_BIM" QualityControl/"$IN_BASENAME".bim
cp "$IN_BED" QualityControl/"$IN_BASENAME".bed


##########################################################
########################### QC ###########################
##########################################################
cd QualityControl || exit

echo "==========================================================="
printf "Perfoming QC Step 0 -- Remove dpulicated SNPs\n"
echo "==========================================================="
plink1.9 \
    --bfile "$IN_BASENAME" \
    --list-duplicate-vars suppress-first \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC

echo "==========================================================="
printf "Perfoming QC Step 1 -- Initial QC\n"
echo "==========================================================="
plink1.9 \
    --bfile "$IN_BASENAME" \
    --exclude "$IN_BASENAME".QC.dupvar \
    --maf "$MAF_QC" \
    --hwe "$HWE_QC" \
    --geno "$GENO_QC" \
    --mind "$MIND_QC" \
    --write-snplist \
    --make-just-fam \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC

cp "$IN_BASENAME".QC.log "$IN_BASENAME".QC.log.txt

# missingness
plink1.9 \
    --bfile "$IN_BASENAME" \
    --missing \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC

tail -n +2 "$IN_BASENAME".QC.imiss | awk -v MIND_QC="$MIND_QC" '$6 > MIND_QC {print $1,$2}' > fail_miss.txt

# MAF
plink1.9 \
    --bfile "$IN_BASENAME" \
    --freq \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC

# H-W
plink1.9 \
    --bfile "$IN_BASENAME" \
    --hardy \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC


echo "==========================================================="
printf "Perfoming QC Step 2 -- Variant pruning\n"
echo "==========================================================="
plink1.9 \
    --bfile "$IN_BASENAME" \
    --keep "$IN_BASENAME".QC.fam \
    --extract "$IN_BASENAME".QC.snplist \
    --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_THRESHOLD" \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC


echo "==========================================================="
printf "Perfoming QC Step 3 -- Heterozygosity\n"
echo "==========================================================="
plink1.9 \
    --bfile "$IN_BASENAME" \
    --extract "$IN_BASENAME".QC.prune.in \
    --het \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC

# fail heterozygosity
Rscript "$SRC_DIR"/analyze_heterozygosity.R \
    "$IN_BASENAME".QC.het \
    het_analysis.txt \
    fail_het.txt \
    "$HETER_SD"


echo "==========================================================="
printf "Perfoming QC Step 4 -- Checking gender information\n"
echo "==========================================================="
if [ $CHECK_SEX ]; then
    plink1.9 \
        --bfile "$IN_BASENAME" \
        --extract "$IN_BASENAME".QC.snplist \
        --check-sex \
        --allow-extra-chr \
        --out "$IN_BASENAME".QC
    
    grep "PROBLEM" "$IN_BASENAME".QC.sexcheck > fail_sex.txt
else
    echo -n "" > fail_sex.txt
fi


echo "==========================================================="
printf "Perfoming QC Step 5 -- Plots of analysis\n"
echo "==========================================================="
echo "Plot for Heterozygosity and Missingness of Samples"
Rscript "$SRC_DIR"/imiss-vs-het.R \
    "$IN_BASENAME".QC.het \
    "$IN_BASENAME".QC.imiss \
    analysis_het_vs_imiss.jpg \
    "$MIND_QC" \
    "$HETER_SD"

echo "Plot Missingness of SNPs"
Rscript "$SRC_DIR"/locus_Missing_plot.R \
    "$IN_BASENAME".QC.lmiss \
    analysis_lmiss.jpg \
    "$GENO_QC"

echo "Plot for Minor Allele Frequency"
Rscript "$SRC_DIR"/MAF_plot.R \
    "$IN_BASENAME".QC.frq \
    analysis_maf.jpg \
    "$MAF_QC"

echo "Plot for Hardy Weinberg Equilibrium"
Rscript "$SRC_DIR"/hwe_plot.R \
    "$IN_BASENAME".QC.hwe \
    analysis_hwe.jpg \
    "$HWE_QC"


echo "==========================================================="
printf "Perfoming QC Step 6 -- Removing failed samples\n"
echo "==========================================================="
# merge failed samples
awk '{print $1, $2}' fail_het.txt >> fail.txt
awk '{print $1, $2}' fail_sex.txt >> fail.txt

plink1.9 \
    --bfile "$IN_BASENAME" \
    --keep "$IN_BASENAME".QC.fam \
    --remove fail.txt \
    --make-just-fam \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC


echo "==========================================================="
printf "Perfoming QC Step 7 -- Checking identity by descent\n"
echo "==========================================================="
plink1.9 \
    --bfile "$IN_BASENAME" \
    --extract "$IN_BASENAME".QC.prune.in \
    --keep "$IN_BASENAME".QC.fam \
    --rel-cutoff "$IBD_QC" \
    --allow-extra-chr \
    --allow-no-sex \
    --out "$IN_BASENAME".QC

awk '{print $1, $2}' "$IN_BASENAME".QC.fam > temp1.txt
awk '{print $1, $2}' "$IN_BASENAME".QC.rel.id > temp2.txt
sort temp1.txt temp2.txt | uniq -u | sort -u > fail_ibd.txt
rm temp*

plink1.9 \
    --bfile "$IN_BASENAME" \
    --make-bed \
    --keep "$IN_BASENAME".QC.rel.id \
    --extract "$IN_BASENAME".QC.snplist \
    --allow-extra-chr \
    --allow-no-sex \
    --keep-allele-order \
    --chr 1-22 \
    --out "$IN_BASENAME".QC


##########################################################
#################Population Stratification################
##########################################################
cd ../PopulationStratification || exit


echo "==========================================================="
printf "Population Stratification Step 1 -- QC for hapmap\n"
echo "==========================================================="
plink1.9 \
    --bed "$HM_BED" \
    --bim "$HM_BIM" \
    --fam "$HM_FAM" \
    --maf "$MAF_QC" \
    --hwe "$HWE_QC" \
    --geno "$GENO_QC" \
    --mind "$MIND_QC" \
    --write-snplist \
    --make-just-fam \
    --allow-extra-chr \
    --chr 1-22 \
    --out "$HM_BASENAME".PS


echo "==========================================================="
printf "Population Stratification Step 2 -- Merging input and hapmap\n"
echo "==========================================================="
awk '{print $2}' ../QualityControl/"$IN_BASENAME".QC.bim | grep -F -x -f "$HM_BASENAME".PS.snplist > common_rsid.txt

plink1.9 \
    --bfile ../QualityControl/"$IN_BASENAME".QC \
    --extract common_rsid.txt \
    --make-just-bim \
    --recode \
    --allow-no-sex \
    --out "$IN_BASENAME".PS

plink1.9 \
    --bed "$HM_BED" \
    --bim "$HM_BIM" \
    --fam "$HM_FAM" \
    --keep "$HM_BASENAME".PS.fam \
    --extract common_rsid.txt \
    --update-map "$IN_BASENAME".PS.map 4 2 \
    --make-just-bim \
    --out "$HM_BASENAME".PS

# Resolve strand issues
## get difference between input and heatmap
awk '{print $2,$5,$6}' "$IN_BASENAME".PS.bim > "$IN_BASENAME".snp.txt
awk '{print $2,$5,$6}' "$HM_BASENAME".PS.bim > "$HM_BASENAME".snp.txt
sort "$HM_BASENAME".snp.txt "$IN_BASENAME".snp.txt | uniq -u | awk '{print $1}' | sort -u > flip.snp.txt

## flip SNPs for resolving difference
plink1.9 \
    --bfile ../QualityControl/"$IN_BASENAME".QC \
    --extract common_rsid.txt \
    --flip flip.snp.txt \
    --make-just-bim \
    --allow-no-sex \
    --out "$IN_BASENAME".PS

## get unresolved SNPs
awk '{print $2,$5,$6}' "$IN_BASENAME".PS.bim > "$IN_BASENAME".snp.txt
sort "$HM_BASENAME".snp.txt "$IN_BASENAME".snp.txt | uniq -u | awk '{print $1}' | sort -u > exclusion.snp.txt

## exclude unresolved SNPs
plink1.9 \
    --bfile ../QualityControl/"$IN_BASENAME".QC \
    --extract common_rsid.txt \
    --exclude exclusion.snp.txt \
    --flip flip.snp.txt \
    --allow-no-sex \
    --keep-allele-order \
    --make-bed \
    --out "$IN_BASENAME".PS

plink1.9 \
    --bed "$HM_BED" \
    --bim "$HM_BIM" \
    --fam "$HM_FAM" \
    --keep "$HM_BASENAME".PS.fam \
    --extract common_rsid.txt \
    --exclude exclusion.snp.txt \
    --update-map "$IN_BASENAME".PS.map 4 2 \
    --keep-allele-order \
    --make-bed \
    --out "$HM_BASENAME".PS

# Merge input and heatmap
plink1.9 \
    --bfile "$IN_BASENAME".PS \
    --bmerge "$HM_BASENAME".PS \
    --allow-no-sex \
    --keep-allele-order \
    --make-bed \
    --out PS


echo "==========================================================="
printf "Population Stratification Step 3 -- Performing MDS\n"
echo "==========================================================="
plink1.9 \
    --bfile "$IN_BASENAME".PS \
    --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_THRESHOLD" \
    --allow-no-sex \
    --out "$IN_BASENAME".PS

plink1.9 \
    --bfile PS \
    --extract "$IN_BASENAME".PS.prune.in \
    --genome \
    --allow-no-sex \
    --out PS

if [ $BY_CLUSTER ];then
    plink1.9 \
        --bfile PS \
        --read-genome PS.genome \
        --cluster \
        --ppc 1e-3 \
        --mds-plot 3 by-cluster \
        --allow-no-sex \
        --out PS
else
    plink1.9 \
        --bfile PS \
        --read-genome PS.genome \
        --cluster \
        --ppc 1e-3 \
        --mds-plot 3 \
        --allow-no-sex \
        --out PS
fi

Rscript "$SRC_DIR"/PCA.R \
    PS.mds \
    "$IN_BASENAME".PS.fam \
    "$HM_INFO" \
    "$POPULATION" \
    ./

awk '{print$1, $2, $4, $5, $6}' PS.mds > "$OUT_DIR"/"$IN_BASENAME".MDS.txt


echo "==========================================================="
printf "Population Stratification Step 4 -- Removing Failed Samples\n"
echo "==========================================================="
if [ -f fail_PCA_list.txt ];then
    awk 'NR>1{print $1"\t"$2}' fail_PCA_list.txt > fail_PCA_samples.txt

    plink1.9 \
        --bfile ../QualityControl/"$IN_BASENAME".QC \
        --remove fail_PCA_samples.txt \
        --allow-no-sex \
        --keep-allele-order \
        --make-bed \
        --out ../"$IN_BASENAME".QC.PS

else
    touch fail_PCA_samples.txt
    cp ../QualityControl/"$IN_BASENAME".QC.bed ../"$IN_BASENAME".QC.PS.bed
    cp ../QualityControl/"$IN_BASENAME".QC.bim ../"$IN_BASENAME".QC.PS.bim
    cp ../QualityControl/"$IN_BASENAME".QC.fam ../"$IN_BASENAME".QC.PS.fam
    cp ../QualityControl/"$IN_BASENAME".QC.log ../"$IN_BASENAME".QC.PS.log
fi

cd ../

# build covariates: only consider SEX and PCAs
if [ $CHECK_SEX ]; then
    awk 'FNR==NR{a[$1 FS $2]=$3 FS $4 FS $5;next}{print $1 FS $2 FS $5 FS a[$1 FS $2]}' "$IN_BASENAME".MDS.txt "$IN_BASENAME".QC.PS.fam > "$IN_BASENAME".COV.txt
    echo -e "FID IID SEX C1 C2 C3" | cat - "$IN_BASENAME".COV.txt > tmp && mv tmp "$IN_BASENAME".COV.txt
else
    awk 'FNR==NR{a[$1 FS $2]=$3 FS $4 FS $5;next}{print $1 FS $2 FS a[$1 FS $2]}' "$IN_BASENAME".MDS.txt "$IN_BASENAME".QC.PS.fam > "$IN_BASENAME".COV.txt
    echo -e "FID IID C1 C2 C3" | cat - "$IN_BASENAME".COV.txt > tmp && mv tmp "$IN_BASENAME".COV.txt
fi


##########################################################
##########################Report##########################
##########################################################
python3 "$SRC_DIR/QC_report.py" \
    --first_log_file QualityControl/"$IN_BASENAME".QC.log.txt \
    --last_log_file "$IN_BASENAME".QC.PS.log \
    --fail_sex_file QualityControl/fail_sex.txt \
    --fail_het_file QualityControl/fail_het.txt \
    --fail_miss_file QualityControl/fail_miss.txt \
    --fail_ibd_file QualityControl/fail_ibd.txt \
    --fail_ps_file PopulationStratification/fail_PCA_samples.txt \
    --fig_het_vs_imiss QualityControl/analysis_het_vs_imiss.jpg \
    --fig_hwe QualityControl/analysis_hwe.jpg \
    --fig_maf QualityControl/analysis_maf.jpg \
    --fig_lmiss QualityControl/analysis_lmiss.jpg \
    --fig_population PopulationStratification/population_stratification_C1_C2.jpg \
    --heter "$HETER_SD" \
    --ind_miss "$MIND_QC" \
    --geno_miss "$GENO_QC" \
    --ibd "$IBD_QC" \
    --maf "$MAF_QC" \
    --hwe "$HWE_QC" \
    --prune_window "$PRUNE_WINDOW" \
    --prune_step "$PRUNE_STEP" \
    --prune_threshold "$PRUNE_THRESHOLD" \
    --population "$POPULATION" \
    --logo_file "$SRC_DIR/logo.png" \
    --work_dir ./ \
    --basename "$IN_BASENAME"
