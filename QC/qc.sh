#!/bin/bash

IN_BED="/volume/rickchang/RCVS_analysis/TWB_only/TWBR10901-05_TWB1.bed"
IN_BIM="/volume/rickchang/RCVS_analysis/TWB_only/TWBR10901-05_TWB1.bim"
IN_FAM="/volume/rickchang/RCVS_analysis/TWB_only/TWBR10901-05_TWB1.fam"

HM_BED="/volume/prsdata/GWAS_REF/hapmap3_r1_b37_fwd_consensus.hapmap3.bed"
HM_BIM="/volume/prsdata/GWAS_REF/hapmap3_r1_b37_fwd_consensus.hapmap3.bim"
HM_FAM="/volume/prsdata/GWAS_REF/hapmap3_r1_b37_fwd_consensus.hapmap3.fam"
HM_INFO="/volume/prsdata/GWAS_REF/relationships_w_pops_121708.hapmap3"

OUT_DIR="/volume/rickchang/RCVS_analysis/TWB_only/OUT/"

MAF_QC=0.05
HWE_QC=5e-7
GENO_QC=0.05
MIND_QC=0.05
IBD_QC=0.1875
HETER_SD=3
PRUNE_WINDOW=50
PRUNE_STEP=5
PRUNE_THRESHOLD=0.2

mkdir -p $OUT_DIR
cd $OUT_DIR
mkdir QualityControl
mkdir PopulationStratification

IN_BASENAME=$(basename ${IN_BED%%.*})

# "_" is not allowed in the fam file
sed 's/_/-/g' $IN_FAM > QualityControl/$IN_BASENAME.fam
cp $IN_BIM QualityControl/$IN_BASENAME.bim
cp $IN_BED QualityControl/$IN_BASENAME.bed


##########################################################
########################### QC ###########################
##########################################################
cd QualityControl

echo "==========================================================="
printf "Perfoming QC Step 1 -- Initial QC\n"
echo "==========================================================="
plink1.9 \
    --bfile $IN_BASENAME \
    --maf $MAF_QC \
    --hwe $HWE_QC \
    --geno $GENO_QC \
    --mind $MIND_QC \
    --write-snplist \
    --make-just-fam \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex

cp $IN_BASENAME.QC.log $IN_BASENAME.QC.log.txt

# missingness
plink1.9 \
    --bfile $IN_BASENAME \
    --missing \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex

tail -n +2 $IN_BASENAME.QC.imiss | awk -v MIND_QC=$MIND_QC '$6 > MIND_QC {print $1,$2}' > fail_miss.txt

# MAF
plink1.9 \
    --bfile $IN_BASENAME \
    --freq \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex

# H-W
plink1.9 \
    --bfile $IN_BASENAME \
    --hardy \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex


echo "==========================================================="
printf "Perfoming QC Step 2 -- Variant pruning\n"
echo "==========================================================="
plink1.9 \
    --bfile $IN_BASENAME \
    --keep $IN_BASENAME.QC.fam \
    --extract $IN_BASENAME.QC.snplist \
    --indep-pairwise $PRUNE_WINDOW $PRUNE_STEP $PRUNE_THRESHOLD \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex


echo "==========================================================="
printf "Perfoming QC Step 3 -- Heterozygosity\n"
echo "==========================================================="
plink1.9 \
    --bfile $IN_BASENAME \
    --extract $IN_BASENAME.QC.prune.in \
    --het \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex

# fail heterozygosity
Rscript /opt/analyze_heterozygosity.R \
    $IN_BASENAME.QC.het \
    het_analysis.txt \
    fail_het.txt \
    $HETER_SD


echo "==========================================================="
printf "Perfoming QC Step 4 -- Checking gender information\n"
echo "==========================================================="
plink1.9 \
    --bfile $IN_BASENAME \
    --extract $IN_BASENAME.QC.snplist \
    --check-sex \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex

grep "PROBLEM" $IN_BASENAME.QC.sexcheck > fail_sex.txt


echo "==========================================================="
printf "Perfoming QC Step 5 -- Plots of analysis\n"
echo "==========================================================="
echo "Plot for Heterozygosity and Missingness of Samples"
Rscript /opt/imiss-vs-het.R \
    $IN_BASENAME.QC.het \
    $IN_BASENAME.QC.imiss \
    analysis_het_vs_imiss.jpg \
    $MIND_QC \
    $HETER_SD

echo "Plot Missingness of SNPs"
Rscript /opt/locus_Missing_plot.R \
    $IN_BASENAME.QC.lmiss \
    analysis_lmiss.jpg \
    $GENO_QC

echo "Plot for Minor Allele Frequency"
Rscript /opt/MAF_plot.R \
    $IN_BASENAME.QC.frq \
    analysis_maf.jpg \
    $MAF_QC

echo "Plot for Hardy Weinberg Equilibrium"
Rscript /opt/hwe_plot.R \
    $IN_BASENAME.QC.hwe \
    analysis_hwe.jpg \
    $HWE_QC


echo "==========================================================="
printf "Perfoming QC Step 6 -- Removing failed samples\n"
echo "==========================================================="
# merge failed samples
cat fail_het.txt | awk '{print $1, $2}' >> fail.txt
cat fail_sex.txt | awk '{print $1, $2}' >> fail.txt

plink1.9 \
    --bfile $IN_BASENAME \
    --keep $IN_BASENAME.QC.fam \
    --remove fail.txt \
    --make-just-fam \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex


echo "==========================================================="
printf "Perfoming QC Step 7 -- Checking identity by descent\n"
echo "==========================================================="
plink1.9 \
    --bfile $IN_BASENAME \
    --extract $IN_BASENAME.QC.prune.in \
    --keep $IN_BASENAME.QC.fam \
    --rel-cutoff $IBD_QC \
    --allow-extra-chr \
    --out $IN_BASENAME.QC \
    --allow-no-sex

awk '{print $1, $2}' $IN_BASENAME.QC.fam > temp1.txt
awk '{print $1, $2}' $IN_BASENAME.QC.rel.id > temp2.txt
sort temp1.txt temp2.txt | uniq -u | sort -u > fail_ibd.txt
rm temp*

plink1.9 \
    --bfile $IN_BASENAME \
    --make-bed \
    --keep $IN_BASENAME.QC.rel.id \
    --extract $IN_BASENAME.QC.snplist \
    --allow-extra-chr \
    --chr 1-22 \
    --out $IN_BASENAME.QC \
    --allow-no-sex


##########################################################
#################Population Stratification################
##########################################################
HM_BASENAME=$(basename ${HM_BED%%.*})

cd ../PopulationStratification


echo "==========================================================="
printf "Population Stratification Step 1 -- QC for hapmap\n"
echo "==========================================================="
plink1.9 \
    --bed $HM_BED \
    --bim $HM_BIM \
    --fam $HM_FAM \
    --maf $MAF_QC \
    --hwe $HWE_QC \
    --geno $GENO_QC \
    --mind $MIND_QC \
    --write-snplist \
    --make-just-fam \
    --allow-extra-chr \
    --chr 1-22 \
    --out $HM_BASENAME.PS \
    --allow-no-sex


echo "==========================================================="
printf "Population Stratification Step 2 -- Merging input and hapmap\n"
echo "==========================================================="
awk '{print $2}' ../QualityControl/$IN_BASENAME.QC.bim | grep -F -x -f $HM_BASENAME.PS.snplist > common_rsid.txt

plink1.9 \
    --bfile ../QualityControl/$IN_BASENAME.QC \
    --extract common_rsid.txt \
    --make-just-bim \
    --recode \
    --out $IN_BASENAME.PS \
    --allow-no-sex

plink1.9 \
    --bed $HM_BED \
    --bim $HM_BIM \
    --fam $HM_FAM \
    --keep $HM_BASENAME.PS.fam \
    --extract common_rsid.txt \
    --update-map $IN_BASENAME.PS.map 4 2 \
    --make-just-bim \
    --out $HM_BASENAME.PS \
    --allow-no-sex

# Resolve strand issues
## get difference between input and heatmap
awk '{print $2,$5,$6}' $IN_BASENAME.PS.bim > $IN_BASENAME.snp.txt
awk '{print $2,$5,$6}' $HM_BASENAME.PS.bim > $HM_BASENAME.snp.txt
sort $HM_BASENAME.snp.txt $IN_BASENAME.snp.txt | uniq -u | awk '{print $1}' | sort -u > flip.snp.txt

## flip SNPs for resolving difference
plink1.9 \
    --bfile ../QualityControl/$IN_BASENAME.QC \
    --extract common_rsid.txt \
    --flip flip.snp.txt \
    --make-just-bim \
    --out $IN_BASENAME.PS \
    --allow-no-sex

## get unresolved SNPs
awk '{print $2,$5,$6}' $IN_BASENAME.PS.bim > $IN_BASENAME.snp.txt
sort $HM_BASENAME.snp.txt $IN_BASENAME.snp.txt | uniq -u | awk '{print $1}' | sort -u > exclusion.snp.txt

## exclude unresolved SNPs
plink1.9 \
    --bfile ../QualityControl/$IN_BASENAME.QC \
    --extract common_rsid.txt \
    --exclude exclusion.snp.txt \
    --flip flip.snp.txt \
    --make-bed \
    --out $IN_BASENAME.PS \
    --allow-no-sex

plink1.9 \
    --bed $HM_BED \
    --bim $HM_BIM \
    --fam $HM_FAM \
    --keep $HM_BASENAME.PS.fam \
    --extract common_rsid.txt \
    --exclude exclusion.snp.txt \
    --update-map $IN_BASENAME.PS.map 4 2 \
    --make-bed \
    --out $HM_BASENAME.PS \
    --allow-no-sex

# Merge input and heatmap
plink1.9 \
    --bfile $IN_BASENAME.PS \
    --bmerge $HM_BASENAME.PS \
    --allow-no-sex \
    --make-bed \
    --out PS


echo "==========================================================="
printf "Population Stratification Step 3 -- Performing MDS\n"
echo "==========================================================="
plink1.9 \
    --bfile $IN_BASENAME.PS \
    --indep-pairwise $PRUNE_WINDOW $PRUNE_STEP $PRUNE_THRESHOLD \
    --out $IN_BASENAME.PS

plink1.9 \
    --bfile PS \
    --extract $IN_BASENAME.PS.prune.in \
    --genome \
    --out PS

plink1.9 \
    --bfile PS \
    --read-genome PS.genome \
    --cluster \
    --ppc 1e-3 \
    --mds-plot 3 \
    --allow-no-sex \
    --out PS

Rscript /opt/PCA.R \
    PS.mds \
    $IN_BASENAME.PS.fam \
    $HM_INFO \
    $Population_Stratification \
    ./

awk '{print$1, $2, $4, $5, $6}' PS.mds > $OUT_DIR/covar_mds.txt
  

echo "==========================================================="
printf "Population Stratification Step 4 -- Removing Failed Samples\n"
echo "==========================================================="
if [ -f fail_PCA_list.txt ];then
    awk 'NR>1{print $1"\t"$2}' fail_PCA_list.txt > fail_PCA_samples.txt

    plink1.9 \
        --bfile ../QualityControl/$IN_BASENAME.QC \
        --remove fail_PCA_samples.txt \
        --allow-no-sex \
        --make-bed \
        --out ../$IN_BASENAME.QC.PS

else
    cp ../QualityControl/$IN_BASENAME.QC.bed ../$IN_BASENAME.QC.PS.bed
    cp ../QualityControl/$IN_BASENAME.QC.bim ../$IN_BASENAME.QC.PS.bim
    cp ../QualityControl/$IN_BASENAME.QC.fam ../$IN_BASENAME.QC.PS.fam
    cp ../QualityControl/$IN_BASENAME.QC.log ../$IN_BASENAME.QC.PS.log
fi
