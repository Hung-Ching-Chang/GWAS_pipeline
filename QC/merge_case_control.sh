#!/bin/bash

IN_BFILE_CASE="/volume/rickchang/RCVS_analysis/merge_association/RCVS_only/OUT/chrall_556case_only.QC.PS.bed"
IN_BFILE_CONTROL="/volume/rickchang/RCVS_analysis/merge_association/TWB_only/OUT/TWBR10901-05_TWB1.QC.PS.bed"

OUT_DIR="/volume/rickchang/RCVS_analysis/merge_association/merge_again/"
# input dir
if [[ $IN_BFILE_CASE =~ ([^[:space:]]+)/ ]]; then IN_CASE_DIR=${BASH_REMATCH[1]}; fi
if [[ $IN_BFILE_CONTROL =~ ([^[:space:]]+)/ ]]; then IN_CONTROL_DIR=${BASH_REMATCH[1]}; fi
mkdir -p $OUT_DIR
# Set Working Directory
mkdir -p WorkingDir
cd WorkingDir

IN_CASE=${IN_BFILE_CASE%.*}
IN_CONTROL=${IN_BFILE_CONTROL%.*}
echo "==========================================================="
printf "STEP 1: Preparing case and control data\n"
echo "==========================================================="
printf "Finding duplicate SNPs ID\n"
printf ".\n.\n.\n"
cat $IN_CASE.bim | awk '{print $2}' > all_snps.snplist
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink1.9 --bfile $IN_CASE \
         --exclude duplicated_snps.snplist \
         --make-bed \
         --out unique_CASE \
         --allow-no-sex \
         --allow-extra-chr \
         --chr 1-22 \


cat $IN_CONTROL.bim | awk '{print $2}' > all_snps.snplist
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink1.9 --bfile $IN_CONTROL \
         --exclude duplicated_snps.snplist \
         --make-bed \
         --out unique_CONTROL\
         --allow-no-sex \
         --allow-extra-chr \
         --chr 1-22 \


printf "Finding Positional and Allelic Duplicates\n"
printf ".\n.\n.\n"
plink1.9 --bfile unique_CASE \
         --list-duplicate-vars ids-only suppress-first \
         --out Dups2Remove_CASE \
         --allow-no-sex
plink1.9 --bfile unique_CONTROL \
         --list-duplicate-vars ids-only suppress-first \
         --out Dups2Remove_CONTROL \
         --allow-no-sex

printf "Removing Positional and Allelic Duplicates if they exist\n"
if [[ $(wc -l <Dups2Remove_CASE.dupvar) > 0 ]]
then 
    printf ".\n.\n.\n"
    plink1.9 --bfile unique_CASE \
             --exclude Dups2Remove_CASE.dupvar \
             --make-bed \
             --out CASE_Dups2Remove \
             --allow-no-sex
    printf "Done \n"
    mv CASE_Dups2Remove.bed unique_CASE.bed
    mv CASE_Dups2Remove.bim unique_CASE.bim
    mv CASE_Dups2Remove.fam unique_CASE.fam
fi

if [[ $(wc -l <Dups2Remove_CONTROL.dupvar) > 0 ]]
then 
    printf ".\n.\n.\n"
    plink1.9 --bfile unique_CONTROL \
             --exclude Dups2Remove_CONTROL.dupvar \
             --make-bed \
             --out CONTROL_Dups2Remove \
             --allow-no-sex
    printf "Done \n"
    mv CONTROL_Dups2Remove.bed unique_CONTROL.bed
    mv CONTROL_Dups2Remove.bim unique_CONTROL.bim
    mv CONTROL_Dups2Remove.fam unique_CONTROL.fam
fi


echo "==========================================================="
plink1.9 --bfile unique_CONTROL \
         --write-snplist \
         --allow-no-sex \
         --out CONTROL
awk '{print $2}' unique_CASE.bim | grep -F -x -f CONTROL.snplist > common_rsid.txt

# Resolve strand issues
awk '{print $2,$5,$6}' unique_CONTROL.bim > unique_CONTROL.snp.txt
awk '{print $2,$5,$6}' unique_CASE.bim > unique_CASE.snp.txt
sort unique_CASE.snp.txt unique_CONTROL.snp.txt | uniq -u | awk '{print $1}' | sort -u > flip.snp.txt
## flip SNPs for resolving difference
plink1.9 \
    --bfile unique_CONTROL \
    --extract common_rsid.txt \
    --flip flip.snp.txt \
    --make-just-bim \
    --allow-no-sex \
    --out flip_CONTROL
## get unresolved SNPs
awk '{print $2,$5,$6}' flip_CONTROL.bim > flip_CONTROL.snp.txt
sort unique_CASE.snp.txt flip_CONTROL.snp.txt | uniq -u | awk '{print $1}' | sort -u > exclusion.snp.txt
## exclude unresolved SNPs
plink1.9 \
    --bfile unique_CONTROL \
    --extract common_rsid.txt \
    --exclude exclusion.snp.txt \
    --flip flip.snp.txt \
    --allow-no-sex \
    --make-bed \
    --out flip_CONTROL

plink1.9 \
    --bfile unique_CASE \
    --extract common_rsid.txt \
    --exclude exclusion.snp.txt \
    --make-bed \
    --out common_CASE

# Merge input and heatmap
plink1.9 \
    --bfile common_CASE \
    --bmerge flip_CONTROL \
    --allow-no-sex \
    --make-bed \
    --out CASE_CONTROL_merged



echo "==========================================================="
printf "STEP 5: exclude INDEL (insertion deletion)\n"
echo "==========================================================="
# extract SNPs with INDEL
cat CASE_CONTROL_merged.bim | awk '{if (length($5)>1 || length($6)>1) print $2}' > INDEL.txt

plink1.9 --bfile CASE_CONTROL_merged \
        --exclude INDEL.txt \
        --make-bed \
        --out QC_DONE \
        --allow-no-sex

mv QC_DONE.bed $OUT_DIR/RCVS_TWB_QC_merged.bed
mv QC_DONE.bim $OUT_DIR/RCVS_TWB_QC_merged.bim
mv QC_DONE.fam $OUT_DIR/RCVS_TWB_QC_merged.fam
