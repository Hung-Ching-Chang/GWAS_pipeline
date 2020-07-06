#!/bin/bash
set -e

# OPTIONS
if [ -z $IN_CASE ] && [ -z $IN_CONTROL ] ; then
    echo "Err: missing input parameters"
    exit 1
fi

if [ -z $REF_PANEL ]; then
    echo "Err: missing reference panel in M3VCF format. \
    Recommmend downloading from 'https://genome.sph.umich.edu/wiki/Minimac3' or using the resource in PUB"
    exit 2
fi

if [ -z $OUT_DIR ] ; then
    echo "Err: missing output parameters"
    exit 3
fi

# input dir and reference dir
if [[ $IN_CASE =~ ([^[:space:]]+)/ ]]; then IN_DIR=${BASH_REMATCH[1]}; fi
if [[ $REF_PANEL =~ ([^[:space:]]+)/ ]]; then REF_DIR=${BASH_REMATCH[1]}; fi
mkdir -p $OUT_DIR
cd $OUT_DIR
mkdir -p WorkingDir
cd WorkingDir

echo "Preparing input file ..."
# CASE
if [ -s "$IN_CASE" ]
then
    Filename="$(ls $IN_CASE | awk -F/ '{print $NF}' | awk -F'.' '{print $1}')"
    mkdir -p CASE_GROUP
    unzip $IN_CASE -d $OUT_DIR/WorkingDir/CASE_GROUP
fi

# CONTROL
if [ -s "$IN_CONTROL" ]
then
    Filename="$(ls $IN_CONTROL | awk -F/ '{print $NF}' | awk -F'.' '{print $1}')"
    mkdir -p CONTROL_GROUP
    unzip $IN_CONTROL -d $OUT_DIR/WorkingDir/CONTROL_GROUP
fi

# prepare ref panel 'm3vcf'
echo "Preparing reference panel ..."
mkdir -p Reference
cp $REF_PANEL $OUT_DIR/WorkingDir/Reference
tar -xzf ./Reference/*.tar.gz -C ./Reference/

# =======================================================================================================
# ========================================= Minimac Imputation ==========================================
# =======================================================================================================
echo "============================"
echo "Using Minimac4 for Imputation"
echo "============================"
echo ""
echo ""
echo "Processing CASE group ..."
if [ -s "$IN_CASE" ]
then
    cd $OUT_DIR/WorkingDir/CASE_GROUP
    # Impute the Autosomes (Chr1-22)
    # Set Chromosome Start and End Parameters
    for chr in {1..22}; do
        echo "Processing Chromosome ${chr} Scripts "
        echo "======================================"
            
        #Search the reference directory for the chromosome specific reference map, legend, and hap files and create their respective variables on the fly
        echo "Looking for Minimac .m3vcf Reference Files ... "
        echo "Found the following references for Chromosome ${chr}: "

        MiniRef="$(ls ../Reference/ | egrep --ignore-case "^${chr}[^[:digit:]]{1}.*\.m3vcf\.gz")"
        printf "    $MiniRef \n"

        # unzip haps
        gunzip ${Filename}_Chr${chr}_Phased.haps.gz
        echo "========================================"
        echo "Convert Chr${chr} Phased .haps to .vcf"
        echo "========================================"
        printf '.\n.\n.\n'
        shapeit -convert \
                --input-haps ${Filename}_Chr${chr}_Phased \
                --output-vcf ${Filename}_Chr${chr}_Phased.vcf.gz \
                --output-log ${Filename}_Chr${chr}_Phased.vcf.log
        printf 'done'
        printf '\n\n'
        # Minimac Command	
        echo "========================================"
        echo "Impute Chr${chr} Using Minimac"
        echo "========================================"
        printf '.\n.\n.\n'
        # save output files here before zip
        mkdir -p OUTPUT
        minimac4 --cpus $REQUEST_CPU \
                 --allTypedSites \
                 --minRatio 0.00001 \
                 --ignoreDuplicates \
                 --refHaps ../Reference/${MiniRef} \
                 --haps ${Filename}_Chr${chr}_Phased.vcf.gz \
                 --prefix OUTPUT/${Filename}_Chr${chr}_Imputed
        printf 'done'
        printf '\n\n'	
    done

    # ============================================================================================
    # ====================== Minimac4 Post Imputation File Cleaning ==============================
    # ============================================================================================
    echo
    echo "================================================================"
    echo "Performing Post Imputation File Clean Up on Minimac4 Files"
    echo "================================================================"
    # Make Directory in which to place merged concatenated files
    #------------------------------------------------------------
    echo
    echo "Creating Concat Cohort Folder ..."
    echo
    mkdir -p ConcatImputation/

    for chr in {1..22}; do	
        # Use Plink to Filter .dosage.vcf.gz with in-built R2 to .dosage.vcf.gz
        # =========================================================================================
        # =========================================================================================	
        echo "Using Plink to Filter the Chr${chr}.dosage.vcf.gz to an R2 filtered .vcf.gz" 
        echo "========================================================================================="
        
        # perform the filtration to make .vcf.gz file
        # Runs Plink to filter the concatenated dosage VCF to an R2 Filtered VCF

        ### Specify an IMPUTE INFO or Minimac4 R2 (they are more or less the same) score threshold. This will be used to filter out poorly imputed variants (OPTIONAL)
        ### Set to "0" to keep everything/ default is 0.3
        plink2 --vcf OUTPUT/${Filename}_Chr${chr}_Imputed.dose.vcf.gz dosage=DS \
               --exclude-if-info "R2<=${IMPUTE_INFO_THRESHOLD}" \
               --export vcf vcf-dosage=GP \
               --threads $REQUEST_CPU \
               --out ConcatImputation/${Filename}_Chr${chr}
            
        # Zip the Resulting File
        echo "Zipping the Plink Filtered .VCF.gz file for Chromosome ${chr}"
        echo "============================================================"
        gzip -f ConcatImputation/${Filename}_Chr${chr}.vcf
    done

    # ======================================================================================================
    #                  Use BCFTools to Merge the VCF Files Created by Plink
    # ======================================================================================================	
    # If there are .vcf.gz files to concatenate then list them (in order) on the screen
    echo
    echo "Concatenating the following VCF.gz files using BCFTools:"
    echo "============================================================"
    ls -1av ConcatImputation/*.vcf.gz
    echo "============================================================"
    echo

    # Set List entries as a variable
    VCF2Merge="$(find ConcatImputation/ -maxdepth 1 -type f -name "*.vcf.gz" |sort -V)"

    # Use BCFTools to Merge the VCF Files Listed in the variable
    bcftools concat \
             --threads ${REQUEST_CPU} ${VCF2Merge} \
             --output-type z \
             --output ConcatImputation/Imputed_${Filename}_Merged.vcf.gz

    # Change Permission so the cat script can access it
    chmod -f 700 ConcatImputation/Imputed_${Filename}_Merged.vcf.gz || true 

    echo "==========================================================="
    printf "VCF to bfile\n"
    echo "==========================================================="
    # Modify fam file and save as phen file (case)
    #cat FAIL_INDS.FAIL_CASE.txt case.fam | sort -u -k1,1 | awk '{ if ($3!="") print  $1, $2, $6}' | sed 's/FAM//g' > Case.phen
    # vcf to bfile 
    plink1.9 --vcf ConcatImputation/Imputed_${Filename}_Merged.vcf.gz \
                --make-bed \
                --out ${Filename}.Imputed_CASE
    
    plink1.9 --allow-no-sex \
             --bfile ${Filename}.Imputed_CASE \
             --geno $CONTROL_GENO_QC \
             --maf $CONTROL_MAF_QC \
             --hwe $CASE_HWE_QC \
             --make-bed \
             --out ${Filename}_QC.Imputed_CASE

    cat ${Filename}_QC.Imputed_CASE.fam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"2}' > ${Filename}_QC_phen.Imputed_CASE.fam
    mv ${Filename}_QC.Imputed_CASE.bed $OUT_BFILE_CASE_BED
    mv ${Filename}_QC.Imputed_CASE.bim $OUT_BFILE_CASE_BIM
    mv ${Filename}_QC_phen.Imputed_CASE.fam $OUT_BFILE_CASE_FAM
fi
echo "CASE group Done ..."

echo ""
echo ""
echo "Processing CONTROL group ..."
if [ -s "$IN_CONTROL" ]
then
    cd $OUT_DIR/WorkingDir/CONTROL_GROUP
    # Impute the Autosomes (Chr1-22)
    # Set Chromosome Start and End Parameters
    for chr in {1..22}; do
        echo "Processing Chromosome ${chr} Scripts "
        echo "======================================"
            
        #Search the reference directory for the chromosome specific reference map, legend, and hap files and create their respective variables on the fly
        echo "Looking for Minimac .m3vcf Reference Files ... "
        echo "Found the following references for Chromosome ${chr}: "

        MiniRef="$(ls ../Reference/ | egrep --ignore-case "^${chr}[^[:digit:]]{1}.*\.m3vcf\.gz")"
        printf "    $MiniRef \n"

        # unzip haps
        gunzip ${Filename}_Chr${chr}_Phased.haps.gz
        echo "========================================"
        echo "Convert Chr${chr} Phased .haps to .vcf"
        echo "========================================"
        printf '.\n.\n.\n'
        shapeit -convert \
                --input-haps ${Filename}_Chr${chr}_Phased \
                --output-vcf ${Filename}_Chr${chr}_Phased.vcf.gz \
                --output-log ${Filename}_Chr${chr}_Phased.vcf.log
        printf 'done'
        printf '\n\n'
        # Minimac Command	
        echo "========================================"
        echo "Impute Chr${chr} Using Minimac"
        echo "========================================"
        printf '.\n.\n.\n'
        # save output files here before zip
        mkdir -p OUTPUT
        minimac4 --cpus $REQUEST_CPU \
                 --allTypedSites \
                 --minRatio 0.00001 \
                 --ignoreDuplicates \
                 --refHaps ../Reference/${MiniRef} \
                 --haps ${Filename}_Chr${chr}_Phased.vcf.gz \
                 --prefix OUTPUT/${Filename}_Chr${chr}_Imputed
        printf 'done'
        printf '\n\n'	
    done

    # ============================================================================================
    # ====================== Minimac4 Post Imputation File Cleaning ==============================
    # ============================================================================================
    echo
    echo "================================================================"
    echo "Performing Post Imputation File Clean Up on Minimac4 Files"
    echo "================================================================"
    # Make Directory in which to place merged concatenated files
    #------------------------------------------------------------
    echo
    echo "Creating Concat Cohort Folder ..."
    echo
    mkdir -p ConcatImputation/

    for chr in {1..22}; do	
        # Use Plink to Filter .dosage.vcf.gz with in-built R2 to .dosage.vcf.gz
        # =========================================================================================
        # =========================================================================================	
        echo "Using Plink to Filter the Chr${chr}.dosage.vcf.gz to an R2 filtered .vcf.gz" 
        echo "========================================================================================="
        
        # perform the filtration to make .vcf.gz file
        # Runs Plink to filter the concatenated dosage VCF to an R2 Filtered VCF

        ### Specify an IMPUTE INFO or Minimac4 R2 (they are more or less the same) score threshold. This will be used to filter out poorly imputed variants (OPTIONAL)
        ### Set to "0" to keep everything/ default is 0.3
        plink2 --vcf OUTPUT/${Filename}_Chr${chr}_Imputed.dose.vcf.gz dosage=DS \
               --exclude-if-info "R2<=${IMPUTE_INFO_THRESHOLD}" \
               --export vcf vcf-dosage=GP \
               --threads $REQUEST_CPU \
               --out ConcatImputation/${Filename}_Chr${chr}
            
        # Zip the Resulting File
        echo "Zipping the Plink Filtered .VCF.gz file for Chromosome ${chr}"
        echo "============================================================"
        gzip -f ConcatImputation/${Filename}_Chr${chr}.vcf
    done

    # ======================================================================================================
    #                  Use BCFTools to Merge the VCF Files Created by Plink
    # ======================================================================================================	
    # If there are .vcf.gz files to concatenate then list them (in order) on the screen
    echo
    echo "Concatenating the following VCF.gz files using BCFTools:"
    echo "============================================================"
    ls -1av ConcatImputation/*.vcf.gz
    echo "============================================================"
    echo

    # Set List entries as a variable
    VCF2Merge="$(find ConcatImputation/ -maxdepth 1 -type f -name "*.vcf.gz" |sort -V)"

    # Use BCFTools to Merge the VCF Files Listed in the variable
    bcftools concat \
             --threads ${REQUEST_CPU} ${VCF2Merge} \
             --output-type z \
             --output ConcatImputation/Imputed_${Filename}_Merged.vcf.gz

    # Change Permission so the cat script can access it
    chmod -f 700 ConcatImputation/Imputed_${Filename}_Merged.vcf.gz || true 

    echo "==========================================================="
    printf "VCF to bfile\n"
    echo "==========================================================="
    # Modify fam file and save as phen file (CONTROL)
    #cat FAIL_INDS.FAIL_CONTROL.txt control.fam | sort -u -k1,1 | awk '{ if ($3!="") print  $1, $2, $6}' | sed 's/FAM//g' > Control.phen
    # vcf to bfile 
    plink1.9 --vcf ConcatImputation/Imputed_${Filename}_Merged.vcf.gz \
             --make-bed \
             --out ${Filename}.Imputed_CONTROL

    plink1.9 --allow-no-sex \
             --bfile ${Filename}.Imputed_CONTROL \
             --geno $CONTROL_GENO_QC \
             --maf $CONTROL_MAF_QC \
             --hwe $CONTROL_HWE_QC \
             --make-bed \
             --out ${Filename}_QC.Imputed_CONTROL

    cat ${Filename}_QC.Imputed_CONTROL.fam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"2}' > ${Filename}_QC_phen.Imputed_CONTROL.fam
    mv ${Filename}_QC.Imputed_CONTROL.bed $OUT_BFILE_CONTROL_BED
    mv ${Filename}_QC.Imputed_CONTROL.bim $OUT_BFILE_CONTROL_BIM
    mv ${Filename}_QC_phen.Imputed_CONTROL.fam $OUT_BFILE_CONTROL_FAM
fi
echo "CONTROL group Done ..."

cd $OUT_DIR
rm -r WorkingDir
