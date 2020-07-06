#!/bin/bash
set -e

# OPTIONS
if [ -z $IN_BFILE_ONE ] ; then
    echo "Err: missing input parameters"
    exit 1
fi

if [ -z $REF_PANEL ]; then
    echo "Err: missing reference panel"
    exit 2
fi

if [ -z $OUT_DIR ] ; then
    echo "Err: missing output parameters"
    exit 3
fi

# input dir and reference dir
if [[ $REF_PANEL =~ ([^[:space:]]+)/ ]]; then RES_DIR=${BASH_REMATCH[1]}; fi
if [[ $IN_BFILE_ONE =~ ([^[:space:]]+)/ ]]; then INPUT_DIR=${BASH_REMATCH[1]}; fi
mkdir -p $OUT_DIR
cd $OUT_DIR
if [ -s "$IN_BFILE_ONE" ] 
then
    BfileName="$(ls $IN_BFILE_ONE | awk -F/ '{print $NF}' | awk -F'.' '{print $1}')"
else
    BfileName="$(ls $IN_BFILE_TWO | awk -F/ '{print $NF}' | awk -F'.' '{print $1}')"
fi

mkdir -p WorkingDir
cd WorkingDir

echo "Preparing input dataset ..."
Extention_one="${IN_BFILE_ONE#*.}"
if [ -s "$IN_BFILE_TWO" ]
then
    if [ -s "$IN_BFILE_ONE" ] 
    then
        echo "Detect Case and Control!"
        mkdir Case_group
        mkdir Control_group
        cp $INPUT_DIR/$BfileName.CASE.* ./Case_group
        cp $INPUT_DIR/$BfileName.CONTROL.* ./Control_group
    else
        echo "Detect Control only!"
        # input extension: .CONTROL.bfile
        mkdir Control_group
        cp $INPUT_DIR/$BfileName.CONTROL.* ./Control_group
    fi
else
    if [ $Extention_one == CASE.bed ] 
    then
        echo "Detect Case only!"
        # input extension: .CASE.bfile
        mkdir Case_group
        cp $INPUT_DIR/$BfileName.CASE.* ./Case_group
    elif [ $Extention_one == CONTROL.bed ]
    then
        echo "Detect Control only!"
        # input extension: .CONTROL.bfile
        mkdir Control_group
        cp $INPUT_DIR/$BfileName.CONTROL.* ./Control_group
    else
        # splite control group
        cat $INPUT_DIR/$BfileName.fam | awk '{if ($6 == 2) print $1"\t"$2}' > case_list.txt
        if [ -s case_list.txt ]
        then
            echo "Detect Case!"
            mkdir Case_group
            plink1.9 --allow-no-sex \
                     --bfile $INPUT_DIR/$BfileName \
                     --keep case_list.txt \
                     --make-bed \
                     --out ./Case_group/${BfileName}.CASE
        fi
        # input extension: .bfile
        # splite case group
        cat $INPUT_DIR/$BfileName.fam | awk '{if ($6 == 1) print $1"\t"$2}' > control_list.txt
        if [ -s control_list.txt ]
        then
            echo "Detect Control!"
            mkdir Control_group
            plink1.9 --allow-no-sex \
                     --bfile $INPUT_DIR/$BfileName \
                     --keep control_list.txt \
                     --make-bed \
                     --out ./Control_group/${BfileName}.CONTROL
        fi
    fi
fi

echo "Preparing 1K Genome Phase 3 Ref Panel"
mkdir -p ./Reference
cp -r $RES_DIR/*.tgz ./Reference/
echo "Unpackaging Ref Panel ..."
tar -xzf ./Reference/1000GP_Phase3.tgz -C ./Reference/
echo "Cleaning Up ..."
mv ./Reference/1000GP_Phase3/* ./Reference/
# Delete the now empty directory and the tgz zipped Ref panel
rm -r ./Reference/1000GP_Phase3/

echo "==================================="
echo "Processing CASE group"
echo "==================================="
if [ -d Case_group ]
then
    cd $OUT_DIR/WorkingDir/Case_group
    echo 
    echo "==========================================================="
    printf "check missingness\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
    plink1.9 --bfile $BfileName.CASE \
             --missing \
             --out missing
    printf "Done \n"
    cat missing.lmiss | awk '{if ($5 == 1) print $0}' > all_miss_loc
    cat missing.imiss | awk '{if ($5 == 1) print $0}' > all_miss_ind
    if [[ $(wc -l <all_miss_loc) > 0 ]] || [[ $(wc -l <all_miss_ind) > 0 ]] 
    then
        echo "SHAPEIT doesn't allow a completely missing SNP or individual (call rate = 0%) included in the dataset."
        exit 4
    fi
    echo 
    echo "==========================================================="
    printf "Finding Positional and Allelic Duplicates\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
    plink1.9 --bfile $BfileName.CASE \
             --list-duplicate-vars ids-only suppress-first \
             --out ./Dups2Remove \
             --allow-no-sex
    printf "Done \n"
    echo 
    echo "==========================================================="
    printf "Removing Positional and Allelic Duplicates if they exist\n"
    echo "==========================================================="
    if [[ $(wc -l <Dups2Remove.dupvar) > 0 ]]
    then 
        printf ".\n.\n.\n"
        plink1.9 --bfile $BfileName.CASE \
                 --exclude ./Dups2Remove.dupvar \
                 --make-bed \
                 --out ${BfileName}_Dups2Remove.CASE \
                 --allow-no-sex
        printf "Done \n"
        mv ${BfileName}_Dups2Remove.CASE.bed $BfileName.CASE.bed
        mv ${BfileName}_Dups2Remove.CASE.bim $BfileName.CASE.bim
        mv ${BfileName}_Dups2Remove.CASE.fam $BfileName.CASE.fam
    fi

    echo "==================================="
    echo "split by chromosome"
    echo "==================================="
    for chr in {1..22}; do
        plink1.9 --bfile $BfileName.CASE \
                 --chr ${chr} \
                 --make-bed \
                 --out $BfileName.chr${chr}
    done

    echo "============================"
    echo "Using Shapeit for Phasing"
    echo "============================"
    mkdir -p ./output
    if [[ $(wc -l <$BfileName.CASE.fam) < 100 ]]; then  REQUEST_CPU=1; fi
    echo "REQUEST_CPU setting: $REQUEST_CPU"
    #Set Chromosome Start and End Parameters
    for chr in {1..22}; do
        #Search the reference directory for the chromosome specific reference map, legend, and hap files and create their respective variables on the fly
        echo "Processing Chromosome ${chr} Script ..."
        echo "======================================"
        echo "Pre-Check:  "
        echo "Found the Following Shapeit References for Chromosome ${chr}: "

        GeneticMap="$(ls ../Reference/ | egrep --ignore-case ".*map.*chr${chr}[^[:digit:]]{1}.*|.*chr${chr}[^[:digit:]]{1}.*map.*")"
            echo "   Genetic Map File: $GeneticMap"	

        ## -----------------------------------------------
        # shapeit function 
        ## -----------------------------------------------  
        # the default values of [burn, prune, main] are [7, 8, 20]
        shapeit --input-bed $BfileName.chr${chr} \
                --input-map ../Reference/$GeneticMap \
                --output-max ./output/${BfileName}_Chr${chr}_Phased.haps.gz ./output/${BfileName}_Chr${chr}_Phased.sample \
                --output-log ./output/${BfileName}_Chr${chr}_Phased.log \
                --thread $REQUEST_CPU \
                --burn 7 \
                --prune 8 \
                --main 20 \
                --force || true
        # run forcely when data contain individual with high rates of missing data (>5%)
    done
    cd ./output
    zip -r $BfileName.phased_CASE.zip ./*
    mv $BfileName.phased_CASE.zip ${OUT_BFILE_CASE}
fi
echo ""

echo "==================================="
echo "Processing CONTROL group"
echo "==================================="
if  [ -d Control_group ]
then
    cd $OUT_DIR/WorkingDir/Control_group
    echo 
    echo "==========================================================="
    printf "check missingness\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
    plink1.9 --bfile $BfileName.CONTROL \
             --missing \
             --out missing
    printf "Done \n"
    cat missing.lmiss | awk '{if ($5 == 1) print $0}' > all_miss_loc
    cat missing.imiss | awk '{if ($5 == 1) print $0}' > all_miss_ind
    if [[ $(wc -l <all_miss_loc) > 0 ]] || [[ $(wc -l <all_miss_ind) > 0 ]] 
    then
        echo "SHAPEIT doesn't allow a completely missing SNP or individual (call rate = 0%) included in the dataset."
        exit 5
    fi
    echo 
    echo "==========================================================="
    printf "Finding Positional and Allelic Duplicates\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
    plink1.9 --bfile $BfileName.CONTROL \
             --list-duplicate-vars ids-only suppress-first \
             --out ./Dups2Remove \
             --allow-no-sex
    printf "Done \n"
    echo 
    echo "==========================================================="
    printf "Removing Positional and Allelic Duplicates if they exist\n"
    echo "==========================================================="
    if [[ $(wc -l <Dups2Remove.dupvar) > 0 ]]
    then 
        printf ".\n.\n.\n"
        plink1.9 --bfile $BfileName.CONTROL \
                 --exclude ./Dups2Remove.dupvar \
                 --make-bed \
                 --out ${BfileName}_Dups2Remove.CONTROL \
                 --allow-no-sex
        printf "Done \n"
        mv ${BfileName}_Dups2Remove.CONTROL.bed $BfileName.CONTROL.bed
        mv ${BfileName}_Dups2Remove.CONTROL.bim $BfileName.CONTROL.bim
        mv ${BfileName}_Dups2Remove.CONTROL.fam $BfileName.CONTROL.fam
    fi

    echo "==================================="
    echo "split by chromosome"
    echo "==================================="
    for chr in {1..22}; do
        plink1.9 --bfile $BfileName.CONTROL \
                --chr ${chr} \
                --make-bed \
                --out $BfileName.chr${chr}
    done

    echo "============================"
    echo "Using Shapeit for Phasing"
    echo "============================"
    mkdir -p ./output
    if [[ $(wc -l <$BfileName.CONTROL.fam) < 100 ]]; then  REQUEST_CPU=1; fi
    #Set Chromosome Start and End Parameters
    for chr in {1..22}; do
        #Search the reference directory for the chromosome specific reference map, legend, and hap files and create their respective variables on the fly
        echo "Processing Chromosome ${chr} Script ..."
        echo "======================================"
        echo "Pre-Check:  "
        echo "Found the Following Shapeit References for Chromosome ${chr}: "

        GeneticMap="$(ls ../Reference/ | egrep --ignore-case ".*map.*chr${chr}[^[:digit:]]{1}.*|.*chr${chr}[^[:digit:]]{1}.*map.*")"
            echo "   Genetic Map File: $GeneticMap"	

        ## -----------------------------------------------
        # shapeit function 
        ## -----------------------------------------------  
        # the default values of [burn, prune, main] are [7, 8, 20]
        shapeit --input-bed $BfileName.chr${chr} \
                --input-map ../Reference/$GeneticMap \
                --output-max ./output/${BfileName}_Chr${chr}_Phased.haps.gz ./output/${BfileName}_Chr${chr}_Phased.sample \
                --output-log ./output/${BfileName}_Chr${chr}_Phased.log \
                --thread $REQUEST_CPU \
                --burn 7 \
                --prune 8 \
                --main 20 \
                --force
        # run forcely when data contain individual with high rates of missing data (>5%)
        done
    cd ./output
    zip -r $BfileName.phased_CONTROL.zip ./*
    mv $BfileName.phased_CONTROL.zip ${OUT_BFILE_CONTROL}
fi
# Clean
cd $OUT_DIR
rm -r WorkingDir
echo "==========================================================="
echo "Phasing for Imputation Done!!!"
echo "==========================================================="

