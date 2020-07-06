#!/bin/bash
set -e
free -m


# OPTIONS
if [ -z $IN_BED ] || [ -z $IN_BIM ] || [ -z $IN_FAM ] ; then
    echo "Err: missing input parameters"
    exit 1
fi
if [ $Fix_Input_with_Reference_Annotation == "TRUE" ]; then
    if [ -z $REF_FASTA ] || [ -z $REF_FASTA_INDEX ]|| [ -z $ANNOTATION ] || [ -z $ANNOTATION_INDEX ] ; then
        echo "Err: missing reference files..."
        echo "Checkout to PUB Resource for the Reference Data from 1K Genomes and the Annotation Files from NCBI"
        exit 2
    fi
fi

if [ -z $OUT_CASE_BED ] || [ -z $OUT_CONTROL_BED ] ; then
    echo "Err: missing output parameters"
    exit 3
fi


# ===============================================================================
# ================ Working Directory / split Case and Control ===================
# ===============================================================================
# input dir and reference dir
if [[ $IN_BED =~ ([^[:space:]]+)/ ]]; then IN_DIR=${BASH_REMATCH[1]}; fi
if [[ $REF_FASTA =~ ([^[:space:]]+)/ ]]; then RES_DIR=${BASH_REMATCH[1]}; fi
mkdir -p $OUT_DIR
cd $OUT_DIR

echo "Build up Working Space ..."
# Set Working Directory
mkdir -p WorkingDir
cd WorkingDir
# cp Rscript and Python script
cp -r /opt/* .
# Get the BaseName of the Data Placed in the input folder -- must have a Plink .bim file in the folder
RawData="$(ls $IN_BED | awk -F/ '{print $NF}' | awk -F'.' '{print $1}')"
# do not allow "-" in the original fam file
sed 's/_/-/g' $IN_FAM > $RawData.fam
cp $IN_BIM $RawData.bim
cp $IN_BED $RawData.bed

# splite case group
cat $RawData.fam | awk '{if ($6 == 1) print $1"\t"$2}' > control_list.txt
if [ -s control_list.txt ]
then
    plink1.9 --allow-no-sex \
             --bfile $RawData \
             --keep control_list.txt \
             --make-bed \
             --out ${RawData}.CONTROL
fi

# splite control group
cat $RawData.fam | awk '{if ($6 == 2) print $1"\t"$2}' > case_list.txt
if [ -s case_list.txt ]
then
    plink1.9 --allow-no-sex \
             --bfile $RawData \
             --keep case_list.txt \
             --make-bed \
             --out ${RawData}.CASE
fi

if [ $Fix_Input_with_Reference_Annotation == "TRUE" ]
then
    # Downloading Reference Data and index files from 1K Genomes and NCBI
    # Download the Reference Data ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
    # Download the annotation files from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
    echo "Preparing Reference Data ..."
    mkdir -p RefAnnotationData
    cp $RES_DIR/*.fasta.* RefAnnotationData
    cp $RES_DIR/*.vcf.gz* RefAnnotationData
    cp $RES_DIR/PlinkChrRename.txt RefAnnotationData
    # unzip Reference Data (add "|| true" to ignore meaningless error: decompression OK, trailing garbage ignored)
    echo "unzip Reference Annotation files ..."
    gunzip -d RefAnnotationData/*.fasta.gz || true
    echo "Success!!!"
    # Rezip Fasta File in bgzip
    bgzip ./RefAnnotationData/*.fasta
fi

# ===============================================================================
# =============================== CASE GROUP ====================================
# ===============================================================================
if [ -f ${RawData}.CASE.bed ]
then
    echo "Processing 'CASE' group..."
    mkdir -p CASE
    mv ${RawData}.CASE.* ./CASE
    # Genotype Missingness 
    GenoQC=$CASE_GENO_QC;
    # Individual Missingness 
    MindQC=$CASE_MIND_QC;
    # Hardy Weinberg Equilibrium
    HweQC=$CASE_HWE_QC;
    # Minor Allele Frequency 
    MafQC=$CASE_MAF_QC;
    # IBD
    IBDQC=$CASE_IBD_QC;
    # standard deviation for heterozygosity
    Heter_SD=$CASE_HETER_SD;

    #Make case Directory in which all case files will populate
    cd ./CASE

    # if case only, fix allele based on reference
    if [ $Fix_Input_with_Reference_Annotation == "TRUE" ]
    then
        # --------------------------------------------------------------------------------------------------
        # Convert the Plink BED/BIM/FAM into a VCF into a BCF so that it may be fixed with BCFtools
        # --------------------------------------------------------------------------------------------------
        echo 
        echo "==========================================================="
        printf "Converting Plink files into VCF format \n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            plink1.9 --allow-no-sex \
                     --bfile ${RawData}.CASE \
                     --recode vcf \
                     --out ./DataFixStep1_${RawData}.CASE
        printf "Done \n"

        echo 
        echo "==========================================================="
        printf "Converting VCF into a BCF with chromosome names that match the reference .fasta annotation\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools annotate -Ob --rename-chrs ../RefAnnotationData/PlinkChrRename.txt ./DataFixStep1_${RawData}.CASE.vcf > ./DataFixStep1_${RawData}.CASE.bcf
        printf "Done \n"
        # --------------------------------------------------------------------------------------------------
        # Align Input File to the Reference Annotation (Fix with BCFtools)
        # --------------------------------------------------------------------------------------------------
        echo 
        echo "==========================================================="
        printf "View the number of correctly annotated/aligned variants to the Reference annotation\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep1_${RawData}.CASE.bcf -- -f ../RefAnnotationData/*.fasta.gz 
        printf "Done \n"
        
        echo 
        echo "==========================================================="
        printf "fix the strand if the VCF is TOP-compatible\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep1_${RawData}.CASE.bcf -Ob -o ./DataFixStep2_${RawData}.CASE.bcf -- -f ../RefAnnotationData/*.fasta.gz -m top
        printf "Done \n"
        echo 
        
        echo "==========================================================="
        printf "Flip the alleles based on the downloaded annotation file\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep2_${RawData}.CASE.bcf -Ob -o ./DataFixStep3_${RawData}.CASE.bcf -- -d -f ../RefAnnotationData/*.fasta.gz -i ../RefAnnotationData/*.vcf.gz
        printf "Done \n"
        
        echo 
        echo "==========================================================="
        printf "Sorting the BCF output since fixing it may have made it unsorted\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools sort ./DataFixStep3_${RawData}.CASE.bcf -Ob -o ./DataFixStep4_${RawData}.CASE.bcf
        printf "Done \n"

        echo 
        echo "==========================================================="
        printf "Swap the alleles based on the downloaded annotation file\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep4_${RawData}.CASE.bcf -Ob -o ./DataFixStep5_${RawData}.CASE.bcf -- -f ../RefAnnotationData/*.fasta.gz -m flip -d
        printf "Done \n"


        # Convert BCF back into Plink .bed/.bim/.fam for Shapeit2 Phasing
        echo 
        echo "==========================================================="
        printf "Converting Fixed and Sorted BCF back into Plink .bed/.bim/.fam\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            plink1.9 --bcf ./DataFixStep5_${RawData}.CASE.bcf \
                     --make-bed \
                     --out ./DataFixStep6_${RawData}.CASE
        printf "Done \n"
        cp ${RawData}.CASE.fam DataFixStep6_${RawData}.CASE.fam

        # --------------------------------------------------------------------------------------------------
        # Add back in the sex information
        # --------------------------------------------------------------------------------------------------
        # echo 
        # echo "==========================================================="
        # printf "Restoring Sample Sex Information\n"
        # echo "==========================================================="	
        # printf ".\n.\n.\n"
        #     plink1.9 --bfile ./DataFixStep6_${RawData}.CASE \
        #              --update-sex ${RawData}.CASE.fam 3 \
        #              --allow-no-sex \
        #              --make-bed \
        #              --out ./DataFixStep7_${RawData}.CASE
        # printf "Done \n"
        mv ./DataFixStep6_${RawData}.CASE.bed ${RawData}.CASE.bed
        mv ./DataFixStep6_${RawData}.CASE.bim ${RawData}.CASE.bim
        mv ./DataFixStep6_${RawData}.CASE.fam ${RawData}.CASE.fam
        rm DataFixStep*
    fi 

    # Finally Remove any positional duplicates 
    # i.e. same position and alleles, but differently named variants since Shapeit will not tolerate these
    echo 
    echo "==========================================================="
    printf "Finding Positional and Allelic Duplicates\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
        plink1.9 --bfile ./${RawData}.CASE \
                 --list-duplicate-vars ids-only suppress-first \
                 --out ./Dups2Remove \
                 --allow-no-sex
    printf "Done \n"
    # --------------------------------------------------------------------------------------------------
    # Report Number of duplicates:
    # --------------------------------------------------------------------------------------------------
    echo 
    echo "==========================================================="
    printf "Removing Positional and Allelic Duplicates if they exist\n"
    echo "==========================================================="
        DuplicateNumber="$(wc ./Dups2Remove.dupvar | awk '{print $1}')"
    printf "Found Number of Duplicate Variant(s):\n${DuplicateNumber}"
    printf "\n.\n.\n.\n"
        plink1.9 --bfile ./${RawData}.CASE \
                 --exclude ./Dups2Remove.dupvar \
                 --make-bed \
                 --out ./Dups2Remove_${RawData}.CASE \
                 --allow-no-sex
    printf "Done \n"

    # --------------------------------------------------------------------------------------------------
    # STEP 6: Pre-imputation QC
    # --------------------------------------------------------------------------------------------------
    ##########################################################
    ####################### Sample QC ########################
    ##########################################################
    mkdir -p ./Dataset_QC-Visualization

    echo "==========================================================="
    printf "Perfoming sample QC Step 1 -- Checking sex information\n"
    echo "==========================================================="
    # require X chromosome
    cat Dups2Remove_${RawData}.CASE.bim | awk '{if ($1==23) print $0}' > X_chromosome.txt
    if [ -s X_chromosome.txt ]
    then
        plink1.9 --bfile Dups2Remove_${RawData}.CASE \
                 --check-sex \
                 --out ./Dataset_QC-Visualization/sexstat
        grep "PROBLEM" ./Dataset_QC-Visualization/sexstat.sexcheck > ./Dataset_QC-Visualization/problem_sex_check.txt || true
    else
        echo -n > ./Dataset_QC-Visualization/problem_sex_check.txt
    fi
    cat ./Dataset_QC-Visualization/problem_sex_check.txt | tail -n +2 | awk '{if ($3 != 0 && $4 !=0) print $1, $2}' > ./Dataset_QC-Visualization/fail_sex_check.txt
    cat ./Dataset_QC-Visualization/problem_sex_check.txt | tail -n +2 | awk '{if ($3 != 0 && $4 !=0) print $1, $2, $3, $4, $6}' > ./Dataset_QC-Visualization/fail_sex_check_with_value.txt

    echo "==========================================================="
    printf "Perfoming sample QC Step 2 -- High missingness\n"
    echo "==========================================================="
    plink1.9 --bfile Dups2Remove_${RawData}.CASE \
             --missing \
             --out  ./Dataset_QC-Visualization/Missing
    cat ./Dataset_QC-Visualization/Missing.imiss | tail -n +2 | awk -v MindQC=$MindQC '$6 > MindQC {print $1, $2}' > ./Dataset_QC-Visualization/fail_missingness.txt

    echo "==========================================================="
    printf "Perfoming sample QC Step 3 -- Calculating heterozygosity\n"
    echo "==========================================================="
    plink1.9 --bfile Dups2Remove_${RawData}.CASE \
             --het \
             --out ./Dataset_QC-Visualization/Het
    Rscript ../analyze_heterozygosity.R \
            ./Dataset_QC-Visualization/Het.het \
            ./Dataset_QC-Visualization/heterozygosity_analysis.txt \
            ./Dataset_QC-Visualization/fail_het.txt \
            $Heter_SD

    # plot Heterozygosity vs Individual Missingness
    Rscript ../imiss-vs-het.R \
            ./Dataset_QC-Visualization/Het.het \
            ./Dataset_QC-Visualization/Missing.imiss \
            ./Dataset_QC-Visualization/het_vs_imiss.jpg \
            $MindQC \
            $Heter_SD

    echo "==========================================================="
    printf "Perfoming sample QC Step 4 -- IBD test\n"
    echo "==========================================================="
    plink1.9 --bfile Dups2Remove_${RawData}.CASE \
             --indep-pairwise 50 5 0.2 \
             --out ./Dataset_QC-Visualization/high_LD
    # calculate pairwise IBS
    plink1.9 --bfile Dups2Remove_${RawData}.CASE \
             --extract ./Dataset_QC-Visualization/high_LD.prune.in \
             --genome \
             --out ./Dataset_QC-Visualization/RelatedTest 
    # plot pairwise IBS
    Rscript ../IBD_plot.R \
            ./Dataset_QC-Visualization/RelatedTest.genome \
            ./Dataset_QC-Visualization/failed_related.txt \
            ./Dataset_QC-Visualization/IBD_heatmap.jpg \
            $IBDQC

    cat ./Dataset_QC-Visualization/failed_related.txt | tail -n +2 | awk '{print $1, $2}' > ./Dataset_QC-Visualization/fail_IBD_check.txt

    echo "==========================================================="
    printf "Perfoming sample QC final step -- Removing failed samples\n"
    echo "==========================================================="
    # join samples (fail_sex_check.txt, fail_missingness.txt, fail_het.txt, fail_IBD_check.txt)
    cd ./Dataset_QC-Visualization
    cat fail_sex_check.txt fail_missingness.txt fail_het.txt fail_IBD_check.txt | sort -k1 | uniq > fail_example_inds.FAIL_CASE.txt
    cd ../
    # Remove samples
    plink1.9 --bfile Dups2Remove_${RawData}.CASE \
             --remove ./Dataset_QC-Visualization/fail_example_inds.FAIL_CASE.txt \
             --make-bed \
             --out ${RawData}_Pre-ImputeQC1.CASE

    ##########################################################
    ######################## SNP QC ##########################
    ##########################################################
    echo "==========================================================="
    printf "Perfoming SNP QC Step 1 -- Minor Allele Frequency\n"
    echo "==========================================================="
    plink1.9 --bfile ${RawData}_Pre-ImputeQC1.CASE \
             --freq \
             --out ./Dataset_QC-Visualization/MAF_freq
    Rscript ../MAF_plot.R \
            ./Dataset_QC-Visualization/MAF_freq.frq \
            ./Dataset_QC-Visualization/MAF_plot.jpg \
            $MafQC

    echo "==========================================================="
    printf "Perfoming SNP QC Step 2 -- High Missingness\n"
    echo "==========================================================="
    plink1.9 --bfile ${RawData}_Pre-ImputeQC1.CASE \
             --missing \
             --out ./Dataset_QC-Visualization/locus_Missing
    Rscript ../locus_Missing_plot.R \
            ./Dataset_QC-Visualization/locus_Missing.lmiss \
            ./Dataset_QC-Visualization/locus_Missing_plot.jpg \
            $GenoQC

    echo "==========================================================="
    printf "Perfoming SNP QC Step 3 -- Hardy Weinberg Equilibrium\n"
    echo "==========================================================="
    plink1.9 --bfile ${RawData}_Pre-ImputeQC1.CASE \
             --hardy \
             --out ./Dataset_QC-Visualization/hwe
    Rscript ../hwe_plot.R \
            ./Dataset_QC-Visualization/hwe.hwe \
            ./Dataset_QC-Visualization/hwe_plot.jpg \
            $HweQC


    echo "==========================================================="
    printf "Perfoming SNP QC Step 4 -- Removing Poorly Genotyped Genotypes\n"
    echo "==========================================================="
    plink1.9 --allow-no-sex \
             --bfile ${RawData}_Pre-ImputeQC1.CASE \
             --geno ${GenoQC} \
             --maf ${MafQC} \
             --hwe ${HweQC} \
             --make-bed \
             --out ${RawData}_Pre-ImputeQC2.CASE


    echo "==========================================================="
    printf "remove the X and Y chromosomes\n"
    echo "==========================================================="	
    plink1.9 --bfile ${RawData}_Pre-ImputeQC2.CASE \
             --chr 1-22 \
             --make-bed \
             --out ${RawData}_Pre-ImputeQC3.CASE


    if [ $Population_Stratification == "None" ]
    then
        echo "==========================================================="
        printf "Save processed output file\n"
        echo "==========================================================="
        mv ${RawData}_Pre-ImputeQC3.CASE.bed $OUT_CASE_BED
        mv ${RawData}_Pre-ImputeQC3.CASE.bim $OUT_CASE_BIM
        mv ${RawData}_Pre-ImputeQC3.CASE.fam $OUT_CASE_FAM
    else
        echo "==========================================================="
        printf "Performing population stratification via IBS values\n"
        echo "==========================================================="	
        echo "Target population: $Population_Stratification"
        echo "Preparing Hapmap3 reference"
        printf ".\n.\n.\n"
        cp $HM_BED ./hapmap.bed
        cp $HM_BIM ./hapmap.bim
        cp $HM_FAM ./hapmap.fam
        cp $HM_INFO ./hapmap_info.txt
        plink1.9 --bfile hapmap \
                 --list-duplicate-vars ids-only suppress-first \
                 --out ./Dups2Remove_hapmap \
                 --allow-no-sex
        if [[ $(wc -l <Dups2Remove_hapmap.dupvar) > 0 ]]
        then 
            plink1.9 --bfile hapmap \
                    --exclude ./Dups2Remove_hapmap.dupvar \
                    --make-bed \
                    --out hapmap_Dups2Remove \
                    --allow-no-sex
            mv hapmap_Dups2Remove.bed hapmap.bed
            mv hapmap_Dups2Remove.bim hapmap.bim
            mv hapmap_Dups2Remove.fam hapmap.fam
        fi
        printf "Done \n"
        echo "==========================================================="
        printf "STEP 1: Bim file modification (SNPs ID)\n"
        echo "==========================================================="
        # for CASE
        cat ${RawData}_Pre-ImputeQC3.CASE.bim | awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > CASE_PCA_Step1.bim
        cp ${RawData}_Pre-ImputeQC3.CASE.bed CASE_PCA_Step1.bed
        cp ${RawData}_Pre-ImputeQC3.CASE.fam CASE_PCA_Step1.fam
        # for hapmap 
        cat hapmap.bim | awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > hapmap_PCA_Step1.bim
        mv hapmap.bed hapmap_PCA_Step1.bed
        mv hapmap.fam hapmap_PCA_Step1.fam

        echo "==========================================================="
        printf "STEP 2: extract common snps (CASE and hapmap)\n"
        echo "==========================================================="
        grep -f CASE_PCA_Step1.bim hapmap_PCA_Step1.bim | cut -f2 > common_snps.txt
        # for CASE
        plink1.9 --bfile CASE_PCA_Step1 \
                 --extract common_snps.txt \
                 --make-bed \
                 --out CASE_PCA_Step2

        # for hapmap
        plink1.9 --bfile hapmap_PCA_Step1 \
                 --extract common_snps.txt \
                 --make-bed \
                 --out hapmap_PCA_Step2

        # merge bfile (keep executing when "exclude 3+ alleles" error occurs)
        plink1.9 --bfile CASE_PCA_Step2 --bmerge hapmap_PCA_Step2 --make-bed --out CASE_hapmap_Step3 --allow-no-sex || true
        echo "==========================================================="
        printf "STEP 3: exclude 3+ alleles (CASE and hapmap)\n"
        echo "==========================================================="
        if [ -s CASE_hapmap_Step3-merge.missnp ]
        then
            # for CASE
            plink1.9 --bfile CASE_PCA_Step2 \
                    --exclude CASE_hapmap_Step3-merge.missnp \
                    --make-bed \
                    --out CASE_PCA_Step3

            # for hapmap
            plink1.9 --bfile hapmap_PCA_Step2 \
                    --exclude CASE_hapmap_Step3-merge.missnp \
                    --make-bed \
                    --out hapmap_PCA_Step3

            # merge bfile
            plink1.9 --bfile CASE_PCA_Step3 \
                    --bmerge hapmap_PCA_Step3 \
                    --make-bed \
                    --out CASE_hapmap_Step3 \
                    --allow-no-sex
        fi
        echo "==========================================================="
        printf "STEP 4: exclude INDEL (insertion deletion)\n"
        echo "==========================================================="
        # extract SNPs with INDEL
        cat CASE_hapmap_Step3.bim | awk '{if (length($5)>1 || length($6)>1) print $2}' > INDEL.txt
        plink1.9 --bfile CASE_hapmap_Step3 \
                 --exclude INDEL.txt \
                 --make-bed \
                 --out CASE_hapmap_Step4 \
                 --allow-no-sex

        echo "==========================================================="
        printf "STEP 5: PCA (CASE)\n"
        echo "==========================================================="
        plink1.9 --bfile CASE_hapmap_Step4 \
                --indep-pairwise 50 10 0.2 \
                --out high_LD

        plink1.9 --bfile CASE_hapmap_Step4 \
                --extract high_LD.prune.in \
                --genome \
                --out RelatedTest

        plink1.9 --bfile CASE_hapmap_Step4 \
                --read-genome RelatedTest.genome \
                --cluster cc \
                --ppc 1e-3 \
                --mds-plot 2 \
                --out pop_strat \
                --allow-no-sex
        # PCA
        Rscript ../PCA.R pop_strat.mds \
                         CASE_PCA_Step3.fam \
                         $HM_INFO \
                         $Population_Stratification \
                         $OUT_DIR/WorkingDir/CASE/

        if [ -f fail_PCA_list.txt ]
        then
            echo "==========================================================="
            printf "STEP 6: exclude samples belonging other populations\n"
            echo "==========================================================="
            cat fail_PCA_list.txt | awk 'NR>1{print $1"\t"$2}' > fail_PCA_samples.txt
            plink1.9 --bfile ${RawData}_Pre-ImputeQC3.CASE \
                    --exclude fail_PCA_samples.txt \
                    --make-bed \
                    --out CASE_hapmap_Step5 \
                    --allow-no-sex
            mv CASE_hapmap_Step5.bed $OUT_CASE_BED
            mv CASE_hapmap_Step5.bim $OUT_CASE_BIM
            mv CASE_hapmap_Step5.fam $OUT_CASE_FAM
            mv population_stratification.jpg $OUT_DIR/population_stratification_CASE.jpg
            mv fail_PCA_list.txt $OUT_DIR/fail_PCA_CASE.txt
        else
            mv ${RawData}_Pre-ImputeQC3.CASE.bed $OUT_CASE_BED
            mv ${RawData}_Pre-ImputeQC3.CASE.bim $OUT_CASE_BIM
            mv ${RawData}_Pre-ImputeQC3.CASE.fam $OUT_CASE_FAM
            mv population_stratification.jpg $OUT_DIR/population_stratification_CASE.jpg
        fi
    fi

    echo 
    echo "==========================================================="
    printf "Generate QC report\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
    # python code####################
    cp ../logo.png .
    python3 ../QC_report.py --Heter $Heter_SD \
                            --Ind_Miss $MindQC \
                            --Geno_Miss $GenoQC \
                            --Maf $MafQC \
                            --HWE $HweQC \
                            --Ind_Dec 0.1875 \
                            --Work_Dir $OUT_DIR/WorkingDir/CASE/ \
                            --Basename $RawData \
                            --CaseControl CASE
    printf "Done \n"
    mv pre_QC_report.CASE.pdf $OUT_DIR
    cd ..
fi

# ===============================================================================
# =============================== CONTROL GROUP =================================
# ===============================================================================
if [ -f ${RawData}.CONTROL.bed ]
then
    echo "Processing 'CONTROL' group..."
    mkdir -p ./CONTROL
    mv ${RawData}.CONTROL.* ./CONTROL
    # Genotype Missingness 
    GenoQC=$CONTROL_GENO_QC;
    # Individual Missingness 
    MindQC=$CONTROL_MIND_QC;
    # Hardy Weinberg Equilibrium
    HweQC=$CONTROL_HWE_QC;
    # Minor Allele Frequency 
    MafQC=$CONTROL_MAF_QC;
    # IBD
    IBDQC=$CONTROL_IBD_QC;
    # standard deviation for heterozygosity
    Heter_SD=$CONTROL_HETER_SD;

    #Make CONTROL Directory in which all CONTROL files will populate
    cd ./CONTROL

    # if CONTROL only, fix allele based on reference
    if [ $Fix_Input_with_Reference_Annotation == "TRUE" ]
    then
        # --------------------------------------------------------------------------------------------------
        # STEP 1: Convert the Plink BED/BIM/FAM into a VCF into a BCF so that it may be fixed with BCFtools
        # --------------------------------------------------------------------------------------------------
        echo 
        echo "==========================================================="
        printf "Converting Plink files into VCF format \n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            plink1.9 --allow-no-sex \
                     --bfile ${RawData}.CONTROL \
                     --recode vcf \
                     --out ./DataFixStep1_${RawData}.CONTROL
        printf "Done \n"

        echo 
        echo "==========================================================="
        printf "Converting VCF into a BCF with chromosome names that match the reference .fasta annotation\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools annotate -Ob --rename-chrs ../RefAnnotationData/PlinkChrRename.txt ./DataFixStep1_${RawData}.CONTROL.vcf > ./DataFixStep1_${RawData}.CONTROL.bcf
        printf "Done \n"
        # --------------------------------------------------------------------------------------------------
        # STEP 2: Align Input File to the Reference Annotation (Fix with BCFtools)
        # --------------------------------------------------------------------------------------------------
        echo 
        echo "==========================================================="
        printf "View the number of correctly annotated/aligned variants to the Reference annotation\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep1_${RawData}.CONTROL.bcf -- -f ../RefAnnotationData/*.fasta.gz 
        printf "Done \n"

        echo 
        echo "==========================================================="
        printf "fix the strand if the VCF is TOP-compatible"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep1_${RawData}.CONTROL.bcf -Ob -o ./DataFixStep2_${RawData}.CONTROL.bcf -- -f ../RefAnnotationData/*.fasta.gz -m top
        printf "Done \n"
        echo 
        
        echo "==========================================================="
        printf "Swap the alleles based on the downloaded annotation file\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep2_${RawData}.CONTROL.bcf -Ob -o ./DataFixStep3_${RawData}.CONTROL.bcf -- -d -f ../RefAnnotationData/*.fasta.gz -i ../RefAnnotationData/*.vcf.gz
        printf "Done \n"
        
        echo 
        echo "==========================================================="
        printf "Sorting the BCF output since fixing it may have made it unsorted\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools sort ./DataFixStep3_${RawData}.CONTROL.bcf -Ob -o ./DataFixStep4_${RawData}.CONTROL.bcf
        printf "Done \n"


        echo 
        echo "==========================================================="
        printf "Check to see if the file has been fixed\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            bcftools +fixref ./DataFixStep4_${RawData}.CONTROL.bcf -Ob -o ./DataFixStep5_${RawData}.CONTROL.bcf -- -f ../RefAnnotationData/*.fasta.gz -m flip -d
        printf "Done \n"


        # Convert BCF back into Plink .bed/.bim/.fam for Shapeit2 Phasing
        echo 
        echo "==========================================================="
        printf "Converting Fixed and Sorted BCF back into Plink .bed/.bim/.fam\n"
        echo "==========================================================="
        printf ".\n.\n.\n"
            plink1.9 --bcf ./DataFixStep5_${RawData}.CONTROL.bcf \
                     --make-bed \
                     --out ./DataFixStep6_${RawData}.CONTROL
        printf "Done \n"
        cp ${RawData}.CONTROL.fam DataFixStep6_${RawData}.CONTROL.fam

        # --------------------------------------------------------------------------------------------------
        # Add back in the sex information
        # --------------------------------------------------------------------------------------------------
        # echo 
        # echo "==========================================================="
        # printf "Restoring Sample Sex Information\n"
        # echo "==========================================================="	
        # printf ".\n.\n.\n"
        #     plink1.9 --bfile ./DataFixStep7_${RawData}.CONTROL \
        #              --update-sex ${RawData}.CONTROL.fam 3 \
        #              --allow-no-sex \
        #              --make-bed \
        #              --out ./DataFixStep8_${RawData}.CONTROL > sex_restoring.log
        # printf "Done \n"
        mv ./DataFixStep6_${RawData}.CONTROL.bed ${RawData}.CONTROL.bed
        mv ./DataFixStep6_${RawData}.CONTROL.bim ${RawData}.CONTROL.bim
        mv ./DataFixStep6_${RawData}.CONTROL.fam ${RawData}.CONTROL.fam
        rm DataFixStep*
    fi 

    # Finally Remove any positional duplicates 
    # i.e. same position and alleles, but differently named variants since Shapeit will not tolerate these
    echo 
    echo "==========================================================="
    printf "Finding Positional and Allelic Duplicates\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
        plink1.9 --bfile ${RawData}.CONTROL \
                 --list-duplicate-vars ids-only suppress-first \
                 --out ./Dups2Remove \
                 --allow-no-sex
    printf "Done \n"

    # --------------------------------------------------------------------------------------------------
    # Report Number of duplicates:
    # --------------------------------------------------------------------------------------------------
    echo 
    echo "==========================================================="
    printf "Removing Positional and Allelic Duplicates if they exist\n"
    echo "==========================================================="
        DuplicateNumber="$(wc ./Dups2Remove.dupvar | awk '{print $1}')"
    printf "Found Number of Duplicate Variant(s):\n${DuplicateNumber}"
    printf ".\n.\n.\n"
        plink1.9 --bfile ${RawData}.CONTROL \
                 --exclude ./Dups2Remove.dupvar \
                 --make-bed \
                 --out Dups2Remove_${RawData}.CONTROL \
                 --allow-no-sex
    printf "Done \n"
    # --------------------------------------------------------------------------------------------------
    # STEP 6: Pre-imputation QC
    # --------------------------------------------------------------------------------------------------
    ##########################################################
    ####################### Sample QC ########################
    #########################################################
    mkdir -p ./Dataset_QC-Visualization

    echo "==========================================================="
    printf "Perfoming sample QC Step 1 -- Checking sex information\n"
    echo "==========================================================="
    # require X chromosome
    cat Dups2Remove_${RawData}.CONTROL.bim | awk '{if ($1==23) print $0}' > X_chromosome.txt
    if [ -s X_chromosome.txt ]
    then
        plink1.9 --bfile Dups2Remove_${RawData}.CONTROL \
                 --check-sex \
                 --out ./Dataset_QC-Visualization/sexstat
        grep "PROBLEM" ./Dataset_QC-Visualization/sexstat.sexcheck > ./Dataset_QC-Visualization/problem_sex_check.txt || true
    else
        echo -n > ./Dataset_QC-Visualization/problem_sex_check.txt
    fi
    cat ./Dataset_QC-Visualization/problem_sex_check.txt | tail -n +2 | awk '{if ($3 != 0 && $4 !=0) print $1, $2}' > ./Dataset_QC-Visualization/fail_sex_check.txt
    cat ./Dataset_QC-Visualization/problem_sex_check.txt | tail -n +2 | awk '{if ($3 != 0 && $4 !=0) print $1, $2, $3, $4, $6}' > ./Dataset_QC-Visualization/fail_sex_check_with_value.txt

    echo "==========================================================="
    printf "Perfoming sample QC Step 2 -- High missingness\n"
    echo "==========================================================="
    plink1.9 --bfile Dups2Remove_${RawData}.CONTROL \
             --missing \
             --out  ./Dataset_QC-Visualization/Missing
    cat ./Dataset_QC-Visualization/Missing.imiss | tail -n +2 | awk -v MindQC=$MindQC '$6 > MindQC {print $1, $2}' > ./Dataset_QC-Visualization/fail_missingness.txt

    echo "==========================================================="
    printf "Perfoming sample QC Step 3 -- Calculating heterozygosity\n"
    echo "==========================================================="
    plink1.9 --bfile Dups2Remove_${RawData}.CONTROL \
             --het \
             --out ./Dataset_QC-Visualization/Het
    Rscript ../analyze_heterozygosity.R \
            ./Dataset_QC-Visualization/Het.het \
            ./Dataset_QC-Visualization/heterozygosity_analysis.txt \
            ./Dataset_QC-Visualization/fail_het.txt \
            $Heter_SD

    # plot Heterozygosity vs Individual Missingness
    Rscript ../imiss-vs-het.R \
            ./Dataset_QC-Visualization/Het.het \
            ./Dataset_QC-Visualization/Missing.imiss \
            ./Dataset_QC-Visualization/het_vs_imiss.jpg \
            $MindQC \
            $Heter_SD

    echo "==========================================================="
    printf "Perfoming sample QC Step 4 -- IBD test\n"
    echo "==========================================================="
    plink1.9 --bfile Dups2Remove_${RawData}.CONTROL \
             --indep-pairwise 50 5 0.2 \
             --out ./Dataset_QC-Visualization/high_LD
    # calculate pairwise IBS
    plink1.9 --bfile Dups2Remove_${RawData}.CONTROL \
             --extract ./Dataset_QC-Visualization/high_LD.prune.in \
             --genome \
             --out ./Dataset_QC-Visualization/RelatedTest 
    # plot pairwise IBS
    Rscript ../IBD_plot.R \
            ./Dataset_QC-Visualization/RelatedTest.genome \
            ./Dataset_QC-Visualization/failed_related.txt \
            ./Dataset_QC-Visualization/IBD_heatmap.jpg \
            $IBDQC

    cat ./Dataset_QC-Visualization/failed_related.txt | tail -n +2 | awk '{print $1, $2}' > ./Dataset_QC-Visualization/fail_IBD_check.txt

    echo "==========================================================="
    printf "Perfoming sample QC final step -- Removing failed samples\n"
    echo "==========================================================="
    # join samples (fail_sex_check.txt, fail_missingness.txt, fail_het.txt, fail_IBD_check.txt)
    cd ./Dataset_QC-Visualization
    cat fail_sex_check.txt fail_missingness.txt fail_het.txt fail_IBD_check.txt | sort -k1 | uniq > fail_example_inds.FAIL_CONTROL.txt
    cd ../
    # Remove samples
    plink1.9 --bfile Dups2Remove_${RawData}.CONTROL \
             --remove ./Dataset_QC-Visualization/fail_example_inds.FAIL_CONTROL.txt \
             --make-bed \
             --out ${RawData}_Pre-ImputeQC1.CONTROL

    ##########################################################
    ######################## SNP QC ##########################
    ##########################################################
    echo "==========================================================="
    printf "Perfoming SNP QC Step 1 -- Minor Allele Frequency\n"
    echo "==========================================================="
    plink1.9 --bfile ${RawData}_Pre-ImputeQC1.CONTROL \
             --freq \
             --out ./Dataset_QC-Visualization/MAF_freq
    Rscript ../MAF_plot.R \
            ./Dataset_QC-Visualization/MAF_freq.frq \
            ./Dataset_QC-Visualization/MAF_plot.jpg \
            $MafQC

    echo "==========================================================="
    printf "Perfoming SNP QC Step 2 -- High Missingness\n"
    echo "==========================================================="
    plink1.9 --bfile ${RawData}_Pre-ImputeQC1.CONTROL \
             --missing \
             --out ./Dataset_QC-Visualization/locus_Missing
    Rscript ../locus_Missing_plot.R \
            ./Dataset_QC-Visualization/locus_Missing.lmiss \
            ./Dataset_QC-Visualization/locus_Missing_plot.jpg \
            $GenoQC

    echo "==========================================================="
    printf "Perfoming SNP QC Step 3 -- Hardy Weinberg Equilibrium\n"
    echo "==========================================================="
    plink1.9 --bfile ${RawData}_Pre-ImputeQC1.CONTROL \
            --hardy \
            --out ./Dataset_QC-Visualization/hwe
    Rscript ../hwe_plot.R \
            ./Dataset_QC-Visualization/hwe.hwe \
            ./Dataset_QC-Visualization/hwe_plot.jpg \
            $HweQC


    echo "==========================================================="
    printf "Perfoming SNP QC Step 4 -- Removing Poorly Genotyped Genotypes\n"
    echo "==========================================================="
    plink1.9 --allow-no-sex \
             --bfile ${RawData}_Pre-ImputeQC1.CONTROL \
             --geno ${GenoQC} \
             --maf ${MafQC} \
             --hwe ${HweQC} \
             --make-bed \
             --out ${RawData}_Pre-ImputeQC2.CONTROL


    echo "==========================================================="
    printf "remove the X and Y chromosomes\n"
    echo "==========================================================="	
    plink1.9 --bfile ${RawData}_Pre-ImputeQC2.CONTROL \
             --chr 1-22 \
             --make-bed \
             --out ${RawData}_Pre-ImputeQC3.CONTROL

    if [ $Population_Stratification == "None" ]
    then
        echo "==========================================================="
        printf "Save processed output file\n"
        echo "==========================================================="
        mv ${RawData}_Pre-ImputeQC3.CONTROL.bed $OUT_CONTROL_BED
        mv ${RawData}_Pre-ImputeQC3.CONTROL.bim $OUT_CONTROL_BIM
        mv ${RawData}_Pre-ImputeQC3.CONTROL.fam $OUT_CONTROL_FAM
    else
        echo "==========================================================="
        printf "Performing population stratification via IBS values\n"
        echo "==========================================================="	
        echo "Target population: $Population_Stratification"
        echo "Preparing Hapmap3 reference"
        printf ".\n.\n.\n"
        cp $HM_BED ./hapmap.bed
        cp $HM_BIM ./hapmap.bim
        cp $HM_FAM ./hapmap.fam
        cp $HM_INFO ./hapmap_info.txt
        plink1.9 --bfile hapmap \
                 --list-duplicate-vars ids-only suppress-first \
                 --out ./Dups2Remove_hapmap \
                 --allow-no-sex
        if [[ $(wc -l <Dups2Remove_hapmap.dupvar) > 0 ]]
        then 
            plink1.9 --bfile hapmap \
                     --exclude ./Dups2Remove_hapmap.dupvar \
                     --make-bed \
                     --out hapmap_Dups2Remove \
                     --allow-no-sex
            mv hapmap_Dups2Remove.bed hapmap.bed
            mv hapmap_Dups2Remove.bim hapmap.bim
            mv hapmap_Dups2Remove.fam hapmap.fam
        fi
        printf "Done \n"
        echo "==========================================================="
        printf "STEP 1: Bim file modification (SNPs ID)\n"
        echo "==========================================================="
        # for CONTROL
        cat ${RawData}_Pre-ImputeQC3.CONTROL.bim | awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > CONTROL_PCA_Step1.bim
        cp ${RawData}_Pre-ImputeQC3.CONTROL.bed CONTROL_PCA_Step1.bed
        cp ${RawData}_Pre-ImputeQC3.CONTROL.fam CONTROL_PCA_Step1.fam
        # for hapmap 
        cat hapmap.bim | awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > hapmap_PCA_Step1.bim
        mv hapmap.bed hapmap_PCA_Step1.bed
        mv hapmap.fam hapmap_PCA_Step1.fam

        echo "==========================================================="
        printf "STEP 2: extract common snps (CONTROL and hapmap)\n"
        echo "==========================================================="
        grep -f CONTROL_PCA_Step1.bim hapmap_PCA_Step1.bim | cut -f2 > common_snps.txt
        # for CONTROL
        plink1.9 --bfile CONTROL_PCA_Step1 \
                 --extract common_snps.txt \
                 --make-bed \
                 --out CONTROL_PCA_Step2

        # for hapmap
        plink1.9 --bfile hapmap_PCA_Step1 \
                 --extract common_snps.txt \
                 --make-bed \
                 --out hapmap_PCA_Step2

        # merge bfile (keep executing when "exclude 3+ alleles" error occurs)
        plink1.9 --bfile CONTROL_PCA_Step2 --bmerge hapmap_PCA_Step2 --make-bed --out CONTROL_hapmap_Step3 --allow-no-sex || true
        echo "==========================================================="
        printf "STEP 3: exclude 3+ alleles (CONTROL and hapmap)\n"
        echo "==========================================================="
        if [ -s CONTROL_hapmap_Step3-merge.missnp ]
        then
            # for CONTROL
            plink1.9 --bfile CONTROL_PCA_Step2 \
                     --exclude CONTROL_hapmap_Step3-merge.missnp \
                     --make-bed \
                     --out CONTROL_PCA_Step3

            # for hapmap
            plink1.9 --bfile hapmap_PCA_Step2 \
                     --exclude CONTROL_hapmap_Step3-merge.missnp \
                     --make-bed \
                     --out hapmap_PCA_Step3

            # merge bfile
            plink1.9 --bfile CONTROL_PCA_Step3 \
                     --bmerge hapmap_PCA_Step3 \
                     --make-bed \
                     --out CONTROL_hapmap_Step3 \
                     --allow-no-sex
        fi
        echo "==========================================================="
        printf "STEP 4: exclude INDEL (insertion deletion)\n"
        echo "==========================================================="
        # extract SNPs with INDEL
        cat CONTROL_hapmap_Step3.bim | awk '{if (length($5)>1 || length($6)>1) print $2}' > INDEL.txt
        plink1.9 --bfile CONTROL_hapmap_Step3 \
                 --exclude INDEL.txt \
                 --make-bed \
                 --out CONTROL_hapmap_Step4 \
                 --allow-no-sex

        echo "==========================================================="
        printf "STEP 5: PCA (CONTROL)\n"
        echo "==========================================================="
        plink1.9 --bfile CONTROL_hapmap_Step4 \
                 --indep-pairwise 50 10 0.2 \
                 --out high_LD

        plink1.9 --bfile CONTROL_hapmap_Step4 \
                 --extract high_LD.prune.in \
                 --genome \
                 --out RelatedTest

        plink1.9 --bfile CONTROL_hapmap_Step4 \
                 --read-genome RelatedTest.genome \
                 --cluster cc \
                 --ppc 1e-3 \
                 --mds-plot 2 \
                 --out pop_strat \
                 --allow-no-sex
        # PCA
        Rscript ../PCA.R pop_strat.mds \
                         CONTROL_PCA_Step3.fam \
                         $HM_INFO \
                         $Population_Stratification \
                         $OUT_DIR/WorkingDir/CONTROL/
        if [ -f fail_PCA_list.txt ]
        then
            echo "==========================================================="
            printf "STEP 6: exclude samples belonging other populations\n"
            echo "==========================================================="
            cat fail_PCA_list.txt | awk '{print $1"\t"$2}' > fail_PCA_samples.txt
            plink1.9 --bfile ${RawData}_Pre-ImputeQC3.CONTROL \
                     --exclude fail_PCA_samples.txt \
                     --make-bed \
                     --out CONTROL_hapmap_Step5 \
                     --allow-no-sex
            mv CONTROL_hapmap_Step5.bed $OUT_CONTROL_BED
            mv CONTROL_hapmap_Step5.bim $OUT_CONTROL_BIM
            mv CONTROL_hapmap_Step5.fam $OUT_CONTROL_FAM
            mv population_stratification.jpg $OUT_DIR/population_stratification_CONTROL.jpg
            mv fail_PCA_list.txt $OUT_DIR/fail_PCA_CONTROL.txt
        else
            mv ${RawData}_Pre-ImputeQC3.CONTROL.bed $OUT_CONTROL_BED
            mv ${RawData}_Pre-ImputeQC3.CONTROL.bim $OUT_CONTROL_BIM
            mv ${RawData}_Pre-ImputeQC3.CONTROL.fam $OUT_CONTROL_FAM
            mv population_stratification.jpg $OUT_DIR/population_stratification_CONTROL.jpg
        fi
    fi

    # --------------------------------------------------------------------------------------------------
    # Visualize genomic data for missingness, heterozygosity, and relatedness
    # --------------------------------------------------------------------------------------------------
    echo 
    echo "==========================================================="
    printf "Generate QC report\n"
    echo "==========================================================="
    printf ".\n.\n.\n"
    # python code####################
    cp ../logo.png .
    python3 ../QC_report.py --Heter $Heter_SD \
                            --Ind_Miss $MindQC \
                            --Geno_Miss $GenoQC \
                            --Maf $MafQC \
                            --HWE $HweQC \
                            --Ind_Dec 0.1875 \
                            --Work_Dir $OUT_DIR/WorkingDir/CONTROL/ \
                            --Basename $RawData \
                            --CaseControl CONTROL
    printf "Done \n"
    mv pre_QC_report.CONTROL.pdf $OUT_DIR

fi

cd $OUT_DIR
rm -r WorkingDir 



