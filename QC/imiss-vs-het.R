args <- commandArgs(trailingOnly = TRUE)
INPUT_HET <- args[1]
INPUT_imiss <- args[2]
OUTPUT <- args[3]
IndivQC <- as.numeric(args[4])
Heter_SD <- as.numeric(args[5])

library("data.table")
library("ggplot2")
# het
het = fread(INPUT_HET, header = T)
## Add a column onto the imported het variable The new column consists of the calculation of the mean Het
het$meanHet = (het$`N(NM)` - het$`O(HOM)`)/het$`N(NM)`

# individual 
imiss <- fread(INPUT_imiss, header = T)

jpeg(OUTPUT, units="in", width=10, height=10, res=600)
# Set Standard Deviation Cutoff Criteria
    LowerSD<-mean(het$meanHet) - (Heter_SD * sd(het$meanHet))
    UpperSD<-mean(het$meanHet) + (Heter_SD * sd(het$meanHet))
# Setup the Heterozygosity Plot
    hetplot <- ggplot(imiss, aes(x=imiss$F_MISS, y=het$meanHet)) + geom_point()
# Format the Heterozygosity v IndivMissingness Plot	
    hetplot +
    scale_x_log10(breaks=c(0.001,0.01,0.1,1), 
                    label=c(0.001,0.01,0.1,1)) +
    geom_bin2d(bins = 250)+
# Visualize the Cutoffs
    geom_hline(yintercept=LowerSD, colour="#BB0000", linetype="dashed")+
    geom_hline(yintercept=UpperSD, colour="#BB0000", linetype="dashed")+
    geom_vline(xintercept=IndivQC, colour="#BB0000", linetype="dashed")+
# Setup Titles	
    labs(title="", x="Proportion of Missing per Individuals", y="Heterozygosity")
dev.off()