args <- commandArgs(trailingOnly = TRUE)
INPUT <- args[1]
OUTPUT_SUMMARY <- args[2]
FAILED_SAMPLE <- args[3]
Heter_SD <- as.numeric(args[4])
library("data.table")
het = fread(INPUT, header = T)
## Add a column onto the imported het variable The new column consists of the calculation of the mean Het
het$meanHet = (het$`N(NM)` - het$`O(HOM)`)/het$`N(NM)`



summstr = paste(summary(het$meanHet), collapse=", ")
cat(paste("Min, 25%, 50%, Mean, 75%, Max = ", summstr, "\n", sep=""), file = OUTPUT_SUMMARY)
cat(paste("SD = ", sd(het$meanHet), "\n", sep=""), file = OUTPUT_SUMMARY, append = T)

sdMeanHet = sd(het$meanHet)
meanMeanHet = mean(het$meanHet)
for (i in 2:6) {
  b1 = meanMeanHet+i*sdMeanHet
  b2 = meanMeanHet-i*sdMeanHet
  outliers = length(which(het$meanHet > b1 | het$meanHet < b2))
  cat(paste("heterozygosity outliers for ", i, " SD: ", outliers, "\n", sep=""), file = OUTPUT_SUMMARY, append = T)
}

b1 = meanMeanHet+Heter_SD*sdMeanHet
b2 = meanMeanHet-Heter_SD*sdMeanHet
outlier_sample <- het[which(het$meanHet > b1 | het$meanHet < b2), c("FID", "IID")]
write.table(outlier_sample, file = FAILED_SAMPLE, col.names = F, row.names = F, quote = F)

