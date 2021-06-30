args <- commandArgs(trailingOnly = TRUE)
INPUT <- args[1]
OUTPUT <- args[2]
cutoff <- as.numeric(args[3])


library("ggplot2")
library("scales")
MAF_dataset <- read.table(INPUT, header=T)
MAF_dataset <- MAF_dataset[complete.cases(MAF_dataset),]
MAF_dataset <- MAF_dataset[MAF_dataset$A1 != 0, ]
MAF_dataset <- MAF_dataset[MAF_dataset$A2 != 0, ]
# SNP and MAF only
SNP_MAF <- data.frame(MAF_dataset$SNP, MAF_dataset$MAF)
names(SNP_MAF) <- c("SNP", "MAF")
SNP_MAF$type <- sapply(SNP_MAF$MAF, function(x) if(x <= cutoff) "exclude" else "include")
bar_count <- 20 #ceiling(max(SNP_MAF$MAF)/cutoff)
# plot
MAF_plot <- ggplot(SNP_MAF, aes(x=MAF, fill = type)) +
                   geom_histogram(color = "#e9ecef", breaks = cutoff*(0:bar_count))
jpeg(OUTPUT, units="in", width=10, height=10, res=600)
MAF_plot + 
  scale_x_continuous(breaks = cutoff*(0:bar_count),
                     limits=c(0, cutoff*bar_count)) + 
  scale_y_continuous(expand = c(0,0),
                     label = unit_format(unit="K", scale=1e-3),
                     limits=c(0, max(ggplot_build(MAF_plot)$data[[1]]$count)*1.3)) +
  scale_color_manual(values = c("red","blue")) +
  labs(fill = "", x="Minor Allele Frequency", y="Number of SNPs")
dev.off()

