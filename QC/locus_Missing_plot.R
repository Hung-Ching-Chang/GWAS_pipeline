args <- commandArgs(trailingOnly = TRUE)
INPUT <- args[1]
OUTPUT <- args[2]
cutoff <- as.numeric(args[3])
library("data.table")
library("ggplot2")
library("scales")
lmiss <- fread(INPUT, header = T)
lmiss_dataset <- data.frame(lmiss$SNP, lmiss$F_MISS)
names(lmiss_dataset) <- c("SNP", "F_MISS")
lmiss_dataset <- lmiss_dataset[complete.cases(lmiss_dataset),]
lmiss_dataset$type <- sapply(lmiss_dataset$F_MISS, function(x) if(x >= cutoff) "exclude" else "include")
bar_breaks <- c(0.001, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1)

missing_plot <- ggplot(lmiss_dataset, aes(x=F_MISS,  fill = type)) +
  geom_histogram(color = "#e9ecef")

jpeg(OUTPUT, units="in", width=10, height=10, res=600)
missing_plot + 
  scale_y_continuous(expand = c(0,0), label = unit_format(unit="K", scale=1e-3),
                     limits=c(0,max(ggplot_build(missing_plot)$data[[1]]$count)/2.5)) +
  scale_x_log10(breaks=bar_breaks, label=bar_breaks) + 
  geom_vline(xintercept=cutoff, colour="#BB0000", linetype="dashed") +
  labs(fill = "", x="Fraction of missing data", y="Number of SNPs") 
dev.off()