args <- commandArgs(trailingOnly = TRUE)
INPUT <- args[1]
OUTPUT <- args[2]
HweQC <- as.numeric(args[3])
library("data.table")
library("ggplot2")
library("scales")
hwe <- fread(INPUT, header = T)
df_hwe <- data.frame(Pvalue = unique(hwe$P), ecdf = ecdf(hwe$P)(unique(hwe$P)))
df_hwe$Pvalue[df_hwe$Pvalue == 0] <- 10^(-99)
jpeg(OUTPUT, units="in", width=10, height=10, res=600)
  ggplot(data=df_hwe, mapping = aes(x = Pvalue)) + 
    geom_area(mapping = aes(y = ifelse(Pvalue >= HweQC & Pvalue <=  max(Pvalue) , ecdf, 0), fill="include")) +
    scale_x_log10(limits = c(min(df_hwe$Pvalue), max(df_hwe$Pvalue)), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_area(mapping = aes(y = ifelse(Pvalue > min(Pvalue) & Pvalue <  HweQC , ecdf, 0), fill="exclude")) +
    geom_vline(xintercept = HweQC, colour="#BB0000", linetype="dashed") +
    scale_colour_manual("", values = c("exclude" = "red", "include" = "lightblue")) +
    labs(fill = "", x = "log P(HWE)", y = "CDF") 
dev.off()