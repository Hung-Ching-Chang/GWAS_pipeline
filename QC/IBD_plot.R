args <- commandArgs(trailingOnly = TRUE)
INPUT_GENOME <- args[1]
OUTPUT_list <- args[2]
OUTPUT_plot <- args[3]
threshold <- as.numeric(args[4])
library("data.table")
library("ggplot2")
Related <- fread(INPUT_GENOME, h=T)
# Sample of failed IBD test
failed.related = Related[which(Related$PI_HAT >= threshold), c("FID1", "IID1", "FID2", "IID2", "PI_HAT")]
write.table(failed.related, file = OUTPUT_list, col.names = T, row.names = F, quote = F, sep = "\t")

Related <- Related[, c(2,4,10)]
Related <- Related[!is.na(Related$PI_HAT),]
Related_warning <- Related[Related$PI_HAT > threshold,]
warning_sample <- unique(c(Related_warning$IID1, Related_warning$IID2))
warning_sample_index1 <- Related$IID1 %in% warning_sample
warning_sample_index2 <- Related$IID2 %in% warning_sample
warning_sample_index_all <- which(warning_sample_index1 + warning_sample_index2 == 2)
Related_report <- Related[warning_sample_index_all,]
Related_report_symmetric <- Related_report[, c(2,1,3)]
colnames(Related_report_symmetric) <- c("IID1", "IID2", "PI_HAT")
Related_report <- rbind(Related_report, Related_report_symmetric)
Related_report$PI_HAT_disc <- cut(Related_report$PI_HAT, breaks = c(0, threshold, 0.5, 0.9, 1.1), right = FALSE)
levels(Related_report$PI_HAT_disc)[levels(Related_report$PI_HAT_disc)=="[0.9,1.1)"] <- "[0.9, 1)"

# heatmap
jpeg(OUTPUT_plot, units="in", width=10, height=10, res=600)
    colfunc<-colorRampPalette(c("white","red"))
    ggplot(Related_report, aes(x = IID1, y = IID2, fill=PI_HAT_disc)) + 
        geom_tile(aes(fill = PI_HAT_disc), colour = "grey") + 
        scale_fill_manual(values = colfunc(4))+
        theme(axis.ticks = element_blank(),
                axis.text.x  = element_text(angle=90),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())  +
        labs(fill="")
dev.off()
