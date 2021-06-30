args <- commandArgs(trailingOnly = TRUE)
In_mds <- args[1]
Input_fam <- args[2]
pop_info <- args[3]
selected_pop <- args[4]
OUT_PATH <- args[5]

library(data.table)
library(plyr)
library(ggplot2)
##############################################
# Step 2: Read example datasets
##############################################
mds <- read.table(In_mds, header=TRUE)
target <- read.table(Input_fam, header = FALSE)
hapmap <- read.table(pop_info, header = TRUE)
target_pop <- data.frame(IID = as.character(target[,"V2"]), population = "target")
hapmap_pop <- hapmap[,c("IID", "population")]
pop <- rbind(target_pop, hapmap_pop)
all_population <- rep(NA, dim(mds)[1])
for(i in 1:dim(mds)[1]){
  pos <- which(as.character(pop[,1]) == as.character(mds[i,"IID"]))
  all_population[i] <- as.character(pop[pos,2])
}
all_population = factor(all_population, levels = c("target", "CEU", "ASW", "MKK", "MEX", "CHD", 
                                                   "CHB", "JPT", "LWK", "TSI", "GIH", "YRI"))
list_colormap_race <- c("red", "#0071BC", "#A7A9AC", "#58595B", "#F26649", "#FFF460", 
                        "#FFF200", "#FFF799", "#808284", "#7DA7D8", "#0DB14B", "#171C24")
shape_list <- rep(1, dim(mds)[1])
shape_list[which(all_population == "target")] <- 4
size_list <- rep(2, dim(mds)[1])
##size_list[which(all_population == "target")] <- 4

jpeg(paste0(OUT_PATH, "population_stratification_C1_C2.jpg"), units="in", width=10, height=10, res=600)
ggplot(mds, aes(x = C1, y = C2)) +
  geom_point(shape = shape_list, size = size_list, aes(color = all_population)) + 
  scale_color_manual(breaks = levels(all_population), values=list_colormap_race)
dev.off()

jpeg(paste0(OUT_PATH, "population_stratification_C1_C3.jpg"), units="in", width=10, height=10, res=600)
ggplot(mds, aes(x = C1, y = C3)) +
  geom_point(shape = shape_list, size = size_list, aes(color = all_population)) + 
  scale_color_manual(breaks = levels(all_population), values=list_colormap_race)
dev.off()

jpeg(paste0(OUT_PATH, "population_stratification_C2_C3.jpg"), units="in", width=10, height=10, res=600)
ggplot(mds, aes(x = C2, y = C3)) +
  geom_point(shape = shape_list, size = size_list, aes(color = all_population)) + 
  scale_color_manual(breaks = levels(all_population), values=list_colormap_race)
dev.off()
#############################################
# Calculate the distance between samples and races
#############################################
if (selected_pop == 'N'){
  # population center
  pop_center <- c()
  for(race in levels(all_population)){
    if(race != "target"){
      pos <- which(all_population == race)
      pop_center <- rbind(pop_center, c(sum(mds[pos, "C1"])/length(pos), sum(mds[pos, "C2"])/length(pos), sum(mds[pos, "C3"])/length(pos)))
    }
  }

  # closest population from samples
  selected_pop_list <- unlist(strsplit(selected_pop, "-"))
  fail_PCA_list <- c()
  for(target_row in which(all_population == "target")){
    dist <- pop_center
    dist[,1] <- dist[,1] - mds[target_row, "C1"]
    dist[,2] <- dist[,2] - mds[target_row, "C2"]
    dist[,3] <- dist[,3] - mds[target_row, "C3"]
    euclidean_dist <- rowSums(dist^2)
    closest_pop <- levels(all_population)[which(euclidean_dist == min(euclidean_dist)) + 1]
    if(!closest_pop %in% selected_pop_list){
      fail_PCA_list <- rbind(fail_PCA_list, c(as.character(mds[target_row, "FID"]),
                                              as.character(mds[target_row, "IID"]),
                                              closest_pop))
    }
  }

  if(length(fail_PCA_list) != 0 ){
    colnames(fail_PCA_list) <- c("FID", "IID", "POP")
    write.table(fail_PCA_list, paste0(OUT_PATH, "fail_PCA_list.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  }
}



