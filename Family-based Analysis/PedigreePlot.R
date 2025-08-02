library(kinship2)
library(tidyr)
library(dplyr)

variant = df2$variant_info[1]
print("Now ploting pedigree plot for ")
print(variant)
row <- result %>% filter(variant_info == variant)
ped_data <- metadata
transcript <- unlist(strsplit(row$HGVSc,split = ":"))[2]
protein <- unlist(strsplit(row$HGVSp,split = ":"))[2]
gene <- row$SYMBOL

kid_gender_col <- colnames(result)[grep("Kid.*Gender",colnames(result))]
for (k in 1:length(kid_genotype_col)) {
  kidid <- unlist(strsplit(kid_gender_col[k], split = "Gender_"))[2]
  genotype = row[1,kid_genotype_col[k]]
  ped_data$Genotype[which(ped_data$Individual == kidid)] = genotype
}

# for father and mother genotype 
f_genotype = row[1,"FatherGenotype"]
ped_data$Genotype[which(ped_data$Father ==0 & ped_data$Mother == 0 & ped_data$Sex == 1)] = f_genotype
m_genotype = row[1,"MotherGenotype"]
ped_data$Genotype[which(ped_data$Father ==0 & ped_data$Mother == 0 & ped_data$Sex == 2)] = m_genotype
ped_data$Father[ped_data$Father == "0"] <- NA
ped_data$Mother[ped_data$Mother == "0"] <- NA
ped_data$Sex <- as.integer(ped_data$Sex)
ped_data <- ped_data %>% separate(Genotype, into = c("affected","avail"),sep = "[/|]")
ped_data$Mother[ped_data$Mother == "0"] <- NA
ped_data$affected <- as.numeric(ped_data$affected)
ped_data$avail <- as.numeric(ped_data$avail)

ped <- pedigree(id = ped_data$Individual, dadid = ped_data$Father, momid = ped_data$Mother, sex = ped_data$Sex,
                affected = cbind(ped_data$affected, ped_data$avail), famid = ped_data$FID)
ped1 <- ped[familyid]

# save figure
pdf("Desktop/pedigree_plot.pdf", width=3, height=3)
plot(ped1, col=ifelse(ped_data$avail, 2, 1), cex=0.6)
title(main="Pedigree analysis",cex.main = 0.7)
mtext(paste(gene, transcript, paste0("(", protein, ")")), 
      side=3, line=0.5, cex=0.5)
dev.off()
