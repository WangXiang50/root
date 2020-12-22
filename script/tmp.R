library(tidyverse)
library(ConnectivityMap) 
library(FactoMineR)
library(grid)
library(gridExtra)
library(Matrix)

#import instance metadata
data("instances")
instances <- instances %>% 
  mutate(instanceID = rownames(instances))

#import amplitude matrix

amp <- cbind(readRDS("./data/amplitudeMatrix_1.rds"),readRDS("./data/amplitudeMatrix_2.rds"))
colnames(amp) <- gsub("X","inst_",colnames(amp))


#import effects table
effects <- read_csv("./data/drug_effects.csv")

#instance_subset
#remove "estradiol"
examine_these <- effects$cmap_name[!effects$cmap_name %in% "estradiol"]
tested_instances <- instances$instanceID[instances$cmap_name %in% examine_these]

#filter for instances (samples) from drugs that were tested
amp_inst <- amp[,tested_instances]

#run PCA
pca <- PCA(t(amp_inst),graph = FALSE,scale.unit = TRUE)

#construct a PCA data frame with metadata 
pca.plot.data <- pca$ind$coord %>% 
  as.data.frame() %>% 
  mutate(instanceID = rownames(pca$ind$coord)) %>% 
  mutate(cmap_name = instances$cmap_name[match(instanceID,instances$instanceID)])


#plot
gg_pca <- ggplot(pca.plot.data,aes(Dim.1,Dim.2,color = cmap_name))+
  geom_point(size = 4)+
  scale_color_brewer(palette = "Set1")+
  labs(title = "PCA",
       x = paste0("PC1 (",round(pca$eig[1,2],2),"%)" ),
       y = paste0("PC2 (",round(pca$eig[2,2],2),"%)" ))


# try the same plot but filter for only transition signature (TS) genes

#get TS genes
up_gns <- as.vector(read.table("data/signature_UP.grp")[,1])
dn_gns <- as.vector(read.table("data/signature_DOWN.grp")[,1])

#filter for instances (samples) from drugs that were tested AND features in TS
amp_inst_ts <- amp[c(up_gns,dn_gns),tested_instances]

#run PCA
pca_ts <- PCA(t(amp_inst_ts),graph = FALSE,scale.unit = TRUE)

#construct a PCA data frame with metadata 
pca_ts.plot.data <- pca_ts$ind$coord %>% 
  as.data.frame() %>% 
  mutate(instanceID = rownames(pca_ts$ind$coord)) %>% 
  mutate(cmap_name = instances$cmap_name[match(instanceID,instances$instanceID)])

#plot
gg_pca_ts <- ggplot(pca_ts.plot.data,aes(Dim.1,Dim.2,color = cmap_name))+
  geom_point(size = 4)+
  scale_color_brewer(palette = "Set1")+
  labs(title = "PCA (only TS genes)",
       x = paste0("PC1 (",round(pca_ts$eig[1,2],2),"%)" ),
       y = paste0("PC2 (",round(pca_ts$eig[2,2],2),"%)" ))

#plot them together 
grid.arrange(gg_pca,gg_pca_ts)

