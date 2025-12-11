## June 2025 - edited to rarefy data
## 12 July 2025 - need to remove outlier growth samples 106, 95 (if present - wasn't in fungi), 142


### updated CH2 bacterial analysis
### recreating NMDS (updated from RCV's PCOAs) and PERMANOVA in adonis2
## possibly also a taxonomic barplot or two

library(ggplot2)
library(vegan)
library(car)

#libraries related to plotting (and data manipulation to prep for plotting)
library(dplyr)
library(tidyverse)
library(grid)
library(ggpubr)

#setwd('C:\\Files\\Box Sync\\PhD\\1 Research projects\\Chapter2 - PSF\\Code\\bacteria') 
setwd('C:\\Files\\PhD\\1 Research projects\\Chapter2 - PSF\\Code\\bacteria')
psfb<-read.table(file = 'ASV_name_psf_all_samples.txt', sep = '\t', header = TRUE)

#meta<-read.table(file = 'CH2_metadata_newIDs3.txt', sep = '\t', header = TRUE)
meta<-read.table(file = 'CH2_metadata_newIDs4.txt', sep = '\t', header = TRUE) #adds species_soil_source concat

## make column to match my sample names with RCV's sample names

psfb$intermediate<-rownames(psfb)
psfb$num <- gsub(x = psfb$intermediate, pattern = "PSF_soil_", replacement = "")
psfb$sample<-paste("Sample",psfb$num) ## there are SPACES in here
psfb$sample <- gsub(x = psfb$sample, pattern = " ", replacement = "")

df<-merge(psfb,meta,by="sample")

#make sum column per sample so I know range of seq counts

write.csv(df,"df.csv") #make sum column in excel
readcounts<-read.csv("df_readcounts.csv")

#########################################################
# rarefaction

rarefy_rows<-function(count_data, n_reads) {
  count_matrix<-as.matrix(count_data)
  
  rarefied<-t(apply(count_matrix,1,function(row){
    total<-sum(row)
    if(total<n_reads){
      return(rep(NA,length(row)))
    }
    
    probs<-row/total
    rmultinom(1,size=n_reads,prob=probs)[,1]
  }))
  colnames(rarefied)<-colnames(count_matrix)
  rownames(rarefied)<-rownames(count_matrix)
  
  return(rarefied)
}

### prep the data for rarefaction
df2<-df
rownames(df2)<-df2[[1]]
df2<-df2[,-1]
df3<-df2[,1:795]
df3_matrix<-as.matrix(df3)
storage.mode(df3_matrix)<-"numeric"

## rarefy

rr<-rarefy_rows(df3_matrix,8000) #8k rarefaction eliminates three samples (those with NAs)

###############################################
## NMDS/PERMANOVA
## I think to do the NMDS/PERMANOVA etc we have to 1) remove the NA rows, and 2) merge back
## with the metadata so they have the same rows

rrdf<-as.data.frame(rr)
rrdf2<-na.omit(rrdf) 
rrdf2$sample<-rownames(rrdf2) ## rrdf2 has the sample ID column


## remove outlier samples
#106, 95 (if present - wasn't in fungi), 142

test<-rrdf2
test2<-test[!(test$sample %in% c("Sample106", "Sample142")), ] #success

rrdf2meta<-merge(test2,meta,by="sample") 
rrdf3<-test2[,1:795]


## calc distance matrices + PERMANOVA

dm <- vegdist(rrdf3, method ="bray",na.rm=TRUE)

## PERMANOVA models
m0<-adonis2(dm~soil*species*soil_source,rrdf2meta,by="margin")
m1<-adonis2(dm~soil*soil_source+species,rrdf2meta,by="margin")


##################################################### only live samples

rrdf2metaL<-rrdf2meta[rrdf2meta$soil=="Live",] # live is CAPITALIZED

dfL<-rrdf2metaL[,1:795]
rownames(dfL)<-dfL[[1]]
dfL2<-dfL[,-1]



dmL <- vegdist(dfL2, method ="bray",na.rm=TRUE)

m5<-adonis2(dmL~species+soil_source+species:soil_source,rrdf2metaL,by="margin")
m6<-adonis2(dmL~species+soil_source,rrdf2metaL,by="margin")



########################################## NMDS


## all samples
test<-metaMDS(rrdf3,distance="bray",k=3) 
scoresk3 <- as.data.frame(scores(test,display="sites"))
scoresk3$sample<-rownames(scoresk3)
NMDS_k3<-left_join(meta, scoresk3,by="sample")
NMDS_k3<-na.omit(NMDS_k3)


NMDS_k3$species_soil<-paste(NMDS_k3$species,NMDS_k3$soil)

### live samples

livek3<-metaMDS(dfL2,distance="bray",k=3) #stress ~0.17

scoresL3 <- as.data.frame(scores(livek3,display="sites"))
scoresL3$sample<-rownames(scoresL3)
liveNMDS<-left_join(meta, scoresL3,by="sample")
liveNMDS<-na.omit(liveNMDS)



fulllabs=c("Monarda (Monarda soil)","Plantago (Monarda soil)","Monarda (Plantago soil)","Plantago (Plantago soil)","Monarda (Monarda soil)","Plantago (Monarda soil)","Monarda (Plantago soil)","Plantago (Plantago soil)")

NMDS_k3$soil_species_soil_source <- paste(NMDS_k3$soil,NMDS_k3$species_soil_source)
level_order_full<-c("Live Monarda_Monarda","Live Plantago_Monarda","Live Monarda_Plantago","Live Plantago_Plantago","Sterile Monarda_Monarda","Sterile Plantago_Monarda","Sterile Monarda_Plantago","Sterile Plantago_Plantago")
NMDS_k3$soil_species_soil_source <- factor(NMDS_k3$soil_species_soil_source, levels = level_order_full) # to change levels


soillabs=c("Monarda (Monarda soil)","Plantago (Monarda soil)","Monarda (Plantago soil)","Plantago (Plantago soil)")
level_order<-c("Monarda_Monarda","Plantago_Monarda","Monarda_Plantago","Plantago_Plantago")
liveNMDS$species_soil_source <- factor(liveNMDS$species_soil_source, levels = level_order) # to change levels



pFull = ggplot(NMDS_k3, aes(x = NMDS1, y = NMDS2,colour=soil_species_soil_source)) + 
  geom_point(size = 5, stroke=1, aes( shape = soil_species_soil_source, colour = soil_species_soil_source))+
  theme_bw()  +scale_shape_manual(labels=fulllabs,values=c(16,17,16,17,1,2,1,2),name="Plant (soil origin)")+ggtitle("All samples (Bacteria)")+
  scale_colour_manual(name="Plant (soil origin)",labels=fulllabs,values=c("#002594","#E0B2CD","#00A86B","#D2C500","#002594","#E0B2CD","#00A86B","#D2C500"))


pLive = ggplot(liveNMDS, aes(x = NMDS1, y = NMDS2,colour=species_soil_source,shape=species_soil_source)) + 
  geom_point(size = 5, stroke=1, aes(shape = species_soil_source, colour = species_soil_source))+
  theme_bw()  + stat_ellipse()+scale_shape_manual(values=c(16,17,16,17),name="Plant (soil origin)",labels=soillabs)+ggtitle("Live samples (Bacteria)")+
  scale_colour_manual(name="Plant (soil origin)", labels=soillabs,values=c("#002594","#E0B2CD","#00A86B","#D2C500"))



statL <- grobTree(textGrob("stress < 0.17                   k=3", x=0.02,  y=0.05, hjust=0,
                           gp=gpar(col="black", fontsize=12)))

statF <- grobTree(textGrob("stress < 0.11   k=3", x=0.4,  y=0.05, hjust=0,
                           gp=gpar(col="black", fontsize=12)))


pFull +annotation_custom(statF)+
  theme(plot.title = element_text(size = 14)) ## default font size in dissertation

pLive +annotation_custom(statL)+
  theme(plot.title = element_text(size = 14))


