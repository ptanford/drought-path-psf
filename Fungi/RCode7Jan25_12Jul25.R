## Updated PSF fungal analysis (removing outlier growth rate samples + clarifying notes for
## PSF manuscript submission)
## based on file "RCode7Jan25" (Chapter 2 sequence analysis 7 January 2025)

library(qiime2R)
library(phyloseq)
library(phyloseq.extended)
library(ggplot2)
library(vegan)
library(car)

#libraries related to plotting (and data manipulation to prep for plotting)
library(dplyr)
library(tidyverse)
library(grid)
library(ggpubr)
library(microViz)

# 1. Import qiime2 files as phyloseq object and prep ps's for analyses

setwd('C:\\Files\\Box Sync\\PhD\\1 Research projects\\Chapter2 - PSF\\Code\\7Jan25') 

rawdata<-qza_to_phyloseq(features='dada2_table_2_newIDs.qza', tree = 'testRooted-tree.qza', 
                         taxonomy='Taxonomy_2.qza', metadata = "CH2_metadata_newIDs4.txt")

## rarefy and subset samples
test <- ggrare(rawdata, step = 1000)
test+scale_x_continuous(breaks = seq(0, 20000, 2000))+xlim(c(0,20000))

psf<-rarefy_even_depth(rawdata, rngseed=1, sample.size=7500, replace=F) ## elims samples 427, negS

### removing outlier samples at this step (7/12/2025) - not done in initial analysis.
### will otherwise rerun code from here

## having trouble subsetting from the sample name column so subsetting via a different
## unique data column

# Sample106 = 48E
#Sample95 is already missing? (likely a depleted DNA tube sample)
#Sample142=12D

psfX<-subset_samples(psf,position !='48E')
psfX<-subset_samples(psfX,position !='12D')
psf<-psfX

## Continuing prior analysis as before

psf2<-subset_samples(psf, soil!='neg') ## remove negs
psfLive<-subset_samples(psf2,soil=="Live")
psfLM<-subset_samples(psfLive,species=="Monarda")
psfLP<-subset_samples(psfLive,species=="Plantago")

## sample data as dfs

psf_df <- as(sample_data(psf), "data.frame")
psf2_df <- as(sample_data(psf2), "data.frame")
psfLive_df <- as(sample_data(psfLive), "data.frame")
psfLP_df <- as(sample_data(psfLP), "data.frame")
psfLM_df <- as(sample_data(psfLM), "data.frame")


#calculate overall alpha diversty & merge with metadata (sample IDs come out as rownames - need their own data column to merge)

Alphas<-estimate_richness(psf, split=TRUE, measures=NULL)
Alphas$sample<-rownames(Alphas)
psf_df$sample<-rownames(psf_df) 
AlphasMeta<-left_join(psf_df, Alphas,by="sample")


p1<- ggplot(AlphasMeta, aes(x=soil,fill=soil, y=Shannon)) + 
  geom_boxplot(notch=TRUE)+theme_bw()+ 
  geom_jitter(color="black", size=0.4, alpha=0.9)

## linear models + anovas of linear models

## diversity Shannon

AlphasMeta2<-AlphasMeta[AlphasMeta$soil!='neg',]#remove negs

m1<-lm(Shannon~soil,data=AlphasMeta2)
m2<-lm(Shannon~species,data=AlphasMeta2)
summary(m1)
summary(m2)

car::Anova(m1,type=2) ## significant live vs. sterile difference
car::Anova(m2,type=2) ## no significant effect of plant sp. on alpha diversity



## create distance matrices
dmBC <- phyloseq::distance(psf2, method ="bray")
dmJC <- phyloseq::distance(psf2, method ="jaccard")

dmBClive <- phyloseq::distance(psfLive, method ="bray")
dmBC_LP <- phyloseq::distance(psfLP, method ="bray")
dmBC_LM <- phyloseq::distance(psfLM, method ="bray")


############# PERMANOVA
## writing explicitly with : for interaction tests. prior Ch1 code writing was redundant with *, which specifies e.g. soil + water AND
## soil:water

## all samples
#x1_all<-adonis2(dmBC~soil*species*soil_source+soil*species+soil_source*soil+soil_source*species,psf2_df,by="margin") # three-way interaction!
q1<-adonis2(dmBC~soil*species*soil_source,psf2_df,by="margin")
# x1_2way<-adonis2(dmBC~soil*species+soil_source*soil+soil_source*species,psf2_df,by="margin") #
# 
# x1<- adonis2(dmBC~soil+species+soil_source+species:soil_source,psf2_df,by="margin") #soil:soil_source interaction
# x1_2<- adonis2(dmBC~species+soil:soil_source,psf2_df,by="margin") 
## looking at all samples (live and sterile), soil:soil_source interaction, no interactions with species.
## species alone is significant. Additionally, 3-way soil:species:soil_source interaction

########################
#LIVE SAMPLES

x2<- adonis2(dmBClive~species*soil_source,psfLive_df,by="margin")



#############################################################
## TAXONOMIC COUNTS AND ANALYSIS #############################

## Plots

#### extracting scores for better plotting 
#############################################################################

OTUs = as(otu_table(psf2), "matrix") # Extract abundance matrix from the phyloseq object
if(taxa_are_rows(psf2)){OTUs <- t(OTUs)} # transpose if necessary
fullMDS=metaMDS(OTUs,k=3) # works with bray

## need to add metadata. this tutorial recommends extracting spp scores, adding metadata, and plotting in ggplot
data.scores = as.data.frame(scores(fullMDS)$sites) #extract NMDS scores (x and y coordinates)
data.scores$sample<-rownames(data.scores) ## merge with metadata (rownames are sample IDS)
psf2_df <- as(sample_data(psf2), "data.frame")

## clean up alphas dfs so they can be merged with sample metadata - the phyloseq f(X) output changes the rownames -
## and merge alpha metrics with sample metadata
psf2_df$sample<-rownames(psf2_df) 
baseNMDS<-left_join(psf2_df, data.scores,by="sample") ## all samples


########### the same for just live samples
OTUsL = as(otu_table(psfLive), "matrix")
if(taxa_are_rows(psfLive)){OTUsL <- t(OTUsL)}
liveMDS=metaMDS(OTUsL,k=3) 
data.scoresL = as.data.frame(scores(liveMDS)$sites)
data.scoresL$sample<-rownames(data.scoresL)
psfLive_df <- as(sample_data(psfLive), "data.frame")
psfLive_df$sample<-rownames(psfLive_df) 
liveNMDS<-left_join(psfLive_df, data.scoresL,by="sample")

liveMDS4=metaMDS(OTUsL,k=4) ## k = 3 has relatively high stress (~0.2); will try 4
data.scoresL4 = as.data.frame(scores(liveMDS4)$sites)
data.scoresL4$sample<-rownames(data.scoresL4)
# psfLive_df <- as(sample_data(psfLive), "data.frame")
# psfLive_df$sample<-rownames(psfLive_df) 
liveNMDS4<-left_join(psfLive_df, data.scoresL4,by="sample")



#### just Plantago / Monarda (live samples)
# 
# OTUsM = as(otu_table(psfLM), "matrix")
# if(taxa_are_rows(psfLM)){OTUsM <- t(OTUsM)}
# monMDS=metaMDS(OTUsM,k=3)
# data.scoresM = as.data.frame(scores(monMDS)$sites)
# data.scoresM$sample<-rownames(data.scoresM)
# psfM_df <- as(sample_data(psfLM), "data.frame")
# psfM_df$sample<-rownames(psfM_df)
# monNMDS<-left_join(psfM_df, data.scoresM,by="sample")
# # 
# OTUsP = as(otu_table(psfLP), "matrix")
# if(taxa_are_rows(psfLP)){OTUsP <- t(OTUsP)}
# plaMDS=metaMDS(OTUsP,k=3)
# data.scoresP = as.data.frame(scores(plaMDS)$sites)
# data.scoresP$sample<-rownames(data.scoresP)
# psfP_df <- as(sample_data(psfLP), "data.frame")
# psfP_df$sample<-rownames(psfP_df)
# plaNMDS<-left_join(psfP_df, data.scoresP,by="sample")

####
#make metadata column for combo species + soil - for symbol changing for full plot

baseNMDS$species_soil<-paste(baseNMDS$species,baseNMDS$soil)
# plaNMDS$species_soil<-paste(plaNMDS$species,plaNMDS$soil)
# monNMDS$species_soil<-paste(monNMDS$species,monNMDS$soil)

##########################################

#new NMDS code
#be sure to relevel first or colors/symbols will be incorrect

level_order<-c("Monarda_Monarda","Plantago_Monarda","Monarda_Plantago","Plantago_Plantago")
liveNMDS$species_soil_source <- factor(liveNMDS$species_soil_source, levels = level_order) # to change levels
liveNMDS4$species_soil_source <- factor(liveNMDS4$species_soil_source, levels = level_order) # to change levels

soillabs=c("Monarda (Monarda soil)","Plantago (Monarda soil)","Monarda (Plantago soil)","Plantago (Plantago soil)")
fulllabs=c("Monarda (Monarda soil)","Plantago (Monarda soil)","Monarda (Plantago soil)","Plantago (Plantago soil)","Monarda (Monarda soil)","Plantago (Monarda soil)","Monarda (Plantago soil)","Plantago (Plantago soil)")
# baseNMDS$species_soil_source <- factor(baseNMDS$species_soil_source, levels = level_order) # to change levels
baseNMDS$soil_species_soil_source <- paste(baseNMDS$soil,baseNMDS$species_soil_source)
level_order_full<-c("Live Monarda_Monarda","Live Plantago_Monarda","Live Monarda_Plantago","Live Plantago_Plantago","Sterile Monarda_Monarda","Sterile Plantago_Monarda","Sterile Monarda_Plantago","Sterile Plantago_Plantago")
baseNMDS$soil_species_soil_source <- factor(baseNMDS$soil_species_soil_source, levels = level_order_full) # to change levels



pFull = ggplot(baseNMDS, aes(x = NMDS1, y = NMDS2,colour=soil_species_soil_source)) + 
  geom_point(size = 5, stroke=1, aes( shape = soil_species_soil_source, colour = soil_species_soil_source))+
  theme_bw()  +scale_shape_manual(labels=fulllabs,values=c(16,17,16,17,1,2,1,2),name="Plant (soil origin)")+ggtitle("All samples (Fungi)")+
  scale_colour_manual(labels=fulllabs,name="Plant (soil origin)",values=c("#002594","#E0B2CD","#00A86B","#D2C500","#002594","#E0B2CD","#00A86B","#D2C500"))

#k = 3
pLive = ggplot(liveNMDS, aes(x = NMDS1, y = NMDS2,colour=species_soil_source,shape=species_soil_source)) + 
  geom_point(size = 5, stroke=1, aes(shape = species_soil_source, colour = species_soil_source))+
  theme_bw()  + stat_ellipse()+scale_shape_manual(values=c(16,17,16,17),name="Plant (soil origin)",labels=soillabs)+ggtitle("Live samples (Fungi)")+
  scale_colour_manual(labels=soillabs,name="Plant (soil origin)", values=c("#002594","#E0B2CD","#00A86B","#D2C500"))

#k=4
pLive4 = ggplot(liveNMDS4, aes(x = NMDS1, y = NMDS2,colour=species_soil_source,shape=species_soil_source)) + 
  geom_point(size = 5, stroke=1, aes(shape = species_soil_source, colour = species_soil_source))+
  theme_bw()  + stat_ellipse()+scale_shape_manual(values=c(16,17,16,17),name="Plant (soil origin)",labels=soillabs)+ggtitle("Live samples (Fungi)")+
  scale_colour_manual(labels=soillabs,name="Plant (soil origin)", values=c("#002594","#E0B2CD","#00A86B","#D2C500"))


# Create text for on-plot labels, then add them with annotation_custom()
# Bray <- grobTree(textGrob("Bray-Curtis", x=0.65,  y=0.95, hjust=0,
#                           gp=gpar(col="black", fontsize=12)))

## 12 July 25 - on-plot NMDS stats remain the same, changed positioning slightly

statL <- grobTree(textGrob("stress < 0.15               k=4", x=0.05,  y=0.05, hjust=0,
                           gp=gpar(col="black", fontsize=12)))

statF <- grobTree(textGrob("stress < 0.15   k=3", x=0.1,  y=0.05, hjust=0,
                           gp=gpar(col="black", fontsize=12)))

plotLive<- pLive4 + annotation_custom(statL)+ 
  theme(plot.title = element_text(size = 14))
plotFull<-pFull + annotation_custom(statF)+ 
  theme(plot.title = element_text(size = 14))


# # Just PlaLan & MonFis (live)
# #################################################
# #################################################
# 
# plantlabs=c("Monarda soil","Plantago soil")
# PlaLan<-ggplot(plaNMDS, aes(x = NMDS1, y = NMDS2,colour=species_soil_source)) + 
#   geom_point(size = 5, stroke=1, aes( shape = species_soil_source, colour = species_soil_source))+
#   theme_bw()+scale_shape_manual(values=c(17,17),labels=plantlabs,name="Soil origin")+ggtitle("Plantago Plant (Fungi)")+
#   scale_colour_manual(name="Soil origin",values=c("#E0B2CD","#D2C500"),labels=plantlabs)+stat_ellipse()
# 
# MonFis<-ggplot(monNMDS, aes(x = NMDS1, y = NMDS2,colour=species_soil_source)) + 
#   geom_point(size = 5, stroke=1, aes( shape = species_soil_source, colour = species_soil_source))+
#   theme_bw()  +scale_shape_manual(values=c(16,16),labels=plantlabs,name="Soil origin")+ggtitle("Monarda Plant (Fungi)")+
#   scale_colour_manual(name="Soil origin",values=c("#002594","#00A86B"),labels=plantlabs)+stat_ellipse()
# 
# statP <- grobTree(textGrob("stress < 0.14   k=3", x=0.53,  y=0.9, hjust=0,
#                            gp=gpar(col="black", fontsize=11)))
# 
# statM <- grobTree(textGrob("stress < 0.16                      k=3", x=0.1,  y=0.05, hjust=0,
#                            gp=gpar(col="black", fontsize=11)))
# 
# plotP<-PlaLan+annotation_custom(statP)
# plotM<-MonFis+annotation_custom(statM)

######################################
## stacked barplots
# https://github.com/david-barnett/microViz/blob/main/README.Rmd
# installing microViz for simpler (?) stacked barplots
library(microViz)

## comparing thetax_fix()ed tax table with the ps table, there is some redundancy in Incertae_sedis
## and NA classifications. Tax_fix codes these as different IDs. A simple solution (?) seems
## to be changing all cells containing Incertae_sedis in the tax table to NA, then letting tax_fix fix them.
## I will need to manually check what happens to taxa identified to genus or lower with uncertain
## higher tax placement but I think this will not be affected


psL<-psfLive #make intermediate ps for editing tax table
psLtax<-tax_table(psL) # extract tax table
df<-as.data.frame(psLtax)

#there is probably a way to apply this to all columns but ~if it ain't broke~

df$Phylum<-str_replace_all(df$Phylum,"Incertae_sedis",NA_character_)
df$Class<-str_replace_all(df$Class,"Incertae_sedis",NA_character_)
df$Order<-str_replace_all(df$Order,"Incertae_sedis",NA_character_)
df$Family<-str_replace_all(df$Family,"Incertae_sedis",NA_character_)
df$Genus<-str_replace_all(df$Genus,"Incertae_sedis",NA_character_)
df$Species<-str_replace_all(df$Species,"Incertae_sedis",NA_character_)

## attempt to recreate ps object with new tax table
taxa_names <- rownames(df)  # Store row names
df_mat <- as.matrix(df)
rownames(df_mat) <- taxa_names  # Ensure row names are intact
newtable <- tax_table(df_mat)# Convert back to a phyloseq tax_table
tax_table(psL) <- newtable # Replace the tax_table in the phyloseq object

### make barplots with psL as phyloseq base
################################################

# custom color palette for 25 most common families

myPal <- tax_palette(
  data = bp, rank = "Family", n = 25, pal = "brewerPlus",
  add = c(Other = "grey88")
)

#adjust colors to custom palette

myPal["Fungi Kingdom"] <- "#A898B0FF"
myPal["Aspergillaceae"] <-"#E8E0F8FF"
myPal["Spizellomycetaceae"] <- "#F8D860FF"
myPal["Cladosporiaceae"] <-"#D8B840FF"
myPal["Glomeraceae"] <-  "#408890FF"
myPal["Rozellomycota Phylum"] <- "#000000FF"
myPal["Plectosphaerellaceae"] <- "#E31A1C"
myPal["Pleosporaceae"] <- "#F8F8F8FF"
myPal["Ascobolaceae"] <-"#9fccd4"
myPal["Nectriaceae"] <- "#332288FF"
myPal["Basidiomycota Phylum"] <- "#c26364"
myPal["Chytridiomycota Phylum"] <-  "#b07ccc"
myPal["Auriculariales Order"] <-"#607039"
myPal["Hypocreales Order"] <-"#a8c462"
myPal["Ceratobasidiaceae"] <-"#d65b1e"

bp<-psL %>% 
  tax_fix() %>%
  phyloseq_validate()

comp_barplot(bp,tax_level='Family',n_taxa=15,palette=myPal) + facet_wrap(~species_soil_source,scales="free_y")+
  coord_flip()+labs(x=NULL,y=NULL)+theme(axis.text.y=element_blank())+scale_color_manual(labels=c("1","2","3"))

comp_barplot(bp,tax_level='Phylum',n_taxa=6) + facet_wrap(~species,scales="free_x") + theme(axis.text.x=element_blank())

comp_barplot(bp,tax_level='Phylum',n_taxa=8) + facet_wrap(~species+soil_source,scales="free_x") + theme(axis.text.x=element_blank())

## mean abundance barplots 
bp %>%
  phyloseq::merge_samples(group = "species_soil_source") %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 6, bar_width = 0.8) +
  coord_flip() + labs(x = NULL, y = NULL)

bp %>%
  phyloseq::merge_samples(group = "species_soil_source") %>%
  comp_barplot(tax_level = "Family", n_taxa = 15, bar_width = 0.8) +
  coord_flip() + labs(x = NULL, y = NULL)


### 

bp %>%
  comp_barplot(tax_level = "Family", n_taxa = 15, bar_width = 0.8) +
  labs(x = NULL, y = NULL)






