## ANCOM
#load main data workspace. ANCOM-BC requires unrarefied data for appropriate
## bias correction so will need to work from original ps object ("rawdata")

BiocManager::install("ANCOMBC")
#install.packages("lme4", type = "source") 
## matrix or lme4 needed to be changed for compatibility when I was using R 4.3.1.
## Not an issue in current R 4.5 version

library(tidyverse)
library(ANCOMBC)
library(lme4)
library(phyloseq)
library(microViz)
library(ggplot2)
library(FUNGuildR)

library(ggpubr)
library(stringr)
library(scales)
##
#setwd('C:\\Files\\Box Sync\\PhD\\1 Research projects\\Chapter2 - PSF\\Code\\7Jan25') 

setwd('C:\\Files\\PhD\\1 Research projects\\Chapter2 - PSF\\Code\\7Jan25')

## "rawdata" ps is straight-from-qiime2 object
## remove negs

rawANCOM<-subset_samples(rawdata,species!="neg")
rl<-subset_samples(rawANCOM,soil=="Live")
ps<-rl #make intermediate ps for editing tax table
pstax<-tax_table(ps) # extract tax table
df2<-as.data.frame(pstax)

#there is probably a way to apply this to all columns but ~if it ain't broke~
df2$Family<-str_replace_all(df2$Family,"Incertae_sedis",NA_character_)
df2$Phylum<-str_replace_all(df2$Phylum,"Incertae_sedis",NA_character_)
df2$Class<-str_replace_all(df2$Class,"Incertae_sedis",NA_character_)
df2$Order<-str_replace_all(df2$Order,"Incertae_sedis",NA_character_)
#fam
df2$Genus<-str_replace_all(df2$Genus,"Incertae_sedis",NA_character_)
df2$Species<-str_replace_all(df2$Species,"Incertae_sedis",NA_character_)

## recreate ps object with new tax table
taxa_names2 <- rownames(df2)  # Store row names
df_mat2 <- as.matrix(df2)
rownames(df_mat2) <- taxa_names2  # Ensure row names are intact
newtable2 <- tax_table(df_mat2)# Convert back to a phyloseq tax_table
tax_table(rl) <- newtable2 # Replace the tax_table in the phyloseq object

####
rlMonarda<-subset_samples(rl,species=="Monarda")
rlPlantago<-subset_samples(rl,species=="Plantago")

## make aggregated/tax_fix()ed phyloseq object
rl2<-rl %>% 
  tax_fix() %>%
  phyloseq_validate()

## remove low abundance reads (must occur in 20% of samples)

rl3<-tax_filter(rl2,min_prevalence=0.2) #prior troubleshooting suggested
# 20% is a good cutoff

#BH/fdr common p-adjust choice for microbiome data. Others too strict or
#correct for things that don't necessarily need correcting

out<-ancombc2(rl3,fix_formula="species_soil_source",p_adj_method="BH") #ASV
outG<-ancombc2(rl3,fix_formula="species_soil_source",p_adj_method="BH",tax_level="Genus")
outF<-ancombc2(rl3,fix_formula="species_soil_source",p_adj_method="BH",tax_level="Family")


ancomRes<-out$res
ancomResGen<-outG$res
ancomResFam<-outF$res


write.csv(ancomRes, "ancomRes.csv")
write.csv(ancomResGen,"ancomResGen.csv")
write.csv(ancomResFam,"ancomResFam.csv")

#simplify column names
colnames(ancomRes) = gsub("species_soil_source", "", colnames(ancomRes))
colnames(ancomRes) = gsub("\\(Intercept\\)", "i", colnames(ancomRes))
colnames(ancomRes) = gsub("Monarda_", "m", colnames(ancomRes))
colnames(ancomRes) = gsub("Plantago_", "p", colnames(ancomRes))
colnames(ancomRes) = gsub("Monarda", "m", colnames(ancomRes))
colnames(ancomRes) = gsub("Plantago", "p", colnames(ancomRes))

colnames(ancomRes)

colnames(ancomResGen) = gsub("species_soil_source", "", colnames(ancomResGen))
colnames(ancomResGen) = gsub("\\(Intercept\\)", "i", colnames(ancomResGen))
colnames(ancomResGen) = gsub("Monarda_", "m", colnames(ancomResGen))
colnames(ancomResGen) = gsub("Plantago_", "p", colnames(ancomResGen))
colnames(ancomResGen) = gsub("Monarda", "m", colnames(ancomResGen))
colnames(ancomResGen) = gsub("Plantago", "p", colnames(ancomResGen))

colnames(ancomResGen)

colnames(ancomResFam) = gsub("species_soil_source", "", colnames(ancomResFam))
colnames(ancomResFam) = gsub("\\(Intercept\\)", "i", colnames(ancomResFam))
colnames(ancomResFam) = gsub("Monarda_", "m", colnames(ancomResFam))
colnames(ancomResFam) = gsub("Plantago_", "p", colnames(ancomResFam))
colnames(ancomResFam) = gsub("Monarda", "m", colnames(ancomResFam))
colnames(ancomResFam) = gsub("Plantago", "p", colnames(ancomResFam))

colnames(ancomResGen)

res20<-ancomRes #for simplification of prior code reuse - res20 is ASV results
resGen<-ancomResGen
resFam<-ancomResFam


##########################

###### assiign funguild functions to ASV results


## Funguild
tt<-tax_table(rl2) #the tax_fixed verson
tt<-as.data.frame(tt)

## create Taxonomy column formatted for Funguild analysis
tt$Taxonomy <- paste(tt$Kingdom,";",tt$Phylum,";",tt$Class,";",tt$Order,";",tt$Family,";",tt$Genus,";",tt$Species)
tt$ASV=rownames(tt)
## Funguild~!
tt2<-funguild_assign(tt,db=get_funguild_db(),tax_col="Taxonomy")
head(tt2)


######## merge funguild data with ancom results
colnames(res20)[1]<-"ASV"
res20tax=left_join(res20,tt2,by="ASV")


#######################################################
#Create CI columns for custom inclusion criteria 

#check distribution of W stat
par(mfrow=c(2,2))
hist(res20$W_mp)
hist(res20$W_pp)
hist(res20$W_pm)


## W for significant results is usually greater than 4. for one 3.5. Will see how many taxa
## more "significant taxa" we get if we use p<0.05 & W > 3 and CI doesn't overlap zero

# create CI cols
#Lower Bound=lfc−1.96×se
#Upper Bound=lfc+1.96×se

res20tax$ciL_mp=res20tax$lfc_mp-(1.96*res20tax$se_mp)
res20tax$ciH_mp=res20tax$lfc_mp+(1.96*res20tax$se_mp)
res20tax$neg_lb_mp = (res20tax$ciL_mp < 0 & res20tax$ciH_mp > 0)

res20tax$ciL_pp=res20tax$lfc_pp-(1.96*res20tax$se_pp)
res20tax$ciH_pp=res20tax$lfc_pp+(1.96*res20tax$se_pp)
res20tax$neg_lb_pp = (res20tax$ciL_pp < 0 & res20tax$ciH_pp > 0)

res20tax$ciL_pm=res20tax$lfc_pm-(1.96*res20tax$se_pm)
res20tax$ciH_pm=res20tax$lfc_pm+(1.96*res20tax$se_pm)
res20tax$neg_lb_pm = (res20tax$ciL_pm < 0 & res20tax$ciH_pm > 0)


#Genera

resGen$ciL_mp=resGen$lfc_mp-(1.96*resGen$se_mp)
resGen$ciH_mp=resGen$lfc_mp+(1.96*resGen$se_mp)
resGen$neg_lb_mp = (resGen$ciL_mp < 0 & resGen$ciH_mp > 0)

resGen$ciL_pp=resGen$lfc_pp-(1.96*resGen$se_pp)
resGen$ciH_pp=resGen$lfc_pp+(1.96*resGen$se_pp)
resGen$neg_lb_pp = (resGen$ciL_pp < 0 & resGen$ciH_pp > 0)

resGen$ciL_pm=resGen$lfc_pm-(1.96*resGen$se_pm)
resGen$ciH_pm=resGen$lfc_pm+(1.96*resGen$se_pm)
resGen$neg_lb_pm = (resGen$ciL_pm < 0 & resGen$ciH_pm > 0)

####### Families

resFam$ciL_mp=resFam$lfc_mp-(1.96*resFam$se_mp)
resFam$ciH_mp=resFam$lfc_mp+(1.96*resFam$se_mp)
resFam$neg_lb_mp = (resFam$ciL_mp < 0 & resFam$ciH_mp > 0)

resFam$ciL_pp=resFam$lfc_pp-(1.96*resFam$se_pp)
resFam$ciH_pp=resFam$lfc_pp+(1.96*resFam$se_pp)
resFam$neg_lb_pp = (resFam$ciL_pp < 0 & resFam$ciH_pp > 0)

resFam$ciL_pm=resFam$lfc_pm-(1.96*resFam$se_pm)
resFam$ciH_pm=resFam$lfc_pm+(1.96*resFam$se_pm)
resFam$neg_lb_pm = (resFam$ciL_pm < 0 & resFam$ciH_pm > 0)


##############################################################
## filter/identify results such that

#q <0.2
#W >2
#CI doesn't overlap zero 

########################################## ASV
res20tax <- res20tax %>%
  mutate(
    pass1 = (W_mp > 2 & q_mp < 0.2 & neg_lb_mp == FALSE)
  )
table(res20tax$pass1) #4

### pp
res20tax <- res20tax %>%
  mutate(
    pass3 = (W_pp > 2 & q_pp < 0.2 & neg_lb_pp == FALSE)
  )
table(res20tax$pass3) #also 4

### pm
res20tax <- res20tax %>%
  mutate(
    pass5 = (W_pm > 2 & q_pm < 0.2 & neg_lb_pm == FALSE)
  )
table(res20tax$pass5) #2

######################################################### Genus
resGen <- resGen %>%
  mutate(
    pass1 = (W_mp > 2 & q_mp < 0.2 & neg_lb_mp == FALSE)
  )
table(resGen$pass1) #8

resGen <- resGen %>%
  mutate(
    pass3 = (W_pp > 2 & q_pp < 0.2 & neg_lb_pp == FALSE)
  )
table(resGen$pass3) #8

resGen <- resGen %>%
  mutate(
    pass5 = (W_pm > 2 & q_pm < 0.2 & neg_lb_pm == FALSE)
  )
table(resGen$pass5) #7


## family
#########################################################
resFam <- resFam %>%
  mutate(
    pass1 = (W_mp > 2 & q_mp < 0.2 & neg_lb_mp == FALSE)
  )
table(resFam$pass1) #7

### pp
resFam <- resFam %>%
  mutate(
    pass3 = (W_pp > 2 & q_pp < 0.2 & neg_lb_pp == FALSE)
  )
table(resFam$pass3) #3

### pm
resFam <- resFam %>%
  mutate(
    pass5 = (W_pm > 2 & q_pm < 0.2 & neg_lb_pm == FALSE)
  )
table(resFam$pass5) #5


################ FILTER "SIGNIFICANT" RESULTS FOR RESULTS TABLE / Prep for plotting
######################################################################################

## make Positive/Negative 

###### make direction (pos or neg LFC) column

## make direction column
res20tax <- res20tax %>%
  mutate(direction_pm = ifelse(lfc_pm >= 0, "Positive", "Negative"))%>%
  mutate(direction_mp=ifelse(lfc_mp>=0,"Positive","Negative"))%>%
  mutate(direction_pp=ifelse(lfc_pp>=0,"Positive","Negative"))


########################################
##FILTER "SIGNIFICANT" RESULTS 

sigMP<-res20tax %>%
  filter(diff_mp == "TRUE" | pass1 == "TRUE") 

sigMP<-sigMP %>%
  arrange(desc(lfc_mp)) %>%  # Sorts by log-fold change (largest to smallest)
  mutate(ASV = factor(ASV, levels = ASV)) 

sigPP<-res20tax %>%
  filter(diff_pp == "TRUE" | pass3 == "TRUE") %>%
  filter(abs(lfc_pp)>0.5)

sigPP<-sigPP %>%
  arrange(desc(lfc_pp)) %>%  # Sorts by log-fold change (largest to smallest)
  mutate(ASV = factor(ASV, levels = ASV)) 

sigPM<-res20tax %>%
  filter(diff_pm == "TRUE" | pass5 == "TRUE") %>%
  filter(abs(lfc_pm)>0.5)

sigPM<-sigPM %>%
  arrange(desc(lfc_pm)) %>%  # Sorts by log-fold change (largest to smallest)
  mutate(ASV = factor(ASV, levels = ASV)) 



##########################
## rbind & make comparison column

sigMP$comp<-"MP"
sigPP$comp<-"PP"
sigPM$comp<-"PM"
sigASV<-rbind(sigMP,sigPP,sigPM)

write.csv(sigASV,"sigASV.csv")

##... possible the above isn't including everything?

#return to the below after ASV results confirmed
####################################################


####### PLOTS

my_theme<-theme(plot.title = element_text(hjust = 0.5), #fromPlotsofPlots0000000000
                panel.grid.minor.y = element_blank(),
                plot.margin = margin(5, 10, 5, 5, "pt"),
                axis.text.x = element_text(angle=45, hjust = 1,size=8),axis.title.x = element_blank(),legend.position="none")


p4<-sigMP %>%
  ggplot(aes(x = ASV, y = lfc_mp,fill=direction_mp)) + 
  scale_fill_manual(values = c("Positive" = "#00A86B", "Negative" = "#002594"))+
  ylab("Log fold change")+ggtitle("MP vs MM")+
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(ymin = lfc_mp - se_mp, ymax = lfc_mp + se_mp), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),legend.position="none")+
  scale_x_discrete(labels = word(sigMP$Species,1))+my_theme



p5<-sigPP %>%
  ggplot(aes(x = ASV, y = lfc_pp,fill=direction_pp)) + 
  scale_fill_manual(values = c("Positive" = "#D2C500", "Negative" = "#002594"))+
  ylab("Log fold change")+ggtitle("PP vs MM")+
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(ymin = lfc_pp - se_pp, ymax = lfc_pp + se_pp), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),legend.position="none")+
  scale_x_discrete(labels = word(sigPP$Species,1))+my_theme


p6<-sigPM %>%
  ggplot(aes(x = ASV, y = lfc_pm,fill=direction_pm)) + 
  scale_fill_manual(values = c("Positive" = "#E0B2CD", "Negative" = "#002594"))+
  ylab("Log fold change")+ggtitle("PM vs MM")+
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(ymin = lfc_pm - se_pm, ymax = lfc_pm + se_pm), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),legend.position="none")+
  scale_x_discrete(labels = word(sigPM$Species,1))+my_theme



### the above but with PP as reference group
#### PP as reference group!!!!! Just ASV for now
############################################################################################
###########################################################################################
rl4<-rl3 #new ps object for releveling

levels(sample_data(rl4)$species_soil_source)
sample_data(rl4)$species_soil_source <- relevel(
  as.factor(sample_data(rl4)$species_soil_source), 
  ref = "Plantago_Plantago"
)

rl4_out <- ancombc2(rl4,fix_formula="species_soil_source",p_adj_method="BH") ## ASV level = default
res<-rl4_out$res

colnames(res) = gsub("species_soil_source", "", colnames(res))
colnames(res) = gsub("\\(Intercept\\)", "i", colnames(res))
colnames(res) = gsub("Monarda_", "m", colnames(res))
colnames(res) = gsub("Plantago_", "p", colnames(res))
colnames(res) = gsub("Monarda", "m", colnames(res))
colnames(res) = gsub("Plantago", "p", colnames(res))

colnames(res)
colnames(res)[1]<-"ASV"

##merge with taxonomy and funguild data
restaxPP=left_join(res,tt2,by="ASV")

# make CI columns

#Upper Bound=lfc+1.96×se


restaxPP$ciL_mm=restaxPP$lfc_mm-(1.96*restaxPP$se_mm)
restaxPP$ciH_mm=restaxPP$lfc_mm+(1.96*restaxPP$se_mm)
restaxPP$neg_lb_mm = (restaxPP$ciL_mm < 0 & restaxPP$ciH_mm > 0)

restaxPP$ciL_mp=restaxPP$lfc_mp-(1.96*restaxPP$se_mp)
restaxPP$ciH_mp=restaxPP$lfc_mp+(1.96*restaxPP$se_mp)
restaxPP$neg_lb_mp = (restaxPP$ciL_mp < 0 & restaxPP$ciH_mp > 0)

restaxPP$ciL_pm=restaxPP$lfc_pm-(1.96*restaxPP$se_pm)
restaxPP$ciH_pm=restaxPP$lfc_pm+(1.96*restaxPP$se_pm)
restaxPP$neg_lb_pm = (restaxPP$ciL_pm < 0 & restaxPP$ciH_pm > 0)

#mm
restaxPP <- restaxPP %>%
  mutate(
    pass1 = (W_mm > 2 & q_mm < 0.2 & neg_lb_mm == FALSE)
  )
table(restaxPP$pass1) #3

### mp
restaxPP <- restaxPP %>%
  mutate(
    pass3 = (W_mp > 2 & q_mp < 0.2 & neg_lb_mp == FALSE)
  )
table(restaxPP$pass3) #3

### pm
restaxPP <- restaxPP %>%
  mutate(
    pass5 = (W_pm > 2 & q_pm < 0.2 & neg_lb_pm == FALSE)
  )
table(restaxPP$pass5) #2


###
## make direction column
restaxPP <- restaxPP %>%
  mutate(direction_mm = ifelse(lfc_mm >= 0, "Positive", "Negative"))%>%
  mutate(direction_mp=ifelse(lfc_mp>=0,"Positive","Negative"))%>%
  mutate(direction_pm=ifelse(lfc_pm>=0,"Positive","Negative"))

####

## filter sig results
sigMM2<-restaxPP %>%
  filter(diff_mm == "TRUE" | pass1 == "TRUE") 

sigMM2<-sigMM2 %>%
  arrange(desc(lfc_mm)) %>%  # Sorts by log-fold change (largest to smallest)
  mutate(ASV = factor(ASV, levels = ASV)) 
########
sigMP2<-restaxPP %>%
  filter(diff_mp == "TRUE" | pass3 == "TRUE") %>%
  filter(abs(lfc_mp)>0.5)

sigMP2<-sigMP2 %>%
  arrange(desc(lfc_mp)) %>%  # Sorts by log-fold change (largest to smallest)
  mutate(ASV = factor(ASV, levels = ASV)) 
##
sigPM2<-restaxPP %>%
  filter(diff_pm == "TRUE" | pass5 == "TRUE") %>%
  filter(abs(lfc_pm)>0.5)

sigPM2<-sigPM2 %>%
  arrange(desc(lfc_pm)) %>%  # Sorts by log-fold change (largest to smallest)
  mutate(ASV = factor(ASV, levels = ASV)) 


## rbind & make comparison column

sigMM2$comp<-"MM"
sigMP2$comp<-"MP"
sigPM2$comp<-"PM"
sigASV2<-rbind(sigMM2,sigMP2,sigPM2)

write.csv(sigASV2,"sigASV2.csv")


#########################################

##Plots

p1<-sigMP2 %>%
  ggplot(aes(x = ASV, y = lfc_mp,fill=direction_mp)) + 
  scale_fill_manual(values = c("Positive" = "#00A86B", "Negative" = "#D2C500"))+
  ylab("Log fold change")+ggtitle("MP vs PP")+
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(ymin = lfc_mp - se_mp, ymax = lfc_mp + se_mp), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),legend.position="none")+
  scale_x_discrete(labels = word(sigMP2$Species,1))+my_theme



p2<-sigMM2 %>%
  ggplot(aes(x = ASV, y = lfc_mm,fill=direction_mm)) + 
  scale_fill_manual(values = c("Positive" = "#002594", "Negative" = "#D2C500"))+
  ylab("Log fold change")+ggtitle("MM vs PP")+
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(ymin = lfc_mm - se_mm, ymax = lfc_mm + se_mm), 
                width = 0.2, position = position_dodge2(0.05), color = "black")+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),legend.position="none")+
  scale_x_discrete(labels = word(sigPP$Species,1))+my_theme


p3<-sigPM2 %>%
  ggplot(aes(x = ASV, y = lfc_pm,fill=direction_pm)) + 
  scale_fill_manual(values = c("Positive" = "#E0B2CD", "Negative" = "#D2C500"))+
  ylab("Log fold change")+ggtitle("PM vs PP")+
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(ymin = lfc_pm - se_pm, ymax = lfc_pm + se_pm), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),legend.position="none")+
  scale_x_discrete(labels = word(sigPM$Species,1))+my_theme


ggarrange(nrow=2,ncol=3,p3,p1,p2,p4,p6,p5) 

