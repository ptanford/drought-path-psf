##

#Decontam edit RCV

##Based on 200817_DADA2_BM by B. Mahmud.
## Edited by RCV 210202

#R ran interactively /home/bmahmud/R/x86_64-pc-linux-gnu-library/3.6 as the directory
#The default directory in RStudio is /home/bmahmud/R/x86_64-pc-linux-gnu-library/3.5
#Packages installed in R interactively are not being recognized here
#Change the default directory to /home/bmahmud/R/x86_64-pc-linux-gnu-library/3.6
myPaths <- .libPaths() #Create a list with all the preset paths
myPaths <- c(myPaths, "/home/rcvargas/R/x86_64-pc-linux-gnu-library/3.6") #Append the new directory
myPaths <- c(myPaths[3], myPaths[2]) #Make the new directory the default one
.libPaths(myPaths)

##Decontam
#Decontam is used on a phyloseq object, an ASV or OTU table
#BiocManager::install("decontam")

library("dada2")


## DADA2 installation
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.10")

#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2")

## All runs (16) from DOME year 1

path <- "/scratch/gdlab/rcvargas/psf/raw_reads"
list.files(path) # entries 217
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

#Check to see that there are equal number of F & R files
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

rand_samples <- sample(size = 20, 1:length(fnFs)) # grab 20 random samples to plot
fwd_qual_plots <- plotQualityProfile(fnFs[rand_samples])
rev_qual_plots <- plotQualityProfile(fnRs[rand_samples])

fwd_qual_plots #Trimming at 248
rev_qual_plots #Trimming at 160

sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
sample.names
duplicated(sample.names, incomparables = FALSE) #Check for duplicates

#Writting the sample names into a cvs file to make a metadata file
#sample.names <- as.data.frame(sample.names)
#write.csv(sample.names, "/scratch/gdlab/bmahmud/DOME/MetaData.csv", row.names = FALSE)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered2", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered2", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(248,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
plotQualityProfile(filtRs[1:20])
plotQualityProfile(filtFs[1:20])

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab_inspect1 <- seqtab[,nchar(colnames(seqtab)) %in% 240] ##inspecting sizes below 250
seqtab_inspect2 <- seqtab[,nchar(colnames(seqtab)) %in% 255:312] ##inspecting larger sizes
(length(unique(substr(getSequences(seqtab),1,230)))) #unique sequences the first 230nt
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 240:257] ##filtering based on size
table(nchar(getSequences(seqtab2)))

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track) <- sample.names
head(track)

#Assigning taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/scratch/gdlab/rcvargas/psf/silva_nr_v132_train_set.fa", multithread=TRUE)
#taxa <- addSpecies(taxa, "/scratch/gdlab/bmahmud/Test/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

metadata <- read.csv("/scratch/gdlab/rcvargas/psf/210217_psf_meta.csv", header = TRUE)
rownames(metadata) <- metadata$seq_name
#metadata <- metadata[-c(2)]

library("phyloseq")

ps<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
             sample_data(metadata), 
             tax_table(taxa))


#Use decontam if you have negative controls
##Decontam#####
#library("decontam")
#library("ggplot2")

##View phyloseq object
#ps 
##21,148 amplicon sequence variants
##View sample metadata info needed for decontam
##Is there is a column that determines controls?
#head(sample_data(ps))

#inspect library sizes
#df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
#df$LibrarySize <- sample_sums(ps)
#df <- df[order(df$LibrarySize),]
#df$Index <- seq(nrow(df))
#ggplot(data=df, aes(x=Index, y=LibrarySize, color=season_type)) + geom_point()
##The library sizes of the positive samples primarily fall from 15,000 to 40,000 reads, 
#but there are some low-read outliers. The negative control samples have fewer reads as expected.

##Two ways to identify contaminants
#First is by frequency, using DNA concentrations
#Second is by Prevalence using known negative controls.

#We’ll summarize the column with control vs sample as a logical variable,
#with TRUE for control samples, as that is the form required by isContaminant.

#sample_data(ps)$is.neg <- sample_data(ps)$season_type == "Control"
#contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
#table(contamdf.prev$contaminant) ## In my case, True = 175 meaning 175 contaminants
#head(which(contamdf.prev$contaminant))

#Can use a more aggressive classification threshold based on probability from 0.1 to 0.5
#contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
#table(contamdf.prev05$contaminant) ## True is now 415
#head(which(contamdf.prev05$contaminant))

#View number of times taxa were observed in positive samples
# Make phyloseq object of presence-absence in negative controls and true samples
##ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
#ps.pa.neg <- prune_samples(sample_data(ps.pa)$season_type == "Control", ps.pa)
#ps.pa.pos <- prune_samples(sample_data(ps.pa)$season_type != "Control", ps.pa)

# Make data.frame of prevalence in positive and negative samples
#df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                   # contaminant=contamdf.prev$contaminant)
#ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  #xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#Is there a clean split between prevalence in pos and neg?

#Remove contaminants
#ps #before
#ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps)
#ps.noncontam #after removal
#################


##Extract ASV tables
###From Akhil's script 08.21.20_OTUTableScript.R
#Creating OTU Tables
mydf <- t(otu_table(ps))

test<-as.data.frame(tax_table(ps),stringsAsFactors = FALSE)
#test[is.na(test)==TRUE]<-"Unclassified"

#test$Family[is.na(test$Family)==TRUE]<-paste("Unclassified", test$Order, sep="_")

for(i in 1:5760){
  if(is.na(test$Phylum[i])==TRUE){
    test$Phylum[i]<-paste("Unclassified", test$Kingdom[i], sep="_")
  } 
  if(is.na(test$Class[i])==TRUE){
    test$Class[i]<-paste("Unclassified", test$Phylum[i], sep="_")
  } 
  if(is.na(test$Order[i])==TRUE){
    test$Order[i]<-paste("Unclassified", test$Class[i], sep="_")
  } 
  if(is.na(test$Family[i])==TRUE){
    test$Family[i]<-paste("Unclassified", test$Order[i], sep="_")
  }
  if(is.na(test$Genus[i])==TRUE){
    test$Genus[i]<-paste("Unclassified", test$Family[i], sep="_")
  }
}

mergedTaxa<-setNames(
  lapply(colnames(test)[1:6],function(taxlevel){
    na.omit(merge(test[c(taxlevel)],mydf,by="row.names"))
  }),
  colnames(test)[1:6]
)

#6 Different CSV Files for Each Taxonomic Rank
Kingdom<-as.data.frame(t(aggregate(.~Kingdom,data=mergedTaxa$Kingdom[2:ncol(mergedTaxa$Kingdom)],sum)))
Phylum<-as.data.frame(t(aggregate(.~Phylum,data=mergedTaxa$Phylum[2:ncol(mergedTaxa$Phylum)],sum)))
Class<-as.data.frame(t(aggregate(.~Class,data=mergedTaxa$Class[2:ncol(mergedTaxa$Class)],sum)))
Order<-as.data.frame(t(aggregate(.~Order,data=mergedTaxa$Order[2:ncol(mergedTaxa$Order)],sum)))
Family<-as.data.frame(t(aggregate(.~Family,data=mergedTaxa$Family[2:ncol(mergedTaxa$Family)],sum)))
Genus<-as.data.frame(t(aggregate(.~Genus,data=mergedTaxa$Genus[2:ncol(mergedTaxa$Genus)],sum)))

TaxaDFs<-setNames(list(Kingdom,Phylum,Class,Order,Family,Genus),colnames(test)[1:6])

TaxaDFs<-lapply(TaxaDFs,function(mydf){
  colnames(mydf) <- as.character(unlist(mydf[1,]))
  return(mydf[-1, ])
})

#The _data.txt name can be anything. Change the location based on where you want the files to be saved.
lapply(colnames(test)[1:6],function(taxlevel){
  write.table(TaxaDFs[[taxlevel]],file=paste0("/scratch/gdlab/rcvargas/psf/ASV_files/",
                                              taxlevel,"_psf.csv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep="\t")
})



#Usually do analysis in seperate rscript for organization
############Analysis with ASV files#################

######From BM script
##Analyses
## Plotting abundance, diversity, PCoA, Adonis, Maaslin, Tukey HSD, MANOVA

#Plotting the abundance data
abund_table<-read.table("/scratch/gdlab/rcvargas/akhil_dome/ASV_all_runs/Family_dome.csv",row.names=1,check.names=FALSE)
abund_table

setwd("/scratch/gdlab/rcvargas/akhil_dome")
metadata.v2 <- read.csv("210104_dome_meta_all.csv", header = TRUE)
rownames(metadata.v2) <- metadata.v2$seq_name
metadata.v2 <- metadata.v2[-c(2)]
meta_table<-metadata.v2[rownames(abund_table),]

x<-abund_table/rowSums(abund_table)
x<-x[,order(colSums(x, na.rm = TRUE),decreasing=TRUE)]

N<-20
taxa_list<-colnames(x)[1:N]

new_x<-data.frame(x[,colnames(x) %in% taxa_list],Others=rowSums(x[,!colnames(x) %in% taxa_list]))

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i])
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405",
             "#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00");
df2<-merge(df,meta_table,by.x="Sample",by.y="row.names")

df3 <- subset(df2, Group=="C")
df3 <- subset(df3, Season!="V")
df3 <- subset(df3, df3$Sample != "3_19-C0402-S001fc01")
df3 <- subset(df3, df3$Sample != "3_19-C0505-S001fc01")

p <- ggplot(df3, aes(x = Subject.ID, y= Value, fill = Taxa)) + geom_bar(stat = "identity") + 
  facet_wrap(~ Season, scales = "free_x") + scale_fill_manual(values=colours[1:(N+1)]) + 
  theme_bw() + ylab("Proportions") + 
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.spacing = unit(0.3, "lines")) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p
pdf("/home/bmahmud/Projects/DOME/Analysis/200414_DOME_16S_fc_nodup/Plots/Abund-Genus-Group_C.pdf",height=36,width=21)
print(p)
dev.off()

#Plotting average abundance for cow samples
agg_cow <- df3[,c(2,3,10)]
agg_cow <- aggregate(agg_cow[,2], by = list(agg_cow$Taxa, agg_cow$Season), mean)
names(agg_cow)[1] <- "Genera"
names(agg_cow)[2] <- "Season"
names(agg_cow)[3] <- "Value"
agg_cow$Season <- factor(agg_cow$Season, levels = c("P", "S", "F", "W"))
ggplot(agg_cow, aes(Season, Value, fill = Genera)) + geom_bar(stat = "identity") +
  scale_fill_manual(values=colours[1:(N+1)]) + theme_bw()+ylab("Proportions") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+
  theme(panel.spacing = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#Generate and plot abundance data for human samples (!= "C") only
abund_table_human <- merge(meta_table, abund_table, by.x = "row.names", by.y = "row.names")
abund_table_human <- abund_table_human[which(abund_table_human$Group != "C"), names(abund_table_human)]
rownames(abund_table_human) <- abund_table_human$Row.names
abund_table_human <- abund_table_human[-c(1)]

meta_table_human <- metadata.v2[rownames(abund_table_human),]

x_hm<-abund_table_human/rowSums(abund_table_human)
x_hm<-x_hm[,order(colSums(x_hm, na.rm = TRUE),decreasing=TRUE)]

N<-20
taxa_list_hm<-colnames(x_hm)[1:N]

new_x_hm<-data.frame(x_hm[,colnames(x_hm) %in% taxa_list_hm],Others=rowSums(x_hm[,!colnames(x_hm) %in% taxa_list_hm]))

df_hm<-NULL
for (i in 1:dim(new_x_hm)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x_hm),Taxa=rep(colnames(new_x_hm)[i],dim(new_x_hm)[1]),Value=new_x_hm[,i])
  if(i==1){df_hm<-tmp} else {df_hm<-rbind(df_hm,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405",
             "#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00");

df2_hm<-merge(df_hm,meta_table_human,by.x="Sample",by.y="row.names")

df3_hm<-subset(df2_hm, Group == "W")

df3_hm$Season <- factor(df3_hm$Season, levels = c("P", "S", "F", "W"))

p<-ggplot(df3_hm,aes(factor(Season),Value,fill=Taxa))+geom_bar(stat="identity") +facet_wrap(~ Subject.ID, scales = "free_x")
p<-p+scale_fill_manual(values=colours[1:(N+1)])
p<-p+theme_bw()+ylab("Proportions")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p
pdf("/home/bmahmud/Projects/DOME/Analysis/200414_DOME_16S_fc_nodup/Plots/Abund-Genus-Group_W.pdf",height=36,width=21)
print(p)
dev.off()

#Plotting average abundance for human samples
agg_hm <- df2_hm[,c(2,3,7,10)]
agg_hm[,2] <- as.numeric(as.character(agg_hm[,2]))
agg_hm <- aggregate(agg_hm[,2], by = list(agg_hm$Taxa, agg_hm$Group, agg_hm$Season), mean)
names(agg_hm)[1] <- "Genera"
names(agg_hm)[2] <- "Group"
names(agg_hm)[3] <- "Season"
names(agg_hm)[4] <- "Value"
agg_hm <- agg_hm[!(agg_hm$Group=="P"),]
agg_hm <- agg_hm[!(agg_hm$Season=="V"),]
agg_hm$Season <- factor(agg_hm$Season, levels = c("P", "S", "F", "W"))
ggplot(agg_hm, aes(Season, Value, fill = Genera)) + geom_bar(stat = "identity") + 
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values=colours[1:(N+1)]) + theme_bw()+ylab("Proportions") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#Calculating Shannon index for all of the samples
ShannDiv<-as.data.frame(diversity(abund_table, index = "shannon"))
ShannDiv_hm <- merge(ShannDiv, meta_table, by.x = "row.names", by.y = "row.names")
names(ShannDiv_hm)[1] <- "Sample"
names(ShannDiv_hm)[2] <- "ShannonIndex"
ShannDiv_hm <- ShannDiv_hm[!(ShannDiv_hm$Group=="P"),]
ShannDiv_hm <- ShannDiv_hm[!(ShannDiv_hm$Season=="V"),]
ggplot(data = ShannDiv_hm, aes(x = Group, y = ShannonIndex, fill = Group)) + 
  geom_boxplot() + 
  geom_point() + facet_wrap(~ Season) +
  theme_light() + labs(y = "Shannon Index", title = "Shannon Diversity of the fecal 16S samples")

ggplot(data = ShannDiv_hm, aes(x = Group, y = ShannonIndex, fill = Group)) + 
  geom_boxplot() + 
  geom_point() + facet_wrap(~ Season) +
  theme_light() + labs(y = "Shannon Index", title = "Shannon Diversity of the fecal 16S samples") +
  stat_compare_means(label.y = 5.0, label.x = 0.7) +
  stat_compare_means(comparisons = my_comparisons, test = "wilcox.test")

ggplot(data = subset(ShannDiv_hm, Season == "P"), aes(x = Group, y = ShannonIndex, fill = Group)) + 
  geom_boxplot() + 
  geom_point() + facet_wrap(~ Season) +
  theme_light() + labs(y = "Shannon Index", title = "Shannon Diversity of the fecal 16S samples") +
  stat_compare_means(label.y = 5.0, label.x = 0.7) +
  stat_compare_means(comparisons = my_comparisons, test = "wilcox.test", 
                     p.adjust.methods = "bonferroni")

my_comparisons = list( c("C", "D"), c("D", "W"), c("C", "W") )

#Checking for normal distribution of the Shannon Indices
ggplot(ShannDiv_hm, aes(sample = ShannonIndex)) + 
  stat_qq_line() + stat_qq_band(alpha = 0.5) + stat_qq_point() + 
  facet_grid(Season ~ Group, scales = "free") +
  theme_light() + labs(y = "Shannon Diversity Index", x = "Theoretical Distribution",
                       title = "Q-Q normalcy test plots")

ggplot(ShannDiv_hm, aes(sample = ShannonIndex)) + 
  stat_pp_line() + stat_pp_band(alpha = 0.5) + stat_pp_point() + 
  facet_grid(Season ~ Group, scales = "free") +
  theme_light() + labs(y = "Shannon Diversity Index", x = "Theoretical Distribution",
                       title = "P-P normalcy test plots")

ggplot(data = subset(ShannDiv_hm, Group == "C" & Season == "S"), aes(sample = ShannonIndex)) + 
  stat_qq_line(distribution = "norm") + stat_qq_band(alpha = 0.5) + stat_qq_point() + 
  theme_light() + labs(y = "Shannon Diversity Index", x = "Theoretical Quatiles",
                       title = "Q-Q normalcy test plots")

ggplot(data = subset(ShannDiv_hm, Group == "C" & Season == "S"), aes(sample = ShannonIndex)) + 
  stat_pp_line() + stat_pp_band(alpha = 0.5) + stat_pp_point() + 
  theme_light() + labs(y = "Shannon Diversity Index", x = "Theoretical Quatiles",
                       title = "P-P normalcy test plots")

ggplot(ShannDiv_hm, aes(x = ShannonIndex, color = Group, fill = Group)) + 
  theme_light() +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.25, binwidth = 0.05) + 
  facet_grid(Group ~ Season, scales = "free") + 
  geom_density(alpha = 0.6)

write.csv(as.data.frame(ddply(ShannDiv_hm, .(Group, Season), summarise, ntest = shapiro.test(ShannonIndex))),
          "/home/bmahmud/Projects/DOME/Analysis/200414_DOME_16S_fc_nodup/Plots/Shapiro_full.csv",
          row.names = FALSE)

ggplot(data = subset(ShannDiv_hm, Group == "C" & Season == "S"), aes(x = ShannonIndex)) + 
  geom_histogram(aes(y = ..density..), color = "red", fill = "white", alpha = 0.5, binwidth = 0.05) + 
  geom_density(color = "red", fill = "red", alpha = 0.25) +
  theme_light(base_line_size = 0) + labs(x = "Shannon Diversity Index", y = "Density",
                       title = "Distribution of Shannon Indeces",
                       subtitle = "Group = C, Season = S")

#Running Kruskal-Wallis test, followed by Wilcox rank test, and BH correction
kruskal.test(ShannonIndex ~ Group, data = test2)
test2 <- subset(ShannDiv_hm, Season == "F")
pairwise.wilcox.test(test2$ShannonIndex, test2$Group, p.adjust.method = "BH")
wilcox.test(test2$ShannonIndex, test2$Group, p.adjust.methods = "BH")
class(test2$ShannonIndex)

#Calculating Bray-Curtis Distance and Plotting PCoA
abund_table <- abund_table[rowSums(abund_table[,-1])>0,] #Removing the rows with sum = 0
beta_div <- vegdist(abund_table, index = "bray")
beta_div
mds <- pco(beta_div, k = 4)
mds$points
mds_data <- as.data.frame(mds$points)
mds_data
df4 <- merge(mds_data, metadata.v2, by.x = "row.names", by.y = "row.names")
df4 <- df4[!(df4$Group=="P"),] #Remove rows with group P.
df4 <- df4[!(df4$Season=="V"),] #Remove raows with season V.
ggplot(data = df4, aes(x = V1, y = V2)) + theme_grey() + geom_point(aes(color=Group)) +
  theme_linedraw(base_line_size = 0) + 
  labs(x = "PCo1", y = "PCo2", title = "PCoA of fecal 16S samples", subtitle = "Index = Bray")

#Plotting the screeplot
scree <- mds$eig*100/sum(mds$eig)
scree <- list(scree)
scree <- do.call(rbind, scree)
scree <- as.data.frame(scree)
scree <- t(scree)
scree <- as.data.frame(scree)
scree$Component <- row.names(scree)
scree$Component <- as.numeric(gsub("[a-zA-z ]", "", scree$Component))
scree <- scree[which(scree$Component < 10), names(scree)]
scree$Component <- as.character(scree$Component)
ggplot(data=scree, aes()) + geom_bar(aes(x=Component, y=V1), stat = "identity") +
  labs(x = "Principal Coordinate", y = "Variance explained") + theme_light()

ggplot(data=df4, aes(x = V1, y = V2)) + theme_test() + facet_grid(rows = vars(Group)) +
  geom_point(aes(color = factor(Season))) #+ labs(title = "PCoA - by day")

df4_cow <- df4[which(df4$Group == "C"), names(df4)]
ggplot(data=df4_cow, aes(x = V1, y = V2)) + theme_test() + facet_wrap(~ df4_cow$Season) +
  geom_point(aes(color = factor(Season)))

df4_hm <- df4[!(df4$Group=="C"),]
ggplot(data=df4_hm, aes(x = V1, y = V2)) + theme_test() +
  geom_point(aes(color = factor(Group))) +
  stat_ellipse(aes(color = factor(Group))) +
  labs(title = "PCoA w/ Bray-Curtis", subtitle = "p = 0.014")

df4_w <- df4[which(df4$Group == "W"), names(df4)]
ggplot(data=df4_w, aes(x = V1, y = V2)) + theme_test() +
  geom_point(aes(color = factor(Season))) +
  stat_ellipse(aes(color = factor(Season)))

df4_d <- df4[which(df4$Group == "D"), names(df4)]
ggplot(data=df4_d, aes(x = V1, y = V2)) + theme_test() +
  geom_point(aes(color = factor(Season)))

#PERMANOVA on seasonal distribution of the worker data
metadata.v2_w <- metadata.v2[which(metadata.v2$Group == "W"), names(metadata.v2)]
metadata.v2_w <- metadata.v2_w[!(metadata.v2_w$Season == "V"),]
abund_table_w <- subset(abund_table, row.names(abund_table) %in% row.names(metadata.v2_w))

beta_div_w <- vegdist(abund_table_w, index = "bray")
beta_div_w
mds_w <- pco(beta_div_w, k = 4)
mds_w$points
#ordiplot(mds_w, type = "n") 
#points(mds_w$points,cex = 0.75, pch = 10, col = metadata.v2_w$Season)
#legend("topright", legend = levels(metadata.v2_w$Season))
mds_data_w <- as.data.frame(mds_w$points)
ordination_w <- merge(mds_data_w, metadata.v2_w, by.x = "row.names", by.y = "row.names")
ggplot(data=ordination_w, aes(x = V1, y = V2)) + theme_test() +
  geom_point(aes(color = factor(Season))) +
  stat_ellipse(aes(color = factor(Season))) +
  theme_light(base_line_size = 0) +
  labs(x = "PCo1", y = "PCo2", title = "PCoA of farmer fecal 16S samples", 
       subtitle = "Index = Bray")
adonis(abund_table_w ~ Season, data = metadata.v2_w,
       permutations = 99999,
       method = "bray")

#Scree plot for the PCoA of the farmer samples only
scree_w <- mds_w$eig*100/sum(mds_w$eig)
scree_w <- list(scree_w)
scree_w <- do.call(rbind, scree_w)
scree_w <- as.data.frame(scree_w)
scree_w <- t(scree_w)
scree_w <- as.data.frame(scree_w)
scree_w$Component <- row.names(scree_w)
scree_w$Component <- as.numeric(gsub("[a-zA-z ]", "", scree_w$Component))
scree_w <- scree_w[which(scree_w$Component < 10), names(scree_w)]
scree_w$Component <- as.character(scree_w$Component)
ggplot(data=scree_w, aes()) + geom_bar(aes(x=Component, y=V1), stat = "identity") +
  labs(x = "Principal Coordinate", y = "Variance explained") + theme_light()

#PERMANOVA on the distribution of different human groups (W vs. D)
metadata.v2_hm <- metadata.v2[!(metadata.v2$Group == "C"),]
metadata.v2_hm <- metadata.v2_hm[!(metadata.v2_hm$Group == "P"),]
metadata.v2_hm <- metadata.v2_hm[!(metadata.v2_hm$Season == "V"),]
abund_table_hm <- subset(abund_table, row.names(abund_table) %in% row.names(metadata.v2_hm))

beta_div_hm <- vegdist(abund_table_hm, index = "bray")
mds_hm <- pco(beta_div_hm, k=4)
mds_data_hm <- as.data.frame(mds_hm$points)
ordination_hm <- merge(mds_data_hm, metadata.v2_hm, by.x = "row.names", by.y = "row.names")
ggplot(data=ordination_hm, aes(x = V1, y = V2)) + theme_test() +
  geom_point(aes(color = factor(Group))) +
  stat_ellipse(aes(color = factor(Group))) + 
  facet_wrap(~ ordination_hm$Season) + 
  theme_light(base_line_size = 0) +
  labs(x = "PCo1", y = "PCo2", title = "PCoA of human fecal 16S samples", 
       subtitle = "Index = Bray, Wrap = Season")
#ordiplot(mds_hm, type = "n") 
#points(mds_hm$points,cex = 0.75, pch = 10, col = metadata.v2_hm$Group)
adonis(abund_table_hm ~ Group, data = metadata.v2_hm,
       permutations = 999,
       method = "bray")

#Scree plot for the PCoA of the human samples only
scree_hm <- mds_hm$eig*100/sum(mds_hm$eig)
scree_hm <- list(scree_hm)
scree_hm <- do.call(rbind, scree_hm)
scree_hm <- as.data.frame(scree_hm)
scree_hm <- t(scree_hm)
scree_hm <- as.data.frame(scree_hm)
scree_hm$Component <- row.names(scree_hm)
scree_hm$Component <- as.numeric(gsub("[a-zA-z ]", "", scree_hm$Component))
scree_hm <- scree_hm[which(scree_hm$Component < 10), names(scree_hm)]
scree_hm$Component <- as.character(scree_hm$Component)
ggplot(data=scree_hm, aes()) + geom_bar(aes(x=Component, y=V1), stat = "identity") +
  labs(x = "Principal Coordinate", y = "Variance explained") + theme_light()

#Comparing W and D for summer samples only
metadata.v2_hm_sum <- metadata.v2_hm[which(metadata.v2_hm$Season == "S"), names(metadata.v2_hm)]
abund_table_hm_sum <- subset(abund_table, row.names(abund_table) %in% row.names(metadata.v2_hm_sum))
adonis(abund_table_hm_sum ~ Group, data = metadata.v2_hm_sum,
       permutations = 99999,
       method = "bray")

#Comparing W and D for fall samples only
metadata.v2_hm_fal <- metadata.v2_hm[which(metadata.v2_hm$Season == "F"), names(metadata.v2_hm)]
abund_table_hm_fal <- subset(abund_table, row.names(abund_table) %in% row.names(metadata.v2_hm_fal))
adonis(abund_table_hm_fal ~ Group, data = metadata.v2_hm_fal,
       permutations = 99999,
       method = "bray")

#Comparing W and D for spring samples only
metadata.v2_hm_spr <- metadata.v2_hm[which(metadata.v2_hm$Season == "P"), names(metadata.v2_hm)]
abund_table_hm_spr <- subset(abund_table, row.names(abund_table) %in% row.names(metadata.v2_hm_spr))
adonis(abund_table_hm_spr ~ Group, data = metadata.v2_hm_spr,
       permutations = 99999,
       method = "bray")

#p-adjustment for seaasonal PERMANOVA comparisons of D & W samples using bonferroni
p_val_hm_seas <- data.frame("Season" = c("Spring", "Summer", "Fall"),
                            "p_raw" = c(0.04297, 0.1822, 0.0056))
p_val_hm_seas$p_adj_bonf <- p.adjust(p_val_hm_seas$p_raw, method = "bonferroni")
p_val_hm_seas$p_adj_BH <- p.adjust(p_val_hm_seas$p_raw, method = "BH")

##################################
