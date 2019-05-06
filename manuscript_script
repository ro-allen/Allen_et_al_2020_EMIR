#Packages ----
library(ggplot2)
library(phyloseq)
library(reshape2)
library(ape)
library(gridExtra)
library(vegan)
library(cowplot)
library(knitr)
library(dada2)
library(vegan)
library(cowplot)

#Data import----
path = "~/geotraces-exp/Geotraces_exp_unzipped"
list.files(path)

#Bioconductor pipeline----
##Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

##Extract sample names. Filenames have format: SAMPLENAME_XXX.fastq
setwd("~/geotraces-exp")
list.files()
sn <- read.csv('geotraces_sample_names.csv', header = FALSE)
sample.names = sn[,2]
sample.names

##Plot quality scores 
qpf = plotQualityProfile(fnFs[1:36]) 
qpr = plotQualityProfile(fnRs[1:36])
ggsave("geo-exp-quality-F.jpeg", qpf, width = 60, height = 40, units = "cm", device = "jpeg")
ggsave("geo-exp-quality-R.jpeg", qpr, width = 60, height = 40, units = "cm", device = "jpeg")

##Filtering and trimming
filt_path = file.path(path, "filtered")
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

##Primers removed at this stage
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,181), trimLeft=c(20,19),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

##Learn and plot errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plot.errF = plotErrors(errF, nominalQ=TRUE) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot.errR = plotErrors(errR, nominalQ=TRUE) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("geo-exp-error-F.jpeg", plot.errF, width = 30, height = 20, units = "cm", device = "jpeg")
ggsave("geo-exp-error-R.jpeg", plot.errR, width = 30, height = 20, units = "cm", device = "jpeg")

##Dereplicate sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

##Merge paired reads
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers) #not sure if this has worked well. 


##Construct sequence table
seqtab = makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

##Excise sequences with appropriate length for target region
seqtab2 = seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)] 
dim(seqtab2)

##Remove chimeras
seqtab2.nochim = removeBimeraDenovo(seqtab2, method = 'consensus', multithread = T, verbose = T)
dim(seqtab2.nochim)

##Fraction of sequences removed as chimeras (~4% total)
sum(seqtab2.nochim)/sum(seqtab2)

##Check pipeline read loss
getN = function(x) sum(getUniques(x))
track2 = cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab2.nochim))
colnames(track2) = c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track2) = sample.names
View(track2)

##Assign taxonomy; SILVA v132
silva.taxa3 = assignTaxonomy(seqtab2.nochim, "~/geotraces-exp/Geotraces_exp_unzipped/silva_nr_v132_train_set.fa", multithread = T)

##Species level assignment 
silva.taxa3 <- addSpecies(silva.taxa3, "~/geotraces-exp/Geotraces_exp_unzipped/silva_species_assignment_v132.fa")

##View taxonomic assignments
silva.taxa3.print = silva.taxa3
rownames(silva.taxa3.print) = NULL
View(silva.taxa3.print)

#Phylogenetic tree construction
source("http://bioconductor.org/biocLite.R")
biocLite("DECIPHER")
install.packages("phangorn")
library(DECIPHER)
library(phangorn)

seqs = getSequences(seqtab2.nochim)
names(seqs) = seqs 
alignment = AlignSeqs(DNAStringSet(seqs), anchor = NA)
?plotBS
?optim.pml

phang.align = phyDat(as(alignment, "matrix"), type = "DNA")
dm = dist.ml(phang.align)
treeNJ = NJ(dm)
fit = pml(treeNJ, data = phang.align)

fitGTR = update(fit, k=4, inv=0.2)
fitGTR = optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace=0))
detach("package:phangorn", unload=TRUE)
quartz()
plot(fitGTR, main = "Neighbor Joining")

##Format sample metadata
sd.geo = read.csv("geo-mapping.csv")
row.names(sd.geo) = sample_names(geo.exp) #sample data row names must align with dada2 rowname outputs
head(sd.geo)
View(sd.geo)
sd.geo = as.data.frame(sd.geo)
sd.geo$sample_treatment = paste(sd.geo$SampleID, " ", sd.geo$Treatment)
sd.geo$treatment_time = paste(sd.geo$Treatment, " ", sd.geo$Time)
View(sd.geo)

##Construct initial phyloseq object
geo.exp = phyloseq(tax_table(silva.taxa3), otu_table(seqtab2.nochim, taxa_are_rows = FALSE), sample_data(sd.geo), phy_tree(fitGTR$tree)) #silva rebuild (Overwriting)

##Rename ASVs (rows)
dim(otu_table(geo.exp))
dim(tax_table(geo.exp))

test.vec = as.vector(1:1937)
test.asv = cbind("asv_", test.vec)
test.asv = as.data.frame(test.asv)
asv.names = paste0(test.asv$V1, test.asv$test.vec)
asv.names = as.data.frame(asv.names)
head(asv.names) 

taxa_names(geo.exp) = asv.names$asv.names
taxa_names(geo.exp)
View(tax_table(geo.exp))



##Root phylogenetic tree 
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}

pick_new_outgroup(phy_tree(geo.exp)) #asv_1783"
rootedTree = ape::root(phy_tree(geo.exp), outgroup="asv_1783", resolve.root=TRUE) 
geo.exp = phyloseq(tax_table(geo.exp), otu_table(geo.exp), sample_data(geo.exp), phy_tree(rootedTree))
phy_tree(geo.exp)

##Check total sequences remaining per sample
rowSums(otu_table(geo.exp))
mean(rowSums(otu_table(geo.exp)))
min(rowSums(otu_table(geo.exp)))
max(rowSums(otu_table(geo.exp)))

##Add sample codes
sample.code = as.vector(c("Technical","CT_D1_A3", "CT_D1_B3", "CT_D1_C3", "OA_D1_A3", "OA_D1_B3", "GH_D1_A3", "GH_D1_B3", "GH_D1_C3", "CT_D5_A3", "CT_D5_B3", "CT_D5_C3", "OA_D5_A3", "OA_D5_B3", "OA_D5_C3", "GH_D5_A3", "GH_D5_B3", "GH_D5_C3", "CT_D1_A1", "CT_D5_A1", "CT_D5_B1", "CT_D5_C1", "OA_D5_A1", "OA_D5_B1", "OA_D5_C1", "GH_D5_A1", "Err", "GH_D5_C1", "CT_D1_B1", "CT_D1_C1", "OA_D1_A1", "OA_D1_B1", "OA_D1_C1", "GH_D1_A1", "GH_D1_B1", "GH_D1_C1"))
sample_data(geo.exp)$sample.code = sample.code

##Restructuring tax_table to include 'best' classification
library(zoo)
library(tibble)
bc.t = t(as.data.frame(tax_table(geo.exp)))
bc.t[bc.t==""] <- NA
bc.fill = na.locf(bc.t, na.rm = TRUE)
t.bc.fill = as.data.frame(t(bc.fill))
head(t.bc.fill)
rnc.bc = rownames_to_column(t.bc.fill, "ASV")

###Creates a column with the best classification and the ASV
rnc.bc$taxa_ASV = paste(rnc.bc$Genus,rnc.bc$ASV) #genus chosen over species for clarity

###Bind this column back onto the original tax_table 
safe.bc = as.data.frame(tax_table(geo.exp))
safe.bc$taxa_ASV = paste(rnc.bc$taxa_ASV)
View(safe.bc)

###Setup object as tax_table
bc.tax = tax_table(safe.bc)
colnames(bc.tax) = colnames(safe.bc)
rownames(bc.tax) = rownames(safe.bc)
View(bc.tax)

###Update phyloseq object with new table
identical(bc.tax[1:1937,1:7], tax_table(geo.exp))
tax_table(geo.exp) = bc.tax
View(tax_table(geo.exp))

##Remove non-bacterial and plastid sequences
gg.x = subset_taxa(geo.exp, Kingdom == "Bacteria")
gg.x = subset_taxa(gg.x, Order != "Chloroplast")
gg.x = subset_taxa(gg.x, Family != "Mitochondria")

##Separate experiments for independent analysis
gg.x.3 = subset_samples(gg.x, Experiment == "E3")
gg.x.3 = filter_taxa(gg.x.3, function(x) sum(x > 5) > (0.1*length(x)), TRUE)

gg.x.1 = subset_samples(gg.x, Experiment == "E1")
#remove RLM8 (phosphate contaminated bottle Day 1)
gg.x.1 = subset_samples(gg.x.1, SampleID != "RLM8")
gg.x.1 = filter_taxa(gg.x.1, function(x) sum(x > 5) > (0.1*length(x)), TRUE)

#Rarefaction----
min(rowSums(otu_table(gg.x.1)))
min(rowSums(otu_table(gg.x.3)))

sample.int(1000, 1) #select random integer for seed; 315
rx1 = rarefy_even_depth(gg.x.1, sample.size = 38691, trimOTUs = TRUE, rngseed = 315)
rx3 = rarefy_even_depth(gg.x.3, sample.size = 9819, trimOTUs = TRUE, rngseed = 315) 

##Rarefaction curves
facet.time = c("T0" = "Day 1", "T5" = "Day 5")

rx1.curve = ggrare(rx1, step = 1000, color = "Treatment", se = FALSE)
rx1.curve = rx1.curve + facet_wrap(~Time,  labeller = as_labeller(facet.time)) + theme_bw()
rx1.curve
quartz()
rx3.curve = ggrare(rx3, step = 1000, color = "Treatment", se = FALSE)
rx3.curve = rx3.curve + facet_wrap(~Time,  labeller = as_labeller(facet.time)) + theme_bw()
rx3.curve

sample_data(gg.x.1)$Treatment = factor(sample_data(gg.x.1)$Treatment, levels = c("Control", "OA", "Greenhouse"))
sample_data(gg.x.3)$Treatment = factor(sample_data(gg.x.3)$Treatment, levels = c("Control", "OA", "Greenhouse"))

gg.1.curve = ggrare(gg.x.1, step = 100, color = "Treatment", se = FALSE)
gg.1.curve = gg.1.curve + facet_wrap(~Time, labeller = as_labeller(facet.time)) + theme_bw() + theme(legend.text=element_text(size=10), legend.title=element_blank(), legend.background = element_rect(color="grey50", size=.3, linetype=1), strip.background = element_blank()) + scale_color_manual(values = c("forestgreen", "violetred","orange"), labels = c("Control", "High CO2", "Greenhouse"))

gg.3.curve = ggrare(gg.x.3, step = 100, color = "Treatment", se = FALSE)
gg.3.curve = gg.3.curve + facet_wrap(~Time, labeller = as_labeller(facet.time)) + theme_bw() + theme(legend.text=element_text(size=10), legend.title=element_blank(), legend.background = element_rect(color="grey50", size=.3, linetype=1), strip.background = element_blank()) + scale_color_manual(values = c("forestgreen", "violetred","orange"), labels = c("Control", "High CO2", "Greenhouse"))

rare.curve.fig = plot_grid(gg.1.curve + theme(axis.title.y = element_text(size = 11), legend.text = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), strip.text = element_text(size = 10)), gg.3.curve + theme(axis.title.y = element_text(size = 11), legend.text = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), strip.text = element_text(size = 10)), ncol = 1, labels = "AUTO") 
save_plot("silva.rare.curve.jpeg", rare.curve.fig, base_width =7, base_height = 7) 

#Beta diversity----
rx1.1 = subset_samples(rx1, Time == "T0")
rx1.5 = subset_samples(rx1, Time == "T5")
rx3.1 = subset_samples(rx3, Time == "T0")
rx3.5 = subset_samples(rx3, Time == "T5")

##PERMANOVA on weighted unifrac 
rx1.1.metadata = as(sample_data(rx1.1), "data.frame")
rx1.1.dist = phyloseq::distance(rx1.1, "wunifrac")
adonis(rx1.1.dist ~ Treatment, data = rx1.1.metadata, permutations = 9999) 

rx1.5.metadata = as(sample_data(rx1.5), "data.frame")
rx1.5.dist = phyloseq::distance(rx1.5, "wunifrac")
adonis(rx1.5.dist ~ Treatment, data = rx1.5.metadata, permutations = 9999) 

rx3.1.metadata = as(sample_data(rx3.1), "data.frame")
rx3.1.dist = phyloseq::distance(rx3.1, "wunifrac")
adonis(rx3.1.dist ~ Treatment, data = rx3.1.metadata, permutations = 9999) 

rx3.5.metadata = as(sample_data(rx3.5), "data.frame")
rx3.5.dist = phyloseq::distance(rx3.5, "wunifrac")
adonis(rx3.5.dist ~ Treatment, data = rx3.5.metadata, permutations = 9999) 

##Betadispersal (Homogeneity of dispersal) analysis
rx1.1.beta = betadisper(rx1.1.dist, rx1.1.metadata$Treatment, sqrt.dist = FALSE) 
permutest(rx1.1.beta, pairwise = T, permutations = 9999)

rx1.5.beta = betadisper(rx1.5.dist, rx1.5.metadata$Treatment, sqrt.dist = FALSE)
permutest(rx1.5.beta, pairwise = T, permutations = 9999)

rx3.1.beta = betadisper(rx3.1.dist, rx3.1.metadata$Treatment, sqrt.dist = FALSE) 
permutest(rx3.1.beta, pairwise = T, permutations = 9999)

rx3.5.beta = betadisper(rx3.5.dist, rx3.5.metadata$Treatment, sqrt.dist = FALSE) 
permutest(rx3.5.beta, pairwise = T, permutations = 9999)

##Ordination
o.rx1.1 = ordinate(rx1.1, method = "PCoA", distance = "wunifrac")
o.rx1.5 = ordinate(rx1.5, method = "PCoA", distance = "wunifrac")
o.rx3.1 = ordinate(rx3.1, method = "PCoA", distance = "wunifrac")
o.rx3.5 = ordinate(rx3.5, method = "PCoA", distance = "wunifrac")

##Plotting results
sample_data(rx1.1)$Treatment = factor(sample_data(rx1.1)$Treatment, levels = c("Control", "OA", "Greenhouse"))
sample_data(rx1.5)$Treatment = factor(sample_data(rx1.5)$Treatment, levels = c("Control", "OA", "Greenhouse"))
sample_data(rx3.1)$Treatment = factor(sample_data(rx3.1)$Treatment, levels = c("Control", "OA", "Greenhouse"))
sample_data(rx3.5)$Treatment = factor(sample_data(rx3.5)$Treatment, levels = c("Control", "OA", "Greenhouse"))

#wunifrac
new.pcoa.1.t0 = plot_ordination(rx1.1, o.rx1.1, type = "sites", color = "Treatment") + theme_bw() + geom_point(size = 4, colour = "black") + geom_point(size = 3.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + scale_color_manual(values = c("white", "grey", "black"))
new.pcoa.1.t5 = plot_ordination(rx1.5, o.rx1.5, type = "sites", color = "Treatment") + theme_bw() + geom_point(size = 4, colour = "black") + geom_point(size = 3.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + scale_color_manual(values = c("white", "grey", "black"))
new.pcoa.3.t0 = plot_ordination(rx3.1, o.rx3.1, type = "sites", color = "Treatment") + theme_bw() + geom_point(size = 4, colour = "black") + geom_point(size = 3.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + scale_color_manual(values = c("white", "grey", "black"))
new.pcoa.3.t5 = plot_ordination(rx3.5, o.rx3.5, type = "sites", color = "Treatment") + theme_bw() + geom_point(size = 4, colour = "black") + geom_point(size = 3.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + scale_color_manual(values = c("white", "grey", "black"))

new.p1 = plot_grid(new.pcoa.1.t0, new.pcoa.1.t5, new.pcoa.3.t0, new.pcoa.3.t5, ncol = 2, labels = "AUTO")
save_plot("silva.new.pcoa.quadplot.png", new.p1, base_width = 7, base_height = 6)
