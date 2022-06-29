#### Figure 1B: Prevalence of C. innocuum in the human gut (based on 16S rRNA surveys)
##### Anna M. Seekatz
##### 4.20.22

##### Files needed:
	- cinnocuum_16S_reps_edited.fasta: fasta files of select 16S rRNA sequences from genomes in this study (n=13), formatted for mothur
	- taxonomy_Ci.only_16S_condensed.txt: taxonomy file, formatted to mothur (matched .fasta file)
	- ciprevclost.ci_16S.wang95.tax.summary: custom classifier results (only directed towards C. innocuum)
	- PrevClost_16S_human_meta.txt: metadata for samples included in comparison
	

##### Analyses to be conducted:
	- create a custom classifier for C. innocuum sequences based on our isolates -->
	- compare how many samples match custom classifier across studies (includes our isolates)
	- 16S rRNA datasets used:
		- MacPherson et al: humans, pre and post-abx (n=138; PRJNA414540)
		- Rao et al: hospitalized patients, with or without abx (n = 82; PRJNA631262)
		- Daquigan et al: ONLY selected healthy individuals (n=211; PRJNA386260)
	- metadata:  PrevClost_16S_meta.txt


####-----------------
	
### Step 1: ran 16S rRNA sequences through mothur (v.1.41.3)

# Used following commands in mothur batch:

```
### mbatch for processing 16S rRNA data for comparison
make.contigs(file=ciprevclost.files, processors=2)
summary.seqs(fasta=ciprevclost.trim.contigs.fasta, processors=2)

## note: these changes were implemented after looking at sequence lengths / quality:
screen.seqs(fasta=ciprevclost.trim.contigs.fasta, group=ciprevclost.contigs.groups, maxambig=0, maxlength=300, processors=2)
unique.seqs(fasta=ciprevclost.trim.contigs.good.fasta)
count.seqs(name=ciprevclost.trim.contigs.good.names, group=ciprevclost.contigs.good.groups)
summary.seqs(count=ciprevclost.trim.contigs.good.count_table, processors=2)
count.groups(count=ciprevclost.trim.contigs.good.count_table)

# alignment to FULL database:
align.seqs(fasta=ciprevclost.trim.contigs.good.unique.fasta, reference=silva.seed_v132.align, processors=2)
summary.seqs(fasta=ciprevclost.trim.contigs.good.unique.align, count=ciprevclost.trim.contigs.good.count_table, processors=2)

# eliminate some random seqs
screen.seqs(fasta=ciprevclost.trim.contigs.good.unique.align, count=ciprevclost.trim.contigs.good.count_table, summary=ciprevclost.trim.contigs.good.unique.summary, start=13862, end=23444, maxhomop=8, processors=2)
#system(mv ciprevclost.trim.contigs.good.unique.good.align ciprevclost.trim.contigs.good.unique.good_s13862_e23444.align)
summary.seqs(fasta=ciprevclost.trim.contigs.good.unique.good.align, count=ciprevclost.trim.contigs.good.good.count_table, processors=2)
count.groups(count=ciprevclost.trim.contigs.good.good.count_table)

# now, continue:
filter.seqs(fasta=ciprevclost.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=2)
system(mv ciprevclost.trim.contigs.good.unique.good.filter.fasta ciprevclost.filtered.fasta)
system(mv ciprevclost.trim.contigs.good.good.count_table ciprevclost.filtered.count_table)
	
# remaining steps
unique.seqs(fasta=ciprevclost.filtered.fasta, count=ciprevclost.filtered.count_table)
pre.cluster(fasta=ciprevclost.filtered.unique.fasta, count=ciprevclost.filtered.unique.count_table, diffs=2, processors=2)
count.groups(count=ciprevclost.filtered.unique.precluster.count_table)
chimera.vsearch(fasta=ciprevclost.filtered.unique.precluster.fasta, count=ciprevclost.filtered.unique.precluster.count_table, dereplicate=t, processors=2)
count.groups(count=ciprevclost.filtered.unique.precluster.denovo.vsearch.pick.count_table)
remove.seqs(fasta=ciprevclost.filtered.unique.precluster.fasta, accnos=ciprevclost.filtered.unique.precluster.denovo.vsearch.accnos)
summary.seqs(fasta=current, count=current, processors=2)
classify.seqs(fasta=ciprevclost.filtered.unique.precluster.pick.fasta, count=ciprevclost.filtered.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset16_022016.rdp.fasta, taxonomy=trainset16_022016.rdp.tax, cutoff=80)
remove.lineage(fasta=ciprevclost.filtered.unique.precluster.pick.fasta, count=ciprevclost.filtered.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=ciprevclost.filtered.unique.precluster.pick.rdp.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
count.seqs(name=current, group=current)
count.groups(count=ciprevclost.filtered.unique.precluster.denovo.vsearch.pick.pick.count_table)

```

####-----------------
	
### Step 2: create custom classifier (using mothur) and align 16S data to this

# Used following commands in mothur:

```
align.seqs(fasta=cinnocuum_16S_reps_edited.fasta, reference=silva.seed_v132.align, processors=2)
	# note: 3 sequences were removed; see cinnocuum_16S_reps_edited.flip.accnos
summary.seqs(fasta=cinnocuum_16S_reps_edited.align)
degap.seqs(fasta=cinnocuum_16S_reps_edited.align)
system(mv ciprevclost.filtered.unique.precluster.pick.pick.only_16S_condensed.wang.tax.summary ciprevclost.ci_16S.wang95.tax.summary)
system(mv ciprevclost.filtered.unique.precluster.pick.pick.only_16S_condensed.wang.taxonomy ciprevclost.ci_16S.wang95.taxonomy)
	# note: similar results retrieved with cutoff = 90 or 80.

# checked classification of selected sequences 
align.seqs(fasta=cinnocuum_16S_reps_edited.align, reference=silva.seed_v132.align, processors=2)
unique.seqs(fasta=cinnocuum_16S_reps_edited.align)
classify.seqs(fasta=cinnocuum_16S_reps_edited.unique.align, name=cinnocuum_16S_reps_edited.names, reference=trainset16_022016.rdp.fasta, taxonomy=trainset16_022016.rdp.tax, cutoff=80)
	# note: all classified to Erysi.

# then, align .fasta file of 16S rRNA sample files to custom classifier files:
classify.seqs(fasta=ciprevclost.filtered.unique.precluster.pick.pick.fasta, count=ciprevclost.filtered.unique.prec
luster.denovo.vsearch.pick.pick.count_table, reference=cinnocuum_16S_reps_edited.ng.fasta, taxonomy=taxonomy_Ci.only_16S_co
ndensed.txt, cutoff=95)

```

####-----------------
	
### Step 3: plot in R


```
### read in and format files
taxo<-read.table(file="ciprevclost.ci_16S.wang95.tax.summary", header=TRUE)

# filter out genus-level assignments only, and assign a phylum:
tax<-taxo
	# let's focus on species level, and select from whatever we want later
tax7<-tax[which(tax$taxlevel==7), ]
tax7[, c("rankID", "taxon", "total")]

# cleanup a bit
subtax6 <- tax7
subtax6<-subset(tax7, select=-c(taxlevel, daughterlevels))
subtax6<-subtax6[order(subtax6$rankID, -subtax6$total), ]
taxmatrix<-subtax6[, c(4:ncol(subtax6))]
rownames(taxmatrix)<-subtax6$taxon

# convert to rel abund
genera<- taxmatrix[, colSums(taxmatrix)>8000,]	# remove samples <8000 just in case
	# get rel. abund fraction:
genmatrix<-as.data.frame(t(genera))
genera.fr<-genmatrix/rowSums(genmatrix)*100
genus.fr<-t(genera.fr)
all.genera<-cbind(subtax6[1:3], genus.fr)
	# just for fun, let's look at means...
cbind(as.character(all.genera$taxon), log10(rowMeans(all.genera[,4:length(all.genera)])))
	# looks ok
#write.table(all.genera, file="prevclost.ci_all.genera.16S.isolates.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


### combine with metadata:
summary(as.factor(all.genera$taxon))		# let's select which ones we want here, based on totals and interest

all.genera[, c("taxon", "total")]
genera2<-all.genera
colors <- c("green4", "grey60")
genera2<-cbind(genera2[2:3], colors, genera2[4:ncol(genera2)])
genera2$taxon <- as.character(genera2$taxon)
genera2$taxon[genera2$taxon == c("innocuum") ] <-"C. innocuum"
genera2$taxon <- as.factor(genera2$taxon)
#write.table(genera2, file="prevclost.Ci_genfrac2p.16S.isolates.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### finally, combine with meta:
meta <- read.table("PrevClost_16S_human_meta.txt", header=T, sep="\t")
genbar<-genera2
#genbar <- read.table("prevclost.Ci_genfrac2p.16S.isolates.txt", header=TRUE, sep="\t", comment.char = "")
rownames(genbar)<-genbar$taxon
	#rm_g<-subset(genbar, select =-c(rankID, taxon, colors, total) )
	rm_g<-subset(genbar, select =-c(taxon, colors, total) )
	barg<-as.data.frame(t(rm_g))
	#barg$other<-100-rowSums(barg)
	#others<-100-colSums(barg)
	barg$sampleID<-rownames(barg)
	#col.gen<-c(as.character(genbar$color), "grey80")
	col.gen<-c(as.character(genbar$color))
	#barg$sampleID<-gsub("X19_", "19_", barg$sampleID)
bar<-merge(meta, barg, by.x=c("seqID"), by.y=c("sampleID"))
#write.table(bar, 'prevclost.ci_genfrac2p.16S.isolates.meta.txt',quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)


### graph it:

### let's make a dotplot of the RA (or, log(RA)) for these species, plotting them as sep figures per Group, but separated for human/mouse:
data <- read.table("prevclost.ci_genfrac2p.16S.isolates.meta.txt", header=TRUE, sep="\t", comment.char = "")
summary(as.factor(data$Host_status))
	## there are a few abx the Daquigan study--let's separate those out too
data$Group <- data$Host_status
data$Group[data$Author %in% c("Daquigan") & data$Abx %in% c("YES")] <-"post-abx"
	## let's use these groups 

# plot it:
# note: couldn't figure out a way to add horizontal points, so will just edit in illustrator
par(mai=c(0.75,0.82,0.42,0.22))
data$Group <- as.factor(data$Group)
group.colors <- hcl.colors(4, palette="Zissou 1")
plot<-plot(C..innocuum ~ Group, data = data, ylab=c("log(RA)"), las = 2,
	xlab="", outline=FALSE, ylim=c(0,2.5), cex=0.8, cex.lab=0.8, cex.axis=0.8)
points(C..innocuum ~ jitter(as.numeric(Group, factor=0)), data = data, bg=group.colors[data$Group], cex=0.8, col="black", pch=21, mgp=c(1,2,2))
	# stats:
kruskal.test(C..innocuum ~ Group, data=data)
pairwise.wilcox.test(data$C..innocuum, data$Group, p.adjust.method = "BH")

# let's add the prevalence on the group labels, as well:
# this looks cool--let's also add how many samples per taxa contain >0 of each taxa
	# calculate prevalence
#lapply(split(data[, c("C..innocuum")], data$Group), function(x) signif(100*sum(x > 0)/length(x), 3))
prev <- unlist(as.character(lapply(split(data[, c("C..innocuum")], data$Group), function(x) signif(100*sum(x > 0)/length(x), 3))))
	# graph, adding prevalence to lavel:
par(mai=c(0.75,0.82,0.42,0.22))
data$Group <- as.factor(data$Group)
group.colors <- hcl.colors(4, palette="BluYl")
plot<-plot(log10(C..innocuum+1) ~ Group, data = data, ylab=c("log(RA)"), type="n", xlab="", outline=FALSE, ylim=c(0,1), xaxt="n", col=NA, cex.lab=0.8, cex.axis=0.8)
points(log10(C..innocuum+1) ~ jitter(as.numeric(Group, factor=0)), data = data, bg=adjustcolor(group.colors[data$Group], alpha=.75), cex=0.8, col=adjustcolor("black", alpha=0.5), pch=21, mgp=c(1,2,2))
plot(log10(C..innocuum+1) ~ Group, data = data, ylab=c("log(RA)"), xaxt="n", col=NA, yaxt="n",
	xlab="", outline=FALSE, cex=0.8, cex.lab=0.8, cex.axis=0.8, add=T)
names<-as.character(levels(as.factor(data$Group)))
text(x =  seq(1,4,by=1), y = par("usr")[3]-0.03, srt = 90, adj = 1, labels = names, xpd = TRUE, cex=0.8)
text(x =  seq(1.3,4.3,by=1), y = par("usr")[3]-0.03, srt = 90, adj = 1, labels = paste(prev, "%", sep=""), xpd = TRUE, cex=0.7)

# add n=?
summary(data$Group)
#     Healthy Hospitalized     post-abx       Sepsis 
#         250           50           96           24 

# mean C. innocuum RA:
tapply(data$C..innocuum, data$Group, mean)
#     Healthy Hospitalized     post-abx       Sepsis 
#  0.16421069   0.08586622   0.40832372   0.46090625 

```