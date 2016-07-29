##Loading libraries

library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

##Reading in Brassica FASTA file

DNA <- readDNAStringSet("ChBrGSP/Brassica_napus_v4.1.chromosomes.fa.gz")
testg <-readDNAStringSet("CHBrGSP/Brassica_napus.annotation_v5.gff3.cds.fa.gz")
#Reading in Brassica gff file

genome<- import("CHBrGSP/Brassica_napus.annotation_v5.gff3.gz")

###Preparing Pdict for Biostrings

sequence <- "ACCGGT,ACCGGT,RAATTY,GGATCC,GCTAGC,GGTCTC,TGTACA,CGGCCG,GAATTC,GAATTC,GATATC,AAGCTT,GGTACC,GGTACC,CAATTG,ACGCGT,CCATGG,GCTAGC,GCTAGC,TCGCGA,ATGCAT,CTGCAG,CGATCG,CAGCTG,GAGCTC,GTCGAC,GTCGAC,AGTACT,ACTAGT,GCATGC,AATATT,CCWWGG,TCTAGA,CTCGAG,GCGGCCGC,GCGGCCGC,TTAATTAA,CCTGCAGG,CACNNNGTG,GGTNACC"
enzyme <- "AgeI-HF®,AgeI-HF® RE-Mix®,ApoI-HF,BamHI-HF®,BmtI-HF®,BsaI-HF®,BsrGI-HF®,EagI-HF®,EcoRI-HF®,EcoRI-HF® RE-Mix®,EcoRV-HF®,HindIII-HF®,KpnI-HF®,KpnI-HF® RE-Mix®,MfeI-HF®,MluI-HF®,NcoI-HF®,NheI-HF®,NheI-HF® RE-Mix®,NruI-HF®,NsiI-HF®,PstI-HF®,PvuI-HF®,PvuII-HF®,SacI-HF®,SalI-HF®,SalI-HF® RE-Mix®,ScaI-HF®,SpeI-HF®,SphI-HF®,SspI-HF®,StyI-HF®,XbaI RE-Mix®,XhoI RE-Mix®,NotI-HF®,NotI-HF® RE-Mix® ,PacI RE-Mix®,SbfI-HF®,DraIII-HF®,BstEII-HF®"
name1 <- strsplit(enzyme, ",")
name2 <- strsplit(sequence, ",")
table <- cbind.data.frame(enzyme = name1, sequence =name2)
colnames(table) <- c("enzyme", "sequence")
pdict <- as.character(table[,2])
names(pdict) <- table[,1]
pdict <- DNAStringSet(pdict)
pdict <- PDict(pdict, tb.start = 2, tb.width = 1)

### 1st attempt at matching the pdict to the DNA Sequences

test <- vmatchPDict(pdict, DNA, fixed = "pattern") ### CANT BE RUN, FUNCTION IN PACKAGE INCOMPLETE NEED TO 
                                              ###RUN WITH AN INDIVIDUAL SUBJECT, SEQUENCE IS 41 PARTS,
                                              ###NEED RUN EACH TIME FOR EACH CHROMOSOME

### Counting occurance of each restriction site in each chromosome.
counting <- lapply(1:41, function(s) {
  countPDict(pdict, DNA[[s]])
})

names(counting)<- names(DNA)

### Putting results in dataframe

count.table <- data.frame(counting, row.names = table[,1])
count.table$totalcuts <- rowSums(count.table)
count.table$roughfraglen <- 1300000000/count.table$totalcuts

### Identifies the locations of each restriction site for each enzyme on each chromsome.
famatching <- lapply(1:41, function(s) {
  matchPDict(pdict, DNA[[s]])
})
names(famatching) <- names(DNA) ## Puts names on each chromosome.


### Identifies the location of each restriction site for each gene in the annotated genome.
annomatching <-lapply(1:length(testg), function(s) {
  matchPDict(pdict,testg[[s]])
})
names(annomatching) <-names(testg)
###when indexing with famatching[[x]][[y]] 
###x is chromosome and y is RE

###Isolates one gene and lists the enzymes that cut it and the locations of the cuts as a test
test <- data.frame(annomatching$GSBRNA2T00121640001)
test$group <- table[test$group,1]
test[,1] <- table[test[,1],1]

###Creates a list of data frames of each gene and the enzymes
###which cut it and where
###
listdf <- list()
for(i in 1:length(annomatching)){
  df <- data.frame(annomatching[i])
  listdf[[i]] <- df
}

names(listdf) <- names(testg)

for(i in 1:length(listdf)){
  listdf[[i]][2] <- NULL
  colnames(listdf[[i]]) = c("enzyme", "start", "end", "width")
  listdf[[i]][,1] <- table[listdf[[i]][,1],1]
}

###
### Scratch work due to problem with trusted bands

test.table <- table[c(-3,-32,-39,-40),]
test.dict <- DNAStringSet(pdict)
test.dict<-test.dict[c(-3, -32,-39,-40)]
test.dict <- PDict(test.dict, tb.start = 1, tb.width = 6)
test.dict

test.counting <- lapply(1:41, function(s) {
  countPDict(test.dict, DNA[[s]])
})
names(test.counting) <- names(DNA)
test.count.table <- data.frame(test.counting, row.names = test.table[,1])
test.count.table$totalcuts <- rowSums(test.count.table)
test.count.table$roughfraglen <- 1300000000/test.count.table$totalcuts

test.annomatching <-lapply(1:length(testg), function(s) {
  matchPDict(test.dict,testg[[s]])
})
names(test.annomatching) <-names(testg)

test.listdf <- list()
for(i in 1:length(test.annomatching)){
  df <- data.frame(test.annomatching[i])
  test.listdf[[i]] <- df
}

names(test.listdf) <- names(testg)

for(i in 1:length(test.listdf)){
  test.listdf[[i]][2] <- NULL
  colnames(test.listdf[[i]]) = c("enzyme", "start", "end", "width")
  test.listdf[[i]][,1] <- table[test.listdf[[i]][,1],1]
}

cuts_per_gene <- sapply(test.listdf, function(x) table(x$enzyme))

cuts_all_genes <- rowSums(cuts_per_gene)

cuts_all_genes

