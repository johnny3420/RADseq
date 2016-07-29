##Loading libraries

library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

##Reading in Brassica FASTA file
##Files from ftp://brassicadb.org/Brassica_napus/

DNA <- readDNAStringSet("ChBrGSP/Brassica_napus_v4.1.chromosomes.fa.gz")
Gene <-readDNAStringSet("CHBrGSP/Brassica_napus.annotation_v5.gff3.cds.fa.gz")

###Preparing Pdict for Biostrings

sequence <- "ACCGGT,ACCGGT,RAATTY,GGATCC,GCTAGC,GGTCTC,TGTACA,CGGCCG,GAATTC,GAATTC,GATATC,AAGCTT,GGTACC,GGTACC,CAATTG,ACGCGT,CCATGG,GCTAGC,GCTAGC,TCGCGA,ATGCAT,CTGCAG,CGATCG,CAGCTG,GAGCTC,GTCGAC,GTCGAC,AGTACT,ACTAGT,GCATGC,AATATT,CCWWGG,TCTAGA,CTCGAG,GCGGCCGC,GCGGCCGC,TTAATTAA,CCTGCAGG,CACNNNGTG,GGTNACC"
enzyme <- "AgeI-HF®,AgeI-HF® RE-Mix®,ApoI-HF,BamHI-HF®,BmtI-HF®,BsaI-HF®,BsrGI-HF®,EagI-HF®,EcoRI-HF®,EcoRI-HF® RE-Mix®,EcoRV-HF®,HindIII-HF®,KpnI-HF®,KpnI-HF® RE-Mix®,MfeI-HF®,MluI-HF®,NcoI-HF®,NheI-HF®,NheI-HF® RE-Mix®,NruI-HF®,NsiI-HF®,PstI-HF®,PvuI-HF®,PvuII-HF®,SacI-HF®,SalI-HF®,SalI-HF® RE-Mix®,ScaI-HF®,SpeI-HF®,SphI-HF®,SspI-HF®,StyI-HF®,XbaI RE-Mix®,XhoI RE-Mix®,NotI-HF®,NotI-HF® RE-Mix® ,PacI RE-Mix®,SbfI-HF®,DraIII-HF®,BstEII-HF®"
name1 <- strsplit(enzyme, ",")
name2 <- strsplit(sequence, ",")
table <- cbind.data.frame(enzyme = name1, sequence =name2)
colnames(table) <- c("enzyme", "sequence")
table <- table[c(-3,-32,-39,-40),] #Removing enzymes with ambiguity codes
pdict <- as.character(table[,2])
names(pdict) <- table[,1]
pdict <- DNAStringSet(pdict)
pdict <- PDict(pdict, tb.start = 1, tb.width = 6)

### 1st attempt at matching the pdict to the DNA Sequences

test <- vmatchPDict(pdict, DNA, fixed = "pattern") ### CANT BE RUN, FUNCTION IN PACKAGE INCOMPLETE NEED TO 
                                              ###RUN WITH AN INDIVIDUAL SUBJECT, SEQUENCE IS 41 PARTS,
                                              ###NEED RUN EACH TIME FOR EACH CHROMOSOME

### Counting occurance of each restriction site in each chromosome.
DNAcounting <- lapply(1:41, function(s) {
  countPDict(pdict, DNA[[s]])
})

names(DNAcounting)<- names(DNA)

### Putting results in dataframe

DNAcount.table <- data.frame(DNAcounting, row.names = table[,1])
DNAcount.table$totalcuts <- rowSums(DNAcount.table)
DNAcount.table$roughfraglen <- 1300000000/DNAcount.table$totalcuts

### Identifies the locations of each restriction site for each enzyme on each chromsome.
DNAmatching <- lapply(1:41, function(s) {
  matchPDict(pdict, DNA[[s]])
})
names(DNAmatching) <- names(DNA) ## Puts names on each chromosome.


### Identifies the location of each restriction site for each gene in the annotated genome.
Genematching <-lapply(1:length(Gene), function(s) {
  matchPDict(pdict,Gene[[s]])
})
names(Genematching) <-names(Gene)
###when indexing with famatching[[x]][[y]] 
###x is chromosome and y is RE

###Isolates one gene and lists the enzymes that cut it and the locations of the cuts as a test
test <- data.frame(Genematching$GSBRNA2T00121640001)
test$group <- table[test$group,1]
test[,1] <- table[test[,1],1]

###Creates a list of data frames of each gene and the enzymes
###which cut it and where. If an enzyme cuts more than once, each cut is listed.
###
genelist <- list()
for(i in 1:length(Genematching)){
  df <- data.frame(Genematching[i])
  genelist[[i]] <- df
}

names(genelist) <- names(Gene)

##Cleaning up the data frames
for(i in 1:length(genelist)){
  genelist[[i]][2] <- NULL
  colnames(genelist[[i]]) = c("enzyme", "start", "end", "width")
  genelist[[i]][,1] <- table[genelist[[i]][,1],1]
}

## Applying the confirmed correct enzyme names to data frames
for(i in 1:length(genelist)){
  genelist[[i]][,1] <- table[genelist[[i]][,1],1]
}
###
###

cuts_per_gene <- sapply(genelist, function(x) table(x$enzyme))

cuts_all_genes <- rowSums(cuts_per_gene)

cuts_all_genes

