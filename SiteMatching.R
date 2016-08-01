##Loading libraries

library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

##Reading in Brassica FASTA file
##Files from ftp://brassicadb.org/Brassica_napus/

DNA <- readDNAStringSet("ChBrGSP/Brassica_napus_v4.1.chromosomes.fa.gz")
Gene <-readDNAStringSet("CHBrGSP/Brassica_napus.annotation_v5.gff3.cds.fa.gz")

###Preparing Pdict for Biostrings

sequence <- "ACCGGT,ACCGGT,GGATCC,GCTAGC,GGTCTC,TGTACA,CGGCCG,GAATTC,GAATTC,GATATC,AAGCTT,GGTACC,GGTACC,CAATTG,ACGCGT,CCATGG,GCTAGC,GCTAGC,TCGCGA,ATGCAT,CTGCAG,CGATCG,CAGCTG,GAGCTC,GTCGAC,GTCGAC,AGTACT,ACTAGT,GCATGC,AATATT,TCTAGA,CTCGAG,GCGGCCGC,GCGGCCGC,TTAATTAA,CCTGCAGG"
enzyme <- "AgeI-HF®,AgeI-HF® RE-Mix®,BamHI-HF®,BmtI-HF®,BsaI-HF®,BsrGI-HF®,EagI-HF®,EcoRI-HF®,EcoRI-HF® RE-Mix®,EcoRV-HF®,HindIII-HF®,KpnI-HF®,KpnI-HF® RE-Mix®,MfeI-HF®,MluI-HF®,NcoI-HF®,NheI-HF®,NheI-HF® RE-Mix®,NruI-HF®,NsiI-HF®,PstI-HF®,PvuI-HF®,PvuII-HF®,SacI-HF®,SalI-HF®,SalI-HF® RE-Mix®,ScaI-HF®,SpeI-HF®,SphI-HF®,SspI-HF®,XbaI RE-Mix®,XhoI RE-Mix®,NotI-HF®,NotI-HF® RE-Mix® ,PacI RE-Mix®,SbfI-HF®"
name1 <- strsplit(enzyme, ",")
name2 <- strsplit(sequence, ",")
table <- cbind.data.frame(enzyme = name1, sequence =name2)
colnames(table) <- c("enzyme", "sequence")
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

###Creates list of data frames of each chromosome and the enzymes which cut there.

dnalist<- list()
for(i in 1:length(DNAmatching)){
  df <- data.frame(DNAmatching[i])
  dnalist[[i]] <- df
}
names(dnalist) <- names(DNA)

###Cleaning up dataframes

for(i in 1:length(dnalist)){
  dnalist[[i]][2] <- NULL
  colnames(dnalist[[i]]) = c("enzyme", "start", "end", "width")
  dnalist[[i]][,1] <- table[dnalist[[i]][,1],1]
}

### Counting number of cuts of each restriction enzyme on each gene in annotated genome

Genecounting <- lapply(1:101040, function(s) {
  countPDict(pdict, Gene[[s]])
})

names(Genecounting)<- names(Gene)

### Putting results in dataframe

Genecount.table <- data.frame(Genecounting, row.names = table[,1])
Genecount.table$totalcuts <- rowSums(Genecount.table)
Genecount.table$roughfraglen <- 1300000000/Genecount.table$totalcuts

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

###
### Table of cuts for each RE on genes

cuts_per_gene <- sapply(genelist, function(x) table(x$enzyme))

cuts_all_genes <- rowSums(cuts_per_gene)

cuts_all_genes

### Table of cuts for each RE on chromosomes
### NUMBER IS INFLATED I THINK. DNA SEQUENCES CONTAIN
### 41 NAMED CHROMOSOMES AND REPEATS ARE MOST LIKELY PRESENT
### NEED TO SUBSET "dnalist" BASED ON WHICH COLLECTION OF 
### CHROMOSOMES IS CHOSEN

cuts_per_chromosome <- sapply(dnalist, function(x) table(x$enzyme))

cuts_all_chromosomes <- rowSums(cuts_per_chromosome)

cuts_all_chromosomes

### Data frame comparing cut comparisons

cut_comparison <- cbind.data.frame(whole = cuts_all_chromosomes, genes = cuts_all_genes)

cut_comparison

cut_comparison$whole2gene <- cut_comparison$whole/cut_comparison$genes
cut_comparison$gene2whole <- cut_comparison$gene/cut_comparison$whole
cut_comparison$rough_frag_length <- 1300000000/cut_comparison$whole

### Prices of RE

units <- c(300, 0, 10000, 300, 1000, 1000, 500, 10000, 0, 4000, 10000, 4000, 0, 500, 1000, 1000, 1000, 0, 500, 0, 1000, 1000, 0, 10000, 500, 5000, 2000, 2000, 0, 500, 1000, 500, 500, 1000, 0, 0)
bulk_units <- c(1500, 0, 50000, 1500, 5000, 5000, 2500, 50000, 0, 20000, 50000, 2000, 0, 2500, 5000, 5000, 5000, 0, 2500, 0, 5000, 5000, 0, 50000, 2500, 25000, 10000, 10000, 0, 2500, 5000, 2500, 2500, 5000, 0, 0)
price <- c(69, 0, 57, 65, 65, 65, 62, 57, 0, 57, 57, 62, 0, 70, 62, 60, 65, 0, 70, 0, 57, 62, 0, 62, 69, 57, 57, 60, 0, 70, 62, 65, 65, 67, 0, 0)
bulk_price <- c(278, 0, 229, 262, 262, 262, 249, 229, 0, 229, 229, 249, 0, 282, 249, 241, 262, 0, 282, 0, 229, 249, 0, 249, 278, 229, 229, 241, 0, 282, 249, 262, 262, 270, 0, 0)

cut_comparison$units <- units
cut_comparison$price <- price
cut_comparison$units_per_dollar <- units/price
cut_comparison$bulk_units <- bulk_units
cut_comparison$bulk_price <- bulk_price
cut_comparison$bulk_units_per_dollar <- bulk_units/bulk_price

cut_comparison

###

