library(Rsamtools)
library(stringr)
### Creating bam files and indexes from sam files
### Only needs to be ran if files not present
### Won't run if files already exist

asBam("fifty.sam", "Fastq/fifty", overwrite=FALSE, indexDestination=TRUE)
asBam("hundred.sam", "Fastq/hundred", overwrite=FALSE, indexDestination=TRUE)
asBam("onefifty.sam", "Fastq/onefifty", overwrite=FALSE, indexDestination=TRUE)

### Necessary Param object to include desired tags

p <- ScanBamParam(tag = "XA", what = c("qname","flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar","mrnm", "mpos","isize","seq"))

### Importing bam files

fifty <- scanBam("Fastq/fifty.bam", index = "Fastq/fifty.bam.bai", param = p)
hundred <- scanBam("Fastq/hundred.bam", index = "Fastq/hundred.bam.bai", param = p)
onefifty <- scanBam("Fastq/onefifty.bam", index = "Fastq/onefifty.bam.bai", param = p)

### Creating data frames from bam files


df_50 <- data.frame(fifty)
df_100 <- data.frame(hundred)
df_150 <- data.frame(onefifty)

### Counting unique mapping read (map to one location)

unique50 <- sum(is.na(df_50$XA))
unique100 <- sum(is.na(df_100$XA))
unique150 <- sum(is.na(df_150$XA))

### Percentage unique mapping overall

unique50/nrow(df_50)*100
unique100/nrow(df_100)*100
unique150/nrow(df_150)*100

### Isolating reads that mapped to more than one site

multi_df_50 <- subset(df_50, is.na(df_50$XA)!= TRUE)
multi_df_100 <-subset(df_100, is.na(df_100$XA)!= TRUE)
multi_df_150 <-subset(df_150, is.na(df_150$XA)!= TRUE)

### Splitting XA column for all 3 sets and place into new data frames

XA50 <- as.character(unlist(multi_df_50$XA))
XA50 <- sub(".*,","",XA50)
XA50 <- data.frame(as.numeric(sub(";","",XA50)))
XA100 <- as.character(unlist(multi_df_100$XA))
XA100 <- sub(".*,","",XA100)
XA100 <- data.frame(as.numeric(sub(";","",XA100)))
XA150 <- as.character(unlist(multi_df_150$XA))
XA150 <- sub(".*,","",XA150)
XA150 <- data.frame(as.numeric(sub(";","",XA150)))

### Looking at 50bp Alternatives

sapply(0:max(XA50), function(k){
  if(k == 0){
    paste("Of multiple mapping reads,", round(sum(XA50 == k)/nrow(XA50)*100, digit = 3), "percent are an exact match")
  }else
    paste("Of multiple mapping reads,", round(sum(XA50 == k)/nrow(XA50)*100, digit = 3), "percent have an Edit Distance of", k)
})

### Looking at 100bp Alternatives

sapply(0:max(XA100), function(k){
  if(k == 0){
    paste("Of multiple mapping reads,", round(sum(XA100 == k)/nrow(XA50)*100, digit = 3), "percent are an exact match")
  }else
    paste("Of multiple mapping reads,", round(sum(XA100 == k)/nrow(XA50)*100, digit = 3), "percent have an Edit Distance of", k)
})

### Looking at 150bp Alternatives

sapply(0:max(XA150), function(k){
  if(k == 0){
    paste("Of multiple mapping reads,", round(sum(XA150 == k)/nrow(XA50)*100, digit = 3), "percent are an exact match")
  }else
    paste("Of multiple mapping reads,", round(sum(XA50 == k)/nrow(XA150)*100, digit = 3), "percent have an Edit Distance of", k)
})

