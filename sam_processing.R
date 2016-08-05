library(Rsamtools)
library(stringr)
library(plyr)
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

### Isolating ED of each multiple mapping read

XA50 <- as.character(unlist(multi_df_50$XA))
XA50 <- (str_split(XA50, ";"))
XA50 <- sapply(XA50, function(x) as.numeric(gsub(".*,","",x)))
for(i in 1:length(XA50)){
  XA50[[i]] <- XA50[[i]][1:(length(XA50[[i]])-1)]
}
XA50 <- ldply(XA50, rbind)

XA100 <- as.character(unlist(multi_df_100$XA))
XA100 <- (str_split(XA100, ";"))
XA100 <- sapply(XA100, function(x) as.numeric(gsub(".*,","",x)))
for(i in 1:length(XA100)){
  XA100[[i]] <- XA100[[i]][1:(length(XA100[[i]])-1)]
}
XA100 <- ldply(XA100, rbind)

XA150 <- as.character(unlist(multi_df_150$XA))
XA150 <- (str_split(XA150, ";"))
XA150 <- sapply(XA150, function(x) as.numeric(gsub(".*,","",x)))
for(i in 1:length(XA150)){
  XA150[[i]] <- XA150[[i]][1:(length(XA150[[i]])-1)]
}
XA150 <- ldply(XA150, rbind)

### Creating new data table with ED for each alternative site for each read

stretch50 <- cbind(multi_df_50, XA50)
stretch100 <- cbind(multi_df_100, XA100)
stretch150 <- cbind(multi_df_150,XA150)


### Percentage of multiple site reads with exact match alternative hits

exact <- nrow(stretch50[which(stretch50$`1` == 0 | stretch50$`2` == 0 
                               | stretch50$`3` == 0 | stretch50$`4` == 0 
                               | stretch50$`5` == 0), ])
paste("Of all 50bp reads", round(exact/nrow(df_50)*100, 3), "percent map to at least 2 sites exactly")
paste("Of all 50bp multiple mapping reads", round(exact/nrow(multi_df_50)*100, 3), "percent have at least one exact alternate hit")

exact <- nrow(stretch100[which(stretch100$`1` == 0 | stretch100$`2` == 0 
            | stretch100$`3` == 0 | stretch100$`4` == 0 
            | stretch100$`5` == 0), ])
paste("Of all 100bp reads", round(exact/nrow(df_100)*100, 3), "percent map to at least 2 sites exactly")
paste("Of all 100bp multiple mapping reads", round(exact/nrow(multi_df_100)*100, 3), "percent have at least one exact alternate hit")


exact <- nrow(stretch150[which(stretch150$`1` == 0 | stretch150$`2` == 0 
                               | stretch150$`3` == 0 | stretch150$`4` == 0 
                               | stretch150$`5` == 0), ])
paste("Of all 150bp reads", round(exact/nrow(df_150)*100, 3), "percent map to at least 2 sites exactly")
paste("Of all 150bp multiple mapping reads", round(exact/nrow(multi_df_150)*100, 3), "percent have at least one exact alternate hit")


####