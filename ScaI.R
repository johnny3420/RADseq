### ScaI fragments produced by cutting the Brassica Napus genome
library(Biostrings)

### Loading in genome

DNA <- readDNAStringSet("ChBrGSP/Brassica_napus_v4.1.chromosomes.fa.gz")

### Preparing PDict for matching

sequence <- "GAGCTC"
enzyme <- "SacI"
ref <- cbind.data.frame(enzyme = enzyme, sequence = sequence)
pdict <- as.character(ref[,2])
names(pdict) <- ref[,1]
pdict <- DNAStringSet(pdict)
pdict <- PDict(pdict)

### Matching restriction sites to genome

sitematches <- lapply(1:length(DNA), function(s) {
  matchPDict(pdict, DNA[[s]])
})
names(sitematches) <- names(DNA) ## Puts names on each chromosome.

### Putting list in list of data frames
chromolist<- list()
for(i in 1:length(sitematches)){
  df <- data.frame(sitematches[i])
  chromolist[[i]] <- df
}
names(chromolist) <- names(DNA)

### Cleaning up data frames
for(i in 1:length(chromolist)){
  chromolist[[i]][c(1,2,5)] <- NULL
  colnames(chromolist[[i]]) = c("start", "end")
}

### Creating 50bp from site fragment list ignoring <200 bp site gaps
### few things need to be tweaked
### CURRENTLY DOES NOT WORK

fiftybplist <- list()
l <- 1
for(i in 1:length(chromolist)){
  for(k in 1:nrow(chromolist[[i]]-1)){
    start1 <- chromolist[[i]][1,1]-49
    end1 <- chromolist[[i]][1,1]
    start2 <- chromolist[[i]][1,2]
    end2 <- chromolist[[i]][1,2]+49
    fiftybplist[[l]] <- DNA[[i]][start1:end1]
    l <- l + 1
    fiftybplist[[l]] <- DNA[[i]][start2:end2]
    l <- l + 1
      if((chromolist[[i]][k+1,1] - chromolist[[i]][k,2]) >= 200){
        start3 <- chromolist[[i]][k+1,1]-49
        end3 <- chromolist[[i]][k+1,1]
        start4 <- chromolist[[i]][k+1,2]
        end4 <- chromolist[[i]][k+1,2]+49
        fiftybplist[[l]] <- DNA[[i]][start3:end3]
        l <- l + 1
        fiftybplist[[l]] <- DNA[[i]][start4:end4]
        l <- l + 1
  }
  }}


### Test to see if any sites are within 200 bp of eachother on each chromosome
### TRUE == minimum gap achieved, FALSE == minimum gap failed

for(i in 1:length(chromolist)){
  for(k in 1:nrow(chromolist[[i]])){
      if(k == 1){
        chromolist[[i]][k,3] <- TRUE
      }else if(k == nrow(chromolist[[i]])){
        chromolist[[i]][k,3] <- TRUE
      }else if((chromolist[[i]][k,1] - chromolist[[i]][k-1,2]) >= 200){
        chromolist[[i]][k,3] <- TRUE
      } else
        chromolist[[i]][k,3] <- FALSE
      }
    
  }
  
### Sets up dataframe for new additions
  
for(i in 1:length(chromolist)){
  chromolist[[i]][,4] <- NA
  chromolist[[i]][,5] <- NA
  colnames(chromolist[[i]]) = c("start", "end", "test", "upstream", "downstream")
}

### Adds upstream and downstream DNA sequences when applicable for 50 bp lengths
fiftybplist <- list()
l <-1
for(i in 1:length(chromolist)){
  for(k in 1:nrow(chromolist[[i]])){
    if(((chromolist[[i]][k+1,3] == TRUE)) & (k == 1)){
      start <- chromolist[[i]][k,1]-49
      end <- chromolist[[i]][k,1]
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
      start <- chromolist[[i]][k,2]
      end <- chromolist[[i]][1,2]+49
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
    } else if(((chromolist[[i]][k+1,3] == FALSE)) & (k == 1)){
      start <- chromolist[[i]][k,1]-49
      end <- chromolist[[i]][k,1]
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == TRUE)){
      start <- chromolist[[i]][k,1]-49
      end <- chromolist[[i]][k,1]
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
      start <- chromolist[[i]][k,2]
      end <- chromolist[[i]][1,2]+49
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == FALSE)){
      start <- chromolist[[i]][k,2]
      end <- chromolist[[i]][1,2]+49
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
    } else if((chromolist[[i]][k,3] == TRUE) & chromolist[[i]][k+1,3] == FALSE){
      start <- chromolist[[i]][k,1]-49
      end <- chromolist[[i]][k,1]
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
    } else if((chromolist[[i]][k,3] == FALSE) & chromolist[[i]][k+1,3] == TRUE){
      start <- chromolist[[i]][k,2]
      end <- chromolist[[i]][1,2]+49
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
    } else if((chromolist[[i]][k,3] == TRUE) & chromolist[[i]][k+1,3] == TRUE){
      start <- chromolist[[i]][k,1]-49
      end <- chromolist[[i]][k,1]
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
      start <- chromolist[[i]][k,2]
      end <- chromolist[[i]][1,2]+49
      fiftybplist[[l]] <- DNA[[i]][start:end]
      l <- l + 1
    }
  }}
