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

### Putting list of MIndex objects into list of data frames
chromolist<- list()
for(i in 1:length(sitematches)){
  df <- data.frame(sitematches[i])
  chromolist[[i]] <- df
}
names(chromolist) <- names(DNA)

### Cleaning up data frames

for(i in 1:length(chromolist)){
  chromolist[[i]][c(1,2,5)] <- NULL
  chromolist[[i]][,3] <- NA
  colnames(chromolist[[i]]) = c("start", "end", "test")
}

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

### Overall though process
### Create a very long string composed of each read
### Take that string and turn into DNAStringSet
### Take that string set and turn into FASTQ file

test <- character(0) ### create string
test <- paste(DNA[[i]][start:end]) ### format for getting needed sequence string
test <- c(test, "") ### Appending string
test <- DNAStringSet(test) ### converting string into DNAStringSet
writeXStringSet(test,"", format = "fastq") ### Creating fastq file

### Successful test on chromosome chrA01 for all acceptable 50 bp reads being written to FASTQ file
### Very slow process

fifty <- character(0)
for(i in 1){
  k <- 1
  while(k <= nrow(chromolist$chrA01)){
    if(((chromolist[[1]][k+1,3] == TRUE)) & (k == 1)){
      start <- chromolist[[1]][k,2]-49
      end <- chromolist[[1]][k,2]
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      start <- chromolist[[1]][k,1]
      end <- chromolist[[1]][k,1]+49
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty, seq)
      k <- k +1
    } else if(((chromolist[[1]][k+1,3] == FALSE)) & (k == 1)){
      start <- chromolist[[1]][k,2]-49
      end <- chromolist[[1]][k,2]
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)             
      k <- k + 1
    } else if((k == nrow(chromolist[[1]])) & (chromolist[[1]][k,3] == TRUE)){
      start <- chromolist[[1]][k,2]-49
      end <- chromolist[[1]][k,2]
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      start <- chromolist[[1]][k,1]
      end <- chromolist[[1]][k,1]+49
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((k == nrow(chromolist[[1]])) & (chromolist[[1]][k,3] == FALSE)){
      start <- chromolist[[1]][k,1]
      end <- chromolist[[1]][k,1]+49
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[1]][k,3] == TRUE) & (chromolist[[1]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[1]])){
      start <- chromolist[[1]][k,2]-49
      end <- chromolist[[1]][k,2]
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[1]][k,3] == FALSE) & (chromolist[[1]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[1]])){
      start <- chromolist[[1]][k,1]
      end <- chromolist[[1]][k,1]+49
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[1]][k,3] == TRUE) & (chromolist[[1]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[1]])){
      start <- chromolist[[1]][k,2]-49
      end <- chromolist[[1]][k,2]
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      start <- chromolist[[1]][k,1]
      end <- chromolist[[1]][k,1]+49
      seq <- paste(DNA[[1]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[1]][k,3] == FALSE) & (chromolist[[1]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[1]])){
      k <- k + 1
    }
  }
}

test <- DNAStringSet(fifty)
names(test) <- 1:length(test)
test
writeXStringSet(test, "chromo1.fq", format ="fastq")
check <- readDNAStringSet("chromo1.fq", format = "fastq")
check

### 50bp FASTQ creation

fifty <- character(0)
for(i in 1:length(chromolist)){
  k <- 1
  while(k <= nrow(chromolist[[i]])){
    if(((chromolist[[i]][k+1,3] == TRUE)) & (k == 1)){
      start <- chromolist[[i]][k,2]-49
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+49
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty, seq)
      k <- k +1
    } else if(((chromolist[[i]][k+1,3] == FALSE)) & (k == 1)){
      start <- chromolist[[i]][k,2]-49
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)             
      k <- k + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == TRUE)){
      start <- chromolist[[i]][k,2]-49
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+49
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == FALSE)){
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+49
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == TRUE) & (chromolist[[i]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,2]-49
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == FALSE) & (chromolist[[i]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+49
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == TRUE) & (chromolist[[i]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,2]-49
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+49
      seq <- paste(DNA[[i]][start:end])
      fifty <- c(fifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == FALSE) & (chromolist[[i]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[i]])){
      k <- k + 1
    }
  }
}

fastq50 <- DNAStringSet(fifty)
names(fastq50) <- 1:length(fastq50)
writeXStringSet(fastq50, "50bp.fq")

### 100bp FASTQ creation

hundred <- character(0)
for(i in 1:length(chromolist)){
  k <- 1
  while(k <= nrow(chromolist[[i]])){
    if(((chromolist[[i]][k+1,3] == TRUE)) & (k == 1)){
      start <- chromolist[[i]][k,2]-99
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+99
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred, seq)
      k <- k +1
    } else if(((chromolist[[i]][k+1,3] == FALSE)) & (k == 1)){
      start <- chromolist[[i]][k,2]-99
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)             
      k <- k + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == TRUE)){
      start <- chromolist[[i]][k,2]-99
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+99
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      k <- k + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == FALSE)){
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+99
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == TRUE) & (chromolist[[i]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,2]-99
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == FALSE) & (chromolist[[i]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+99
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == TRUE) & (chromolist[[i]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,2]-99
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+99
      seq <- paste(DNA[[i]][start:end])
      hundred <- c(hundred,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == FALSE) & (chromolist[[i]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[i]])){
      k <- k + 1
    }
  }
}

fastq100 <- DNAStringSet(hundred)
names(fastq100) <- 1:length(fastq100)
writeXStringSet(fastq100, "100bp.fq")


### 150bp FASTQ creation

onefifty <- character(0)
for(i in 1:length(chromolist)){
  k <- 1
  while(k <= nrow(chromolist[[i]])){
    if(((chromolist[[i]][k+1,3] == TRUE)) & (k == 1)){
      start <- chromolist[[i]][k,2]-149
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+149
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty, seq)
      k <- k +1
    } else if(((chromolist[[i]][k+1,3] == FALSE)) & (k == 1)){
      start <- chromolist[[i]][k,2]-149
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)             
      k <- k + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == TRUE)){
      start <- chromolist[[i]][k,2]-149
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+149
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      k <- k + 1
    } else if((k == nrow(chromolist[[i]])) & (chromolist[[i]][k,3] == FALSE)){
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+149
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == TRUE) & (chromolist[[i]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,2]-149
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == FALSE) & (chromolist[[i]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+149
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == TRUE) & (chromolist[[i]][k+1,3] == TRUE)
              & k != 1 & k != nrow(chromolist[[i]])){
      start <- chromolist[[i]][k,2]-149
      end <- chromolist[[i]][k,2]
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      start <- chromolist[[i]][k,1]
      end <- chromolist[[i]][k,1]+149
      seq <- paste(DNA[[i]][start:end])
      onefifty <- c(onefifty,seq)
      k <- k + 1
    } else if((chromolist[[i]][k,3] == FALSE) & (chromolist[[i]][k+1,3] == FALSE)
              & k != 1 & k != nrow(chromolist[[i]])){
      k <- k + 1
    }
  }
}

fastq150 <- DNAStringSet(onefifty)
names(fastq150) <- 1:length(fastq150)
writeXStringSet(fastq150, "150bp.fq")



###