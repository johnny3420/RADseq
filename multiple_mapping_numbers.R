### Looking at alternative hit numbers

###Importing data tables, X1:5 are ED of each alternative hit

fifty <- read.table("50_data")
hundred <- read.table("100_data")
onefifty <- read.table("150_data")

###Overall percentage of multiple mapping reads

paste("Of all 50bp reads", round(nrow(fifty)/252892*100, 3), "percent map to more than one site")
paste("Of all 100bp reads", round(nrow(hundred)/252892*100, 3), "percent map to more than one site")
paste("Of all 150bp reads", round(nrow(onefifty)/252893*100, 3), "percent map to more than one site")

###Percentage of multiple mapping reads which do not have exact alternative hits
exact <- sum(fifty$exact_matches == 0)
paste("of all 50bp multiple mapping reads", round(exact/nrow(fifty)*100, 3), "percent do not have alternative sites which are exact matches")
exact <- sum(hundred$exact_matches == 0)
paste("of all 100bp multiple mapping reads", round(exact/nrow(hundred)*100, 3), "percent do not have alternative sites which are exact matches")
exact <- sum(onefifty$exact_matches == 0)
paste("of all 150bp multiple mapping reads", round(exact/nrow(onefifty)*100, 3), "percent do not have alternative sites which are exact matches")

###Percentages of multiple mapping reads with exact matches
results50 <- character(0)
for(i in 1:5){
  exact <- sum(fifty$exact_matches == i)
  results50 <-c(results50, paste("Of all 50 bp reads", round(exact/252892*100,3), "percent map to", i, "alternative sites exactly"))
  results50 <-c(results50, paste("Of all 50 bp multiple mapping reads", round(exact/nrow(fifty)*100, 3), "percent map to", i, "alternative sites exactly"))
}
results100 <- character(0)
for(i in 1:5){
  exact <- sum(hundred$exact_matches == i)
  results100 <-c(results100, paste("Of all 100 bp reads", round(exact/252892*100,3), "percent map to", i, "alternative sites exactly"))
  results100 <-c(results100, paste("Of all 100 bp multiple mapping reads", round(exact/nrow(hundred)*100, 3), "percent map to", i, "alternative sites exactly"))
}
results150 <- character(0)
for(i in 1:5){
  exact <- sum(onefifty$exact_matches == i)
  results150 <-c(results150, paste("Of all 150 bp reads", round(exact/252893*100,3), "percent map to", i, "alternative sites exactly"))
  results150 <-c(results150, paste("Of all 150 bp multiple mapping reads", round(exact/nrow(onefifty)*100, 3), "percent map to", i, "alternative sites exactly"))
}
results50
results100
results150

###Reads containing alternative hits with Edit distance of 1 to 5+
edresults <- character(0)
for(i in 1:5){
  if(i == 5){
    edresults <- c(edresults, paste("Of all 50 bp multiple mapping reads", round(sum(fifty[,19+i] >= 1)/252892*100,3), "percent have alternative hit sites with edit distance greater than or equal to", i))
    edresults <- c(edresults, paste("Of all 100 bp multiple mapping reads", round(sum(hundred[,19+i] >= 1)/252892*100,3), "percent have alternative hit sites with edit distance greater than or equal to", i))
    edresults <- c(edresults, paste("Of all 150 bp multiple mapping reads", round(sum(onefifty[,19+i] >= 1)/252892*100,3), "percent have alternative hit sites with edit distance greater than or equal to", i))
  }else if(i <= 4){
    edresults <- c(edresults, paste("Of all 50 bp multiple mapping reads", round(sum(fifty[,19+i] >= 1)/252892*100,3), "percent have alternative hit sites with edit distance of", i))
    edresults <- c(edresults, paste("Of all 100 bp multiple mapping reads", round(sum(hundred[,19+i] >= 1)/252892*100,3), "percent have alternative hit sites with edit distance of", i))
    edresults <- c(edresults, paste("Of all 150 bp multiple mapping reads", round(sum(onefifty[,19+i] >= 1)/252892*100,3), "percent have alternative hit sites with edit distance of", i))
}}
edresults
###Reads containg alternative 

