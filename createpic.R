drawboxplot <- function(m, ingene, tmpp, n, mut, gname, outliers){
	png(paste(tmpp, gname, "png", sep = "."))
	if(m == "genotype"){
		muts <- mut[mut$genes == ingene,] == 1
		wts  <- mut[mut$genes == ingene,] == 0
		boxplot(as.numeric(n[n$genes == gname, muts]), as.numeric(n[n$genes == gname, wts]), main = paste(gname, "expression", sep=" "), names = c("Mutant", "WT"), outline = outliers)
	}
	else{
		muts <- mut[mut$genes == gname,] == 1
		wts  <- mut[mut$genes == gname,] == 0
		boxplot(as.numeric(n[n$genes == ingene, muts]), as.numeric(n[n$genes == ingene, wts]), main = paste(ingene, "expression", sep = " "), names = c("Mutant", "WT"), xlab = paste(gname, "mutation status", sep = " "), outline = outliers)
	}
	dev.off()
}

argv       <- commandArgs(trailing = T)
mode       <- argv[1]
inputgene  <- argv[2]
tmpprefix  <- argv[3]
resultfile <- argv[4]
expfile    <- argv[5]
mutfile    <- argv[6]
plotnum    <- as.numeric(argv[7])
outliers   <- as.numeric(argv[8]) == 1

normexp <- read.table(expfile, header = T, check.names = F)
mut     <- read.table(mutfile, header = T, check.names = F)
result  <- read.table(resultfile, header = T, check.names = F)

for(index in 1:plotnum){
	genename <- result[index,1]
	drawboxplot(mode, inputgene, tmpprefix, normexp, mut, genename, outliers)
	cat(paste("MESSAGE: ",index * 100 / plotnum,"%\n",sep=""))
}
