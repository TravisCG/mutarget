drawboxplot <- function(m, tmpp, n, mut, gname){
	png(paste(tmpp, gname, "png", sep = "."))
	if(m == "genotype"){
		boxplot()
	}
	else{
		boxplot()
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

normexp <- read.table(expfile, header = T, check.names = F)
mut     <- read.table(mutfile, header = T, check.names = F)
result  <- read.table(resultfile, header = T, check.names = F)

for(index in 1:plotnum){
	genename <- result[index,1]
	drawboxplot(mode, tmpprefix, normexp, mut, genename)
}
