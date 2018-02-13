# Script for mutarget target analysis

suppressMessages(library(edgeR))
suppressMessages(library(RMySQL))
suppressMessages(library(partykit))

runWilcox <- function(x, winp){
        samples <- names(x)[x > 0]
        index <- grep(paste(samples, collapse="|"), rownames(winp))
	if(length(index) == 0){
		return(c(1.0, 0))
	}
        winp[index,]$mutant <- 1
        testres <- wilcox.test(exp~mutant, data = winp)
        expmu <- median(winp[winp$mutant == 1,1])
        expwt <- median(winp[winp$mutant == 0,1])
        foldchange <- max(expmu, expwt) / min(expmu, expwt)
        return(c(testres$p.value, foldchange))
}
source("common.R")

# Command line arguments
argv         <- commandArgs(trailing = T)
tmpprefix    <- argv[1]
cancerid     <- argv[2]
effect       <- argv[3]
testgene     <- argv[4]
qvalcutoff   <- as.numeric(argv[5])
foldchcutoff <- as.numeric(argv[6])
mutprev      <- as.numeric(argv[7])
dbsrc        <- as.numeric(argv[8])
dodt         <- argv[9]                  # Create decesion tree?
filtergene   <- argv[10]
filterout    <- argv[11]
proc.time()
cat("MESSAGE: Initialisation\n")

# Expression matrix
count <- getExpMatrix(con, cancerid, dbsrc)
proc.time()
cat("MESSAGE: Expression matrix\n")

# Column data
coldata <- data.frame(gene = rep("WT", ncol(count)))
rownames(coldata) <- colnames(count)

# Sample filtering
if(!is.na(filtergene)){
	if(filterout == 'include'){
		query <- paste("select distinct(name) as samples from individual inner join (mutation,genetable) on (individual_patientid = patientid and genetable_geneid = geneid) where cancer_cancerid = ",cancerid," and genename = '",filtergene,"';",sep="")
	} else {
		query <- paste("select distinct(name) as samples from individual where cancer_cancerid = ",cancerid," and name not in ( select distinct(name) from individual inner join (mutation,genetable) on (individual_patientid = patientid and genetable_geneid = geneid) where cancer_cancerid = ",cancerid," and genename = '",filtergene,"' );",sep="")
	}
	raw     <- fetchDB(con, query)
	index   <- grep(paste(raw$samples, collapse="|"), rownames(coldata))
	coldata <- coldata[index,,drop=F]
	count   <- count[,index]
}
proc.time()
cat("MESSAGE: Sample filtering\n")

# Get all the genes
if (dbsrc == 2) {
	# Microarray data
	shortnames <- sub("BRCA-", "", rownames(coldata))
	exp <- count # Microarray data already normalised
} else {
	shortnames <- sub("(^[^-]+-[^-]+-[^-]+).+", "\\1", rownames(coldata))
	edge <- DGEList(counts = count, group = coldata$gene)
	edge <- calcNormFactors(edge)
	exp  <- cpm(edge)
}
proc.time()
cat("MESSAGE: Normalisation\n")

# Select gene of interest
winp <- data.frame(exp = exp[testgene,], mutant = 0)
maxcount   <- length(unique(shortnames)) # maximum number of mutation should be less than all the samples (prevent one group syndrome)
mincount   <- trunc(maxcount * mutprev / 100) # minimum number of mutation calculated from mutation prevalence

mutmatrix <- read.table(paste("mut",cancerid,dbsrc,"tsv",sep="."), check.names = F, sep="\t")
mutmatrix <- mutmatrix[rowSums(mutmatrix) > mincount & rowSums(mutmatrix) < maxcount,]
result_table <- data.frame(foldchange = rep(0, nrow(mutmatrix)), pvalue = rep(1,nrow(mutmatrix)), adj.pval = rep(1,nrow(mutmatrix)))
rownames(result_table) <- rownames(mutmatrix)

# Iterate through genes and calculate Wilcox-test
proc.time()
cat("MESSAGE: Start iteration\n")
r <- apply(mutmatrix, 1, runWilcox, winp)
r <- t(r)
result_table$pvalue <- r[,1]
result_table$foldchange <- r[,2]
proc.time()
cat("MESSAGE: End iteration\n")

# Multiple test correction
result_table$adj.pval <- p.adjust(result_table$pvalue, "BH")

#Filtering and producing table
result_table <- result_table[result_table$adj.pval < qvalcutoff & result_table$foldchange > foldchcutoff,]
write.table(format(result_table, digits = 2), paste(tmpprefix, "tsv",sep="."), quote=F, sep ="\t")
proc.time()
cat("MESSAGE: Create decision tree\n")

# Decision tree
if (dodt == "1") {
	mutmatrix2 <- data.frame(lapply(mutmatrix, factor, levels = c(0,1), labels=c("WT","Mut")), check.names = F)
	rownames(mutmatrix2) <- rownames(mutmatrix)
	comm <- intersect(shortnames, colnames(mutmatrix)) # Get common elements to avoid missing rows between mutation and expression matrix
	reg <- regexpr(paste(comm, collapse = "|"), rownames(winp))
	comm2 <- regmatches(rownames(winp), reg) # Get short names again from expression matrix, because some times a person has more than one expression data
	ctree <- cbind(winp[grep(paste(comm, collapse="|"),rownames(winp)),1,drop=F], t(mutmatrix2[,comm2]))
	resulttree <- ctree(exp ~ ., data = ctree, , control = ctree_control(maxdepth=3, minprob=0.05,testtype=c("Univariate")))
}
proc.time()
cat("MESSAGE: Decision tree created\n")

png(paste(tmpprefix,"dectree.png",sep="."))
plot(resulttree)
dev.off()

# Boxplots
proc.time()
cat("MESSAGE: Start creating boxplots\n")
for(gene in rownames(result_table)){
	f <- as.formula(paste("exp ~ ", gene, sep=""))
	png(paste(tmpprefix,gene,"png",sep="."))
	boxplot(f, data = ctree, main = paste(testgene, "expression", sep = " "), xlab = paste(gene, "mutation status", sep = " "), outline = F)
	dev.off()
}

proc.time()
cat("MESSAGE: End script\n")
