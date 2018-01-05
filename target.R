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

# Command line arguments
argv         <- commandArgs(trailing = T)
tmpprefix    <- argv[1]
cancerid     <- argv[2]
effect       <- argv[3]
testgene     <- argv[4]
qvalcutoff   <- as.numeric(argv[5])
foldchcutoff <- as.numeric(argv[6])
mutprev      <- as.numeric(argv[7])
filtergene   <- argv[8]
filterout    <- argv[9]
proc.time()
print("MESSAGE: Start")
# MySQL connection
con  <- dbConnect(MySQL(), user="XXXX", password="XXXX", dbname="mutarget", host="localhost")

# Expression matrix
count <- as.matrix(read.table(paste(cancerid, "tsv", sep = "."), check.names = F, sep = "\t"))
proc.time()
print("MESSAGE: Exp matrix")

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
	rs      <- dbSendQuery(con, query)
	raw     <- fetch(rs, n=-1)
	index   <- grep(paste(raw$samples, collapse="|"), rownames(coldata))
	coldata <- coldata[index,,drop=F]
	count   <- count[,index]
}
proc.time()
print("MESSAGE: Sample filtering")

# Normalise raw counts
edge <- DGEList(counts = count, group = coldata$gene)
edge <- calcNormFactors(edge)
exp  <- cpm(edge)
proc.time()
print("MESSAGE: Normalisation")

# Select the gene of interest
winp <- data.frame(exp = exp[testgene,], mutant = 0)

# Get all the genes
shortnames <- sub("-...-...-....-..$","",rownames(coldata))
maxcount   <- length(unique(shortnames)) # maximum number of mutation should be less than all the samples (prevent one group syndrome)
mincount   <- trunc(maxcount * mutprev / 100) # minimum number of mutation calculated from mutation prevalence

query <- paste("select name,genename from mutation inner join (genetable,individual) on (genetable_geneid = geneid and patientid = individual_patientid) where individual_cancerid = ",cancerid," and muteffect_effectid = ",effect,";", sep = "")

rs <- dbSendQuery(con, query)
mutmatrix <- fetch(rs, n=-1)
mutmatrix <- as.data.frame.matrix(table(mutmatrix$genename, mutmatrix$name))
mutmatrix[mutmatrix > 0] <- 1 # We need a binary matrix
mutmatrix <- mutmatrix[rowSums(mutmatrix) > mincount & rowSums(mutmatrix) < maxcount,]
result_table <- data.frame(foldchange = rep(0, nrow(mutmatrix)), pvalue = rep(1,nrow(mutmatrix)), adj.pval = rep(1,nrow(mutmatrix)))
rownames(result_table) <- rownames(mutmatrix)

# Iterate through genes and calculate Wilcox-test
proc.time()
print("MESSAGE: Start iteration")
r <- apply(mutmatrix, 1, runWilcox, winp)
r <- t(r)
result_table$pvalue <- r[,1]
result_table$foldchange <- r[,2]
proc.time()
print("MESSAGE: End iteration")

# Multiple test correction
result_table$adj.pval <- p.adjust(result_table$pvalue, "BH")

#Filtering and producing table
result_table <- result_table[result_table$adj.pval < qvalcutoff & result_table$foldchange > foldchcutoff,]
write.table(format(result_table, digits = 2), paste(tmpprefix, "tsv",sep="."), quote=F, sep ="\t")

# Decision tree
mutmatrix2 <- data.frame(lapply(mutmatrix, factor, levels = c(0,1), labels=c("WT","Mut")), check.names = F)
rownames(mutmatrix2) <- rownames(mutmatrix)
ctree <- cbind(winp[,1,drop=F], t(mutmatrix2[,shortnames]))
resulttree <- ctree(exp ~ ., data = ctree, , control = ctree_control(maxdepth=3, minprob=mutprev/100,testtype=c("Univariate")))
proc.time()
print("MESSAGE: Decision tree created")

png(paste(tmpprefix,"dectree.png",sep="."))
plot(resulttree)
dev.off()

# Boxplots
proc.time()
print("MESSAGE: Start creating boxplots")
for(gene in rownames(result_table)){
	f <- as.formula(paste("exp ~ ", gene, sep=""))
	png(paste(tmpprefix,gene,"png",sep="."))
	boxplot(f, data = ctree, main = paste(testgene, "expression", sep = " "), xlab = paste(gene, "mutation status", sep = " "))
	dev.off()
}

proc.time()
print("MESSAGE: End script")
