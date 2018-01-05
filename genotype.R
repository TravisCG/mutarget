suppressMessages(library(RMySQL))
suppressMessages(library(edgeR))
suppressMessages(library(limma))

argv       <- commandArgs(trailing = T)
tmpprefix  <- argv[1]
genes      <- argv[2]
muttype    <- argv[3]
cancerid   <- argv[4]
pvalue     <- as.numeric(argv[5])
foldchange <- as.numeric(argv[6])
genetable  <- argv[7]
plotnum    <- argv[8]
dbsrc      <- argv[9]
filtergene <- argv[10]
filterout  <- argv[11]

proc.time()
print("MESSAGE: Start")
# MySQL connection
con  <- dbConnect(MySQL(), user="XXXX", password="XXXX", dbname="mutarget", host="localhost")

# Expression matrix
count <- as.matrix(read.table(paste(cancerid, "tsv", sep = "."), check.names = F, sep = "\t"))
proc.time()
print("MESSAGE: Expression matrix")

# Column data
coldata <- data.frame(gene = factor(rep("WT", ncol(count)), levels = c("WT","Mut")))
rownames(coldata) <- colnames(count)

# Filtering samples
if(!is.na(filtergene)){
	if(filterout == "include"){
		query.samples <- paste("select distinct(submitid) as samples from individual inner join (mutation,genetable) on (individual_patientid = patientid and genetable_geneid = geneid) where cancer_cancerid = ",cancerid," and genename = '",filtergene,"';",sep="")
	} else {
		query.samples <- paste("select distinct(submitid) as samples from individual
						where cancer_cancerid = ",cancerid," and 
						name not in (
							     select distinct(name) from individual
							     	inner join (mutation,genetable) on (individual_patientid = patientid and genetable_geneid = geneid) 
							     where cancer_cancerid = ",cancerid," and
							     genename = '",filtergene,"'
							     );",sep="")
	}

	rs        <- dbSendQuery(con, query.samples)
	raw       <- fetch(rs, n=-1)
	coldata   <- coldata[rownames(coldata) %in% raw$samples,,drop=F]
	count     <- count[,colnames(count) %in% raw$samples]
}

# Get mutant samples for coldata
genes <- strsplit(genes, ",")
genes <- paste("genename = '",gsub(",","' or genename = '",genes),"'",sep="")
query <- paste("select distinct(name) as samples from individual inner join (mutation,genetable) on (individual_patientid = patientid and genetable_geneid = geneid) where cancer_cancerid = ",cancerid," and (",genes,") and muteffect_effectid = ",muttype,";",sep="")
rs    <- dbSendQuery(con, query)
raw   <- fetch(rs, n=-1)
index <- grep(paste(raw$samples, collapse="|"), rownames(coldata))
coldata$gene[index] <- "Mut"
if(length(coldata[coldata$gene == "Mut",]) == 0){
	print("MESSAGE: Not enough samples with alteration to split into cohorts")
	quit(save="no")
}
proc.time()
print("MESSAGE: Get mutant samples")

# Differential expression using edgeR
edge <- DGEList(counts = count, group = coldata$gene)
keep <- rowSums(cpm(edge) > 1) > 2
edge <- edge[keep, , keep.lib.sizes=FALSE]
edge <- calcNormFactors(edge)
edge <- estimateDisp(edge)
des  <- exactTest(edge)
des  <- as.data.frame(topTags(des, n = 32000, p.value = pvalue))
des  <- des[abs(des$logFC) > foldchange,]
# I need to remove log, because "Our users cannot understand it..."
# No further comment
des$logFC <- exp(des$logFC)
colnames(des)[1] <- "foldchange"

proc.time()
print("MESSAGE: Differential expression")

# If we would like to process Metabric data
if(cancerid == 1 && dbsrc == 'TCGAandMetabric'){
	metaexp <- as.matrix(read.table("data_expression.txt", sep = "\t", check.names = F))
	query   <- paste("select distinct(submitid) as samples from individual inner join (mutation, genetable) on (individual_patientid = patientid and genetable_geneid = geneid) where cancer_cancerid = ",cancerid," and (",genes,") and muteffect_effectid = ",muttype,";", sep = "")
	rs      <- dbSendQuery(con, query)
	raw     <- fetch(rs, n=-1)
	# Creating coldata
	design  <- rep("WT",ncol(metaexp))
	names(design) <- colnames(metaexp)
	design[names(design) %in% raw$samples] <- "Mut"
	# Check if there is two groups
	l <- length(design[design == "Mut"])
	if(l == 0 || l == length(design)){
		print("There is no two groups for Metabric")
	} else {
		mm <- model.matrix( ~0 + factor(design,levels=c("Mut","WT")))
		colnames(mm) <- c("Mut", "WT")
		# Limma differential expression
		fit <- lmFit(metaexp, mm)
		contr <- makeContrasts(Mut-WT, levels = mm)
		fit <- contrasts.fit(fit,contr)
		fit <- eBayes(fit)
		des2 <- topTable(fit, adjust.method="BH",p.value=pvalue,number=80000) #FIXME filtering by foldchange is missing
		if(nrow(des2) == 0) {
			print("MESSAGE: No significant result in Metabric analysis")
			quit(save="no")
		}
		# Removing log for stupid users
		write.table(des2, "/tmp/checkit.tsv", sep = "\t", quote = F)
		des2$logFC <- exp(des2$logFC)
		colnames(des2)[1] <- "foldchange"
		# Create a hybrid table with duplicated fold change
		commongenes <- intersect(rownames(des), rownames(des2))
		des <- des[commongenes,c(1,3,4)]
		des2 <- des2[commongenes,c(1,4,5)]
		des <- cbind(des,des2)
		colnames(des) <- c('TCGA foldchange', 'TCGA Pvalue', 'TCGA FDR', 'Metabric foldchange', 'Metabric Pvalue', 'Metabric FDR')
	}
	proc.time()
	print('MESSAGE: Processing Metabric')
}

# Filtering results
if(genetable != "all"){
	if(genetable == 'action'){
		query <- "select genename from genetable where actionable = 1;"
	} else {
		query <- "select genename from genetable where fda = 1;"
	}
	rs    <- dbSendQuery(con, query)
	raw   <- fetch(rs, n=-1)
	index <- grep(paste(raw$genename, collapse="|"), rownames(des))
	des   <- des[index,]
}

write.table(format(des, digits = 2), paste(tmpprefix, "tsv", sep = "."), quote = F, sep = "\t")

proc.time()
print("MESSAGE: Creating plots")
if(plotnum == "all" | as.numeric(plotnum) > nrow(des)){
	plotnum <- nrow(des)
} else {
	plotnum <- as.numeric(plotnum)
}

for(index in 1:plotnum){
   gene <- rownames(des)[index]
   png(paste(tmpprefix,gene,"png",sep="."))
   boxplot(cpm(edge)[gene, rownames(coldata[coldata$gene == "Mut",,drop=F])], cpm(edge)[gene,rownames(coldata[coldata$gene == "WT",,drop=F])], main = paste(gene, "expression", sep=" "), names = c("Mutant", "WT"))
   dev.off()
}

proc.time()
print("MESSAGE: End of script")
