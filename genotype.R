suppressMessages(library(RMySQL))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(limma))

source("common.R")

argv       <- commandArgs(trailing = T)
tmpprefix  <- argv[1]
genes      <- argv[2]
muttype    <- argv[3] #FIXME This is a comma separated list
cancerid   <- argv[4]
pvalue     <- as.numeric(argv[5])
foldchange <- as.numeric(argv[6])
genetable  <- argv[7]
plotnum    <- argv[8]
dbsrc      <- argv[9]
diffexp    <- argv[10]
outliers   <- argv[11]
filtergene <- argv[12]
filterout  <- argv[13]

proc.time()
print("MESSAGE: Start")

# Expression matrix
count <- getExpMatrix(con, cancerid, dbsrc)
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

	raw       <- fetchDB(con, query.samples);
	coldata   <- coldata[rownames(coldata) %in% raw$samples,,drop=F]
	count     <- count[,colnames(count) %in% raw$samples]
}

# Get mutant samples for coldata
genes   <- paste("genename in ('",gsub(",","','",genes),"')",sep="")
muttype <- paste("muteffect_effectid in (",muttype,")",sep="")
query   <- paste("select distinct(name) as samples from individual inner join (mutation,genetable) on (individual_patientid = patientid and genetable_geneid = geneid) where cancer_cancerid = ",cancerid," and ",genes," and ",muttype,";",sep="")
raw     <- fetchDB(con, query)
index   <- grep(paste(raw$samples, collapse="|"), rownames(coldata))
coldata$gene[index] <- "Mut"
if(length(coldata[coldata$gene == "Mut",]) == 0){
	print("MESSAGE: Not enough samples with alteration to split into cohorts")
	quit(save="no")
}
proc.time()
print("MESSAGE: Get mutant samples")

# Differential expression using edgeR
if(dbsrc != 2){
	# rna-seq
	if(diffexp == "edgeR"){
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
		des  <- des[,c(1,3,4)]
		normexp <- cpm(edge)
	} else if(diffexp == "DESeq2") {
		des <- DESeqDataSetFromMatrix(count, colData = coldata, design =~gene)
		des <- DESeq(des)
		normexp <- assay(vst(des)) # Later we will overwrite variable des
		des <- results(des, contrast = c("gene", "Mut", "WT"))
		des <- des[!is.na(des$padj) & des$padj < pvalue & abs(des$log2FoldChange) > foldchange,]
		# Removing logarithm
		des$log2FoldChange <- exp(des$log2FoldChange)
		colnames(des)[2] <- "foldchange"
		des <- as.matrix(des)[,c(2,5,6)]
	}
} else {
	# microarray
	mm <- model.matrix( ~0 + coldata$gene)
	colnames(mm) <- c("Mut", "WT")
	fit <- lmFit(count, mm)
	contr <- makeContrasts(Mut-WT, levels = mm)
	fit <- contrasts.fit(fit, contr)
	fit <- eBayes(fit)
	des <- topTable(fit, adjust.method="BH",p.value=pvalue,number=80000) #FIXME filtering by foldchange is missing
	des$logFC <- exp(des$logFC)
	colnames(des)[1] <- "foldchange"
	normexp <- count # microarray already normalised
}

proc.time()
print("MESSAGE: Differential expression")

# Filtering results
if(genetable != "all"){
	if(genetable == 'action'){
		query <- "select genename from genetable where actionable = 1;"
	} else {
		query <- "select genename from genetable where fda = 1;"
	}
	raw   <- fetchDB(con, query)
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
outl <- FALSE
if(outliers == 1){
	outl = TRUE
}

for(index in 1:plotnum){
   gene <- rownames(des)[index]
   png(paste(tmpprefix,gene,"png",sep="."))
   boxplot(normexp[gene, rownames(coldata[coldata$gene == "Mut",,drop=F])], normexp[gene,rownames(coldata[coldata$gene == "WT",,drop=F])], main = paste(gene, "expression", sep=" "), names = c("Mutant", "WT"), outline = outl)
   dev.off()
}

proc.time()
print("MESSAGE: End of script")
