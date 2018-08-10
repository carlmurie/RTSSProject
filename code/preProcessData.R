library(Biobase)

ROOT <- "/shared/silo_researcher/Gottardo_R/cmurie_working/RTSS/"

pheno <- read.csv(paste0(ROOT, "RTSS_Full/data/RTSSFullPheno_20180515.csv"))
dat <- read.table(paste0(ROOT, "RTSS_Full/data/allData.txt"), sep="\t")

## filtering cutoffs
geneCut <- 20
sampCut <- 75000

## separate controls from biological samples
controls <- setdiff(colnames(dat), pheno$colID)
conInd <- match(controls, colnames(dat))
bDat <- dat[,-conInd]

## ensure correct ordering
ordInd <- match( colnames(bDat), pheno$colID)
bPheno <- pheno[ordInd,]
if(sum(colnames(bDat) != bPheno$colID) != 0) {
  stop("ERROR: column names of data don't match annotation file")
}

### remove low read samples
totReads <- apply(bDat, 2, sum)

## get idea of low to mid sample reads
mInd <- totReads > sampCut & totReads <= 300000
mids <- sum(mInd)
mRan <- range(totReads[mInd])

sampInd <- totReads > sampCut
totReadsF <- uRan <- range(totReads[sampInd])
fDat <- bDat[,sampInd]

## filter genes with low read counts. Gene needs at least 20 samples with reads >= 5
summ <- apply(fDat, 1, function(x) sum(x >= 5))
summInd <- summ > geneCut

uDat <- fDat[summInd,]
uPheno <- bPheno[sampInd,]

## set up factors, relevel, rename
uPheno$age <- factor(uPheno$agec, levels=unique(uPheno$agec), labels=c("old", "young"))
uPheno$age  <- relevel(uPheno$age, "young")
uPheno$agec <- NULL
uPheno$stimulation <- relevel(uPheno$stimulation, "dmso")
uPheno$vaccine <- relevel(uPheno$vaccine, "comparator")
uPheno$case <- factor(uPheno$rna_seq_casecon.m12, levels=c("case", "control", ""),
                      labels=c("case", "control", "neither"))
uPheno$rna_seq_casecon.m12 <- NULL
uPheno$case <- relevel(uPheno$case, "control")
uPheno$match <- uPheno$rna_seq_match_id_m12
uPheno$rna_seq_match_id_m12 <- NULL
uPheno$TotalReads <- apply(uDat, 2, sum)
row.names(uPheno) <- uPheno$colID

counts_Eset <- ExpressionSet(as.matrix(uDat), phenoData=AnnotatedDataFrame(uPheno))
saveRDS(counts_Eset, paste0(ROOT, "RTSSProject/data/RawReadsFiltered.rds"))

