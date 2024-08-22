library("readr")
library("dplyr")
library("topGO")

rm(list=ls())

dir_edger <- "/castle/rchavez/Mary/edgeR/edgeR_Mary_fastqs_filtered/"
dir_go <- "/elpool/data/TAIR/GO_20220401/"
diro <- "/castle/rchavez/Mary/topGO/topGO_Mary_fastqs_filtered/"
cat <- "BP"

# This file_obo is optional and has to be created manually.
# It's a tsv containing GO ids and GO names, as topGO will
# truncate long GO names. If this file is not present just comment
# lines 52 and 53 below
file_obo <- paste0(dir_go, "obo_table_full.tsv")
obo <- read.delim(file_obo)

file_annot <- "/castle/rchavez/Mary/transcriptomes/transcriptome__drap_ALL__Mary_fastqs/final_tpm1_Mary/funcdesc/topGO/topGO_P_ppar.txt"
geneID2GO <- readMappings(file=file_annot)

file_genes <- "/castle/rchavez/Mary/transcriptomes/transcriptome__drap_ALL__Mary_fastqs/final_tpm1_Mary/loci_ppar.txt"
geneUniverse <- scan(file_genes, what="character", sep="\n")
length(geneUniverse)

files <- list.files(dir_edger, pattern="\\.tsv$", full.names=TRUE)
files <- files[grepl("fdr", files)]

updown <- "down"
for (file_query in files)
{
	#file_query <- files[1]
	query <- read.delim(file_query)
	#genesOfInterest <- query[, "target_id"]
	if (updown == "up"){genesOfInterest <- query[query$logFC > 0, "target_id"]}
	if (updown == "down"){genesOfInterest <- query[query$logFC < 0, "target_id"]}
	#
	name <- sub("\\.\\w{3}$", "", basename(file_query))
	#
	geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
	names(geneList) <- geneUniverse
	#
	myGOdata <- new("topGOdata", description=name, ontology=cat, nodeSize=10, allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
	#
	result <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
	#top20 <- GenTable(myGOdata, classicFisher=result, topNodes=20)
	allGO <- usedGO(object=myGOdata)
	allRes <- GenTable(myGOdata, classicFisher=result, topNodes=length(allGO))
	#
	colnames(allRes)[colnames(allRes) == "GO.ID"] <- "go_id"
	allRes$Term <- NULL
	allRes <- left_join(allRes, obo[, c("go_id", "go_name")], by="go_id")
	allRes <- dplyr::select(allRes, go_id, go_name, everything())
	#
	file_out <- paste0(diro, name, "_topGO_", cat, "_", updown, ".tsv")
	write_tsv(allRes, file_out)
}
