library("readr")
library("dplyr")
library("tximport")
library("edgeR")
library("RColorBrewer")
library("scico")

rm(list=ls())

dir <- "/home/rchavez/TexasTech/Mary/transcriptomes/transcriptome__drap_ALL__Mary_fastqs/final_tpm1_Mary/kallisto_Mary_fastqs_filtered/" # kallisto output directory
diro <- "/home/rchavez/TexasTech/Mary/manuscript_final/Wisecaver/edgeR/" # edgeR output directory
fdr <- 0.01

# metadata
file_metadata <- "/home/rchavez/TexasTech/Mary/files/metadata_pparvum.txt"
metadata <- read.delim(file_metadata)

sample_id <- list.dirs(path=dir, full.names=FALSE, recursive=FALSE)
files <- file.path(dir, sample_id, "abundance.tsv")
names(files) <- metadata[match(sample_id, metadata$sample), "replicate"] # optional: change to replicate2 for MDS plot
# sample glyph_20_r3 is outlier
#files <- files[!grepl("glyph_20_r3", names(files))]

txi.kallisto <- tximport(files, type="kallisto", txOut=TRUE, dropInfReps=TRUE)
cts <- txi.kallisto$counts
#colnames(cts) <- metadata$replicate[match(colnames(cts), metadata$fastq)]

groups <- factor(sub("_r\\d$", "", colnames(cts)))
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
#
normMat <- txi.kallisto$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
keep <- filterByExpr(y, design=design)
y <- y[keep, ]
#keep <- rowSums(cpm(y)>1) >= 2
#y <- y[keep, , keep.lib.sizes=FALSE]
y2 <- estimateDisp(y, design, robust=TRUE)

#plotBCV(estimateDisp(y, design, robust=TRUE, prior.df=0))
#plotBCV(y2)
file_out <- paste0(diro, "plot_MDS_dim12.svg")
svg(file_out)
plotMDS(y2, col=brewer.pal(nlevels(groups), 'Set1')[as.integer(groups)]) #plotMDS(y2, col=scico(nlevels(groups), palette = 'lapaz')[as.integer(groups)])
dev.off()
file_out <- paste0(diro, "plot_MDS_dim34.svg")
svg(file_out)
plotMDS(y2, col=brewer.pal(nlevels(groups), 'Set1')[as.integer(groups)], dim=c(3,4))
dev.off()

fit <- glmQLFit(y2, design, robust=TRUE)
plotQLDisp(fit)

contVector <- c("glyph_01" = "glyph_01-control",
				"glyph_20" = "glyph_20-control",
				"glyph_40" = "glyph_40-control")

contMatrix <- makeContrasts(contrasts=contVector, levels=design)
colnames(contMatrix) <- names(contVector)

r <- data.frame(filename=character(), contrast=character())
for (contrast in names(contVector))
{
	print(contrast)
	r[nrow(r)+1, ] <- c(contrast, contVector[contrast])
	#
	dif <- glmQLFTest(fit, contrast=contMatrix[ ,contrast])
	#
	tab <- as.data.frame(topTags(dif, n=nrow(dif), adjust.method="BH", sort.by="p.value"))
	tab$target_id <- rownames(tab)
	tab <- dplyr::select(tab, target_id, dplyr::everything())
	file_out <- paste0(diro, contrast, ".tsv")
	write_tsv(tab, file_out)
	#
	tab <- tab[tab$FDR <= fdr, ]
	file_out <- paste0(diro, contrast, "_fdr", sub("\\.", "", as.character(fdr)), ".tsv")
	if (nrow(tab) > 0)
	{
		write_tsv(tab, file_out)
	}
	else
	{
		writeLines("no_DEGs", file_out)
	}
}
write_tsv(r, paste0(diro, "contrasts.txt"))
