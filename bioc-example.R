library(Biobase)
library(GEOquery)
library(limma)

gset <- getGEO("GSE25724", GSEMatrix =TRUE, destdir=".", AnnotGPL=TRUE)
gset <- gset[[1]]

annotation(gset)

# make proper column names to match toptable 
# Removes spaces, colons, other weird symbols, inserts "." instead.
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
View(pData(gset))
rep("ctrl", 7)
rep("diab", 6)
condition <- c(rep("ctrl", 7), rep("diab", 6))
condition <- as.factor(condition)

# log2 transform
exprs(gset) <- log2(exprs(gset))

# set up the data and proceed with analysis
gset$description <- condition
design <- model.matrix(~ description + 0, data=gset)
colnames(design) <- levels(condition)
design
cont.matrix <- makeContrasts(diab-ctrl, levels=design)
cont.matrix

fit <- lmFit(gset, design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit, proportion=0.01)
tt  <- topTable(fit, adjust="fdr", sort.by="p", number=250)

View(tt)

write.table(tt, file='top-250-genes.csv', row.names=F, sep=",")

## Boxplot
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
palette(c("#f4dfdf","#f2cb98", "#AABBCC"))
boxplot(exprs(gset), boxwex=0.6, notch=TRUE, outline=FALSE, col=condition, las=2)
legend("topleft", labels, fill=palette(), bty="n")