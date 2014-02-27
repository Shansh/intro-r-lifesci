## Introduction to R and Bioconductor for Cell Biology 8401

# Introduction to R -------------------------------------------------------

# Introduce R, Rstudio, layout, workspace, etc.
# R is a glorified calculator
# Start this in console, move to editor, use run button
2+2
5*4
2^3

# Source the file with echo

# R Knows order of operations
2+3*4

# built in functions
# functions
sqrt(144)
log(1000)

# Get help
?log

# Functions with arguments
log(1000) # note syntax highlighting
log(1000, base=10)
sqrt(log(1000, base=10)) #nested functions

# vectors
1:5
6:10
1:5 + 6:10 # vectorized operations
1:100 # the bracketed number in the gutter is just a counter of the number of elements
?seq
seq(from=2, to=200, by=2)

# Variables
x <- 5
x

y <- 42
y

y-x
z <- y-x
z

x <- 1:5
y <- 6:10
x
y
x+y
x*y
x^y

# See what's in environment by looking at Environment tab in Rstudio.
# Or run ls()
ls()
rm(x)
ls()
x
attrm(y,z)

# Data types
class(42)
class(log)
name <- "Stephen"
name
class(name)
toupper(name) # name is an object of class character. methods or functions are associated with certain classes.
toupper(log) # can't run a function that expects character on an object of class function

# combine values into vector
?c # built in function called c
x <- c(1,3,5)
x
class(x)
length(x)

y <- c("My", "name", "is", "Stephen")
y
class(y)
length(y)
y <- c(y, "Turner")
y
length(y)


x
sum(x)
y
z <- c(x,y)
z 
class(z) #combining characters with numbers results in coercing to character
sum(z)
?c

# subsetting
x <- 101:150
x
x[1]
x[50]
x[20:30]
x[40:60] #NA is missing value!


# data frames
# load motor trend car road tests from 1974
data(mtcars)
mtcars
head(mtcars)
length(mtcars) # the number of features
dim(mtcars)
dim(mtcars)[1] # number of cars (samples)
str(mtcars)

#access elements
mtcars$mpg
mtcars$cyl
# subset
head(mtcars)
mtcars[1:4, 1:2]
# see only 6 cylinder cars, etc
subset(mtcars, cyl==6)
subset(mtcars, cyl>6)
subset(mtcars, cyl==6, select=c(mpg, disp))

# with
mtcars$mpg / mtcars$cyl * mtcars$disp
with(mtcars, mpg/cyl*disp)

# plot
plot(mtcars$mpg)
with(mtcars, plot(mpg))
hist(mtcars$mpg)
hist(mtcars$mpg, breaks=10)
hist(mtcars$mpg, breaks=10, col="black")
hist(mtcars$mpg, breaks=10, col="gray50")
plot(mtcars$disp, mtcars$mpg)
with(mtcars, plot(disp, mpg))
with(mtcars, plot(disp, mpg, col="green"))
with(mtcars, plot(disp, mpg, pch=16))
with(mtcars, plot(disp, mpg, pch=16,  col="red"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement", 
                  ylab="Miles per Gallon", xlab="Displacement (cu. in.)"))

# exercise: plot horsepower vs displacement only for 8 cylinder vehicles. 
# give it a title, label the x and y axes
with(subset(mtcars, cyl==8), plot(disp, hp, pch=16, xlab="Displacement", ylab="Horsepower", main="HP vs Disp"))
with(mtcars, plot(disp, hp, pch=16, xlab="Displacement", ylab="Horsepower", main="HP vs Disp"))

# saving
mtcars_8cyl <- subset(mtcars, cyl==8)
mtcars_8cyl
?write.table
?write.csv
getwd()
# go to session->set working directory, change to desktop/r
getwd()
setwd("/Users/sdt5z/Desktop/R")
write.csv(mtcars_8cyl, file="cars8.csv")

## remove the dataset, load again using rstudio tools menu, or manual command
rm(mtcars_8cyl)
mtcars_8cyl
?read.csv
cars8 <- read.table(file="cars8.csv", header=TRUE, sep=",", row.names=1)


# BioC --------------------------------------------------------------------

## http://www.bioconductor.org/packages/2.11/bioc/vignettes/SPIA/inst/doc/SPIA.pdf
## Download data from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4107
## http://people.virginia.edu/~sdt5z/GSE4107_RAW.zip
## Extract tar file, gunzip all the .CEL.gz files

# Do this in folder on desktop
setwd("~/R")

biocLite(c("affy", "AnnotationDbi", "hgu133plus2cdf", "hgu133plus2.db", "genefilter", "DBI", "annotate", "arrayQualityMetrics", "limma", "GOstats", "Category", "GO.db", "KEGG.db")

?list.celfiles
library(affy)
?list.celfiles
list.celfiles("GSE4107_RAW", full.names=T)
myfiles <- list.celfiles("GSE4107_RAW", full.names=T)

?ReadAffy
myaffybatch <- ReadAffy(filenames=myfiles)
myaffybatch
class(myaffybatch)
?"AffyBatch-class"

?rma
eset <- rma(myaffybatch)
class(eset)
str(eset)
?"ExpressionSet-class"
head(exprs(eset))
dim(exprs(eset))

## Annotate the samples
pData(eset)
class(pData(eset))
pData(eset)$condition
rep("cancer", 12)
rep("healthy", 10)
c(rep("cancer", 12), rep("healthy", 10))
factor(c(rep("cancer", 12), rep("healthy", 10)))
pData(eset)$condition <- factor(c(rep("cancer", 12), rep("healthy", 10)))
pData(eset)


# filter the dataset
library(genefilter)
eset
?featureFilter
eset <- featureFilter(eset, require.entrez=T, remove.dupEntrez=T, feature.exclude="^AFFX")
eset

## Annotate the features
fData(eset)
annotation(eset)
library(annotate)
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ID     <- featureNames(eset)
class(ID)
head(ID)
head(exprs(eset))
?lookUp
#trust me on this one
Symbol <- as.character(lookUp(ID, "hgu133plus2.db", "SYMBOL"))
Name   <- as.character(lookUp(ID, "hgu133plus2.db", "GENENAME"))
Entrez <- as.character(lookUp(ID, "hgu133plus2.db", "ENTREZID"))
tmp <- data.frame(ID=ID, Entrez=Entrez, Symbol=Symbol, Name=Name, stringsAsFactors=F)
head(tmp)
fData(eset) <- tmp
head(fData(eset))
rm(ID, Symbol, Name, Entrez, tmp)

# Quality Assessment ------------------------------------------------------

# Load the arrayQualityMetrics package
library(arrayQualityMetrics)

# Run arrayQualityMetrics on your eset, giving it the name of an ouput directory.
# intgroup is the name of the sample covariate(s) used to draw a colour side bar next to the heatmap. character matching column in pData(eset)
# also give it a title
?arrayQualityMetrics
arrayQualityMetrics(eset, outdir="aqm", intgroup="condition", reporttitle="Quality Assessment Report")
# Now go open the current working directory, and open aqm/index.html in your browser


# Data Analysis -----------------------------------------------------------

## The limma package will be used for analysis
library(limma)

## Make the design and contrasts matrix
pData(eset)
design <- model.matrix(~0+condition, data=pData(eset))
design
class(design)
colnames(design)
colnames(design) <- c("cancer", "healthy")
design
colnames(design) <- sub("condition", "", colnames(design))
design

## make contrast matrix
contrast.matrix <- makeContrasts(cancer_v_healthy=cancer-healthy, levels=design)
contrast.matrix

## Fit model
fit <- lmFit(eset,design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
class(fit)

## Get statistics from that model fit
topTable(fit)
?topTable
nrow(fit)
tt <- topTable(fit, number=nrow(fit))

head(tt)
class(tt)
dim(tt)

## First, write the entire results table to a csv file and open in Excel
## Assignment: how many genes are significant at an adjusted p-value of 0.005?
## Write this table to a CSV file and open it in excel
write.csv(tt, file="diffexpr-all.csv")
sig <- subset(tt, adj.P.Val<0.005)
dim(sig)
write.csv(sig, file="diffexpr-significant.csv")

## volcano plot
with(tt, plot(logFC, -log10(adj.P.Val), pch=16, col="black"))

## heatmap
heatmap(exprs(eset[featureNames(eset) %in% sig$ID, ]))

library(GOstats)
library(Category)
?hyperGTest
?"GOHyperGParams-class"
myuniverse <- tt$Entrez
mysiggenes <- sig$Entrez
myanno <- annotation(eset)

?"KEGGHyperGParams-class"
library(KEGG.db)
KEGGparams <- new("KEGGHyperGParams", geneIds=mysiggenes, universeGeneIds=myuniverse, annotation=myanno, testDirection="over", pvalueCutoff=.05)
KEGGres <- hyperGTest(KEGGparams)
KEGGres
class(KEGGres)
?"KEGGHyperGResult-class"
?"HyperGResult-accessors"
summary(KEGGres)
htmlReport(KEGGres, file="KEGG report.html")

?"GOHyperGParams-class"
library(GO.db)
GOparams <- new("GOHyperGParams", geneIds=mysiggenes, universeGeneIds=myuniverse, annotation=myanno, conditional=TRUE, ontology="MF")
GOres <- hyperGTest(GOparams)
htmlReport(GOres, file="GO MF terms.html")
