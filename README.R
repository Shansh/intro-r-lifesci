
## ----, results='hide'----------------------------------------------------
2+2
5*4
2^3


## ----, results='hide'----------------------------------------------------
2+3*4/(5+3)*15/2^2+3*4^2


## ----, results='hide'----------------------------------------------------
# Notice that this is a comment. 
# Anything behind a # is "commented out" and is not run.
sqrt(144)
log(1000)


## ----, results='hide'----------------------------------------------------
log(1000)
log(1000, base=10)
sqrt(log(1000, base=10))


## ----, eval=FALSE--------------------------------------------------------
## # Some simple numeric vectors:
## 1:5
## 6:10
## 1:5 + 6:10
## 1:100
## 
## # Get some help with the seq() function, then create a vector from 2 to 200 by 2s.
## # Notice how the seq() function works -- the `to` argument will never be exceeded.
## ?seq
## seq(from=2, to=200, by=4)


## ----, results='hide'----------------------------------------------------
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


## ----, results='hide'----------------------------------------------------
ls()
rm(x)
ls()
x # oops! you should get an error because x no longer exists!
rm(y,z)


## ----, results='hide'----------------------------------------------------
class(42)
class(log)
name <- "Stephen"
name
class(name)


## ------------------------------------------------------------------------
toupper(name) # name is an object of class character. methods or functions are associated with certain classes.
toupper(log) # can't run a function that expects character on an object of class function


## ----, results='hide'----------------------------------------------------
# Get some help with ?c
x <- c(1,3,5)
x
class(x)
length(x)


## ----, results='hide'----------------------------------------------------
y <- c("My", "name", "is", "Stephen")
y
class(y)
length(y)
y <- c(y, "Turner")
y
length(y)


## ----, results='hide'----------------------------------------------------
sum(x)


## ----, results='hide'----------------------------------------------------
y
z <- c(x,y)
z 
class(z)


## ------------------------------------------------------------------------
z
sum(z)


## ----, results='hide'----------------------------------------------------
# Create the vector.
x <- 101:150

# Get the first element.
x[1]

# Get the 42nd element.
x[42]

# Get the 20th through the 25th elements. 
x[20:25]

# If you try to access elements that don't exist, you'll return missing values.
# Missing values are represented as NA
x[45:55] #NA is missing value!


## ----, results='hide'----------------------------------------------------
data(mtcars)
class(mtcars)
mtcars


## ----, results='hide'----------------------------------------------------
head(mtcars)
length(mtcars)
dim(mtcars)
dim(mtcars)[1] # number of rows (individual cars in the survey)
dim(mtcars)[2] # number of columns (number of variables measured)
str(mtcars)


## ----, results='hide'----------------------------------------------------
# display the number of cylinders for each car. 
mtcars$cyl
# first display MPG for all vehicles, then calculate the average.
mtcars$mpg
mean(mtcars$mpg)


## ----, results='hide'----------------------------------------------------
head(mtcars)
mtcars[1:4, 1:2]


## ----, results='hide'----------------------------------------------------
subset(mtcars, cyl==6)
subset(mtcars, cyl>6)
subset(mtcars, mpg>=20 | disp<100)
subset(mtcars, cyl==6, select=c(mpg, disp))
subset(mtcars, cyl>=6 & mpg>=15, select=c(mpg, cyl, qsec))


## ----, results='hide'----------------------------------------------------
# Display the number of cylinders.
mtcars$cyl
with(mtcars, cyl)

# Compute the senseless value described above. Both return the same results.
mtcars$mpg * mtcars$cyl / mtcars$disp
with(mtcars, mpg*cyl/disp)


## ------------------------------------------------------------------------
plot(mtcars$mpg)


## ------------------------------------------------------------------------
hist(mtcars$mpg)
hist(mtcars$mpg, breaks=10)
hist(mtcars$mpg, breaks=10, col="black")


## ------------------------------------------------------------------------
# This would also work, but let's use with().  
# plot(mtcars$disp, mtcars$mpg)
with(mtcars, plot(disp, mpg))


## ------------------------------------------------------------------------
with(mtcars, plot(disp, mpg, pch=16))
with(mtcars, plot(disp, mpg, pch=16,  col="red"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement", 
                  ylab="Fuel Economy (MPG)", xlab="Displacement (cu. in.)"))


## ----, echo=FALSE--------------------------------------------------------
with(subset(mtcars, cyl>4), plot(disp, hp, pch=16, col="blue",
                                 xlab="Displacement (cu. in.)", ylab="Gross Horsepower", 
                                 main="Horsepower vs Displacement for 6 and 8-cylinder vehicles"))


## ----, results='hide'----------------------------------------------------
mtcars_8cyl <- subset(mtcars, cyl==8)
mtcars_8cyl


## ----, eval=FALSE--------------------------------------------------------
## getwd()
## setwd("/Users/sdt5z/Desktop/R")


## ----, eval=FALSE--------------------------------------------------------
## write.csv(mtcars_8cyl, file="cars8.csv")


## ----, eval=FALSE--------------------------------------------------------
## rm(mtcars_8cyl)
## mtcars_8cyl
## cars8 <- read.table(file="cars8.csv", header=TRUE, sep=",", row.names=1)
## cars8
## rm(cars8)
## cars8 <- read.csv(file="cars8.csv", header=TRUE, row.names=1)
## cars8


## ----, eval=FALSE--------------------------------------------------------
## # Install only once.
## install.packages("ggplot2")
## 
## # Load the package every time you want to use it.
## library(ggplot2)


## ----, eval=FALSE--------------------------------------------------------
## # Download the installer script
## source("http://bioconductor.org/biocLite.R")
## 
## # biocLite() is the bioconductor installer function.
## # run it without any arguments to install the core packages or update any installed packages.
## biocLite()


## ----, eval=FALSE--------------------------------------------------------
## # Do only once
## source("http://bioconductor.org/biocLite.R")
## biocLite("limma")
## 
## # Every time you need to use the limma package
## library(limma)


## ----, eval=FALSE, echo=FALSE--------------------------------------------
## # BioC --------------------------------------------------------------------
## 
## ## http://www.bioconductor.org/packages/2.11/bioc/vignettes/SPIA/inst/doc/SPIA.pdf
## ## Download data from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4107
## ## http://people.virginia.edu/~sdt5z/GSE4107_RAW.zip
## ## Extract tar file, gunzip all the .CEL.gz files
## 
## # Do this in folder on desktop
## setwd("~/R")
## 
## biocLite(c("affy", "AnnotationDbi", "hgu133plus2cdf", "hgu133plus2.db", "genefilter", "DBI", "annotate", "arrayQualityMetrics", "limma", "GOstats", "Category", "GO.db", "KEGG.db")
## 
## ?list.celfiles
## library(affy)
## ?list.celfiles
## list.celfiles("GSE4107_RAW", full.names=T)
## myfiles <- list.celfiles("GSE4107_RAW", full.names=T)
## 
## ?ReadAffy
## myaffybatch <- ReadAffy(filenames=myfiles)
## myaffybatch
## class(myaffybatch)
## ?"AffyBatch-class"
## 
## ?rma
## eset <- rma(myaffybatch)
## class(eset)
## str(eset)
## ?"ExpressionSet-class"
## head(exprs(eset))
## dim(exprs(eset))
## 
## ## Annotate the samples
## pData(eset)
## class(pData(eset))
## pData(eset)$condition
## rep("cancer", 12)
## rep("healthy", 10)
## c(rep("cancer", 12), rep("healthy", 10))
## factor(c(rep("cancer", 12), rep("healthy", 10)))
## pData(eset)$condition <- factor(c(rep("cancer", 12), rep("healthy", 10)))
## pData(eset)
## 
## 
## # filter the dataset
## library(genefilter)
## eset
## ?featureFilter
## eset <- featureFilter(eset, require.entrez=T, remove.dupEntrez=T, feature.exclude="^AFFX")
## eset
## 
## ## Annotate the features
## fData(eset)
## annotation(eset)
## library(annotate)
## library(hgu133plus2.db)
## ls("package:hgu133plus2.db")
## ID     <- featureNames(eset)
## class(ID)
## head(ID)
## head(exprs(eset))
## ?lookUp
## #trust me on this one
## Symbol <- as.character(lookUp(ID, "hgu133plus2.db", "SYMBOL"))
## Name   <- as.character(lookUp(ID, "hgu133plus2.db", "GENENAME"))
## Entrez <- as.character(lookUp(ID, "hgu133plus2.db", "ENTREZID"))
## tmp <- data.frame(ID=ID, Entrez=Entrez, Symbol=Symbol, Name=Name, stringsAsFactors=F)
## head(tmp)
## fData(eset) <- tmp
## head(fData(eset))
## rm(ID, Symbol, Name, Entrez, tmp)
## 
## # Quality Assessment ------------------------------------------------------
## 
## # Load the arrayQualityMetrics package
## library(arrayQualityMetrics)
## 
## # Run arrayQualityMetrics on your eset, giving it the name of an ouput directory.
## # intgroup is the name of the sample covariate(s) used to draw a colour side bar next to the heatmap. character matching column in pData(eset)
## # also give it a title
## ?arrayQualityMetrics
## arrayQualityMetrics(eset, outdir="aqm", intgroup="condition", reporttitle="Quality Assessment Report")
## # Now go open the current working directory, and open aqm/index.html in your browser
## 
## 
## # Data Analysis -----------------------------------------------------------
## 
## ## The limma package will be used for analysis
## library(limma)
## 
## ## Make the design and contrasts matrix
## pData(eset)
## design <- model.matrix(~0+condition, data=pData(eset))
## design
## class(design)
## colnames(design)
## colnames(design) <- c("cancer", "healthy")
## design
## colnames(design) <- sub("condition", "", colnames(design))
## design
## 
## ## make contrast matrix
## contrast.matrix <- makeContrasts(cancer_v_healthy=cancer-healthy, levels=design)
## contrast.matrix
## 
## ## Fit model
## fit <- lmFit(eset,design)
## fit <- contrasts.fit(fit, contrast.matrix)
## fit <- eBayes(fit)
## class(fit)
## 
## ## Get statistics from that model fit
## topTable(fit)
## ?topTable
## nrow(fit)
## tt <- topTable(fit, number=nrow(fit))
## 
## head(tt)
## class(tt)
## dim(tt)
## 
## ## First, write the entire results table to a csv file and open in Excel
## ## Assignment: how many genes are significant at an adjusted p-value of 0.005?
## ## Write this table to a CSV file and open it in excel
## write.csv(tt, file="diffexpr-all.csv")
## sig <- subset(tt, adj.P.Val<0.005)
## dim(sig)
## write.csv(sig, file="diffexpr-significant.csv")
## 
## ## volcano plot
## with(tt, plot(logFC, -log10(adj.P.Val), pch=16, col="black"))
## 
## ## heatmap
## heatmap(exprs(eset[featureNames(eset) %in% sig$ID, ]))
## 
## library(GOstats)
## library(Category)
## ?hyperGTest
## ?"GOHyperGParams-class"
## myuniverse <- tt$Entrez
## mysiggenes <- sig$Entrez
## myanno <- annotation(eset)
## 
## ?"KEGGHyperGParams-class"
## library(KEGG.db)
## KEGGparams <- new("KEGGHyperGParams", geneIds=mysiggenes, universeGeneIds=myuniverse, annotation=myanno, testDirection="over", pvalueCutoff=.05)
## KEGGres <- hyperGTest(KEGGparams)
## KEGGres
## class(KEGGres)
## ?"KEGGHyperGResult-class"
## ?"HyperGResult-accessors"
## summary(KEGGres)
## htmlReport(KEGGres, file="KEGG report.html")
## 
## ?"GOHyperGParams-class"
## library(GO.db)
## GOparams <- new("GOHyperGParams", geneIds=mysiggenes, universeGeneIds=myuniverse, annotation=myanno, conditional=TRUE, ontology="MF")
## GOres <- hyperGTest(GOparams)
## htmlReport(GOres, file="GO MF terms.html")


