
## ----, eval=FALSE--------------------------------------------------------
## install.packages("ggplot2")
## 
## source("http://bioconductor.org/biocLite.R")
## biocLite()
## biocLite("Biobase")
## biocLite("GEOquery")
## biocLite("limma")


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
## help(seq)
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
## # Run it without any arguments to install the core packages or update any installed packages.
## # This requires internet connectivity and will take some time!
## biocLite()


## ----, eval=FALSE--------------------------------------------------------
## # Do only once
## source("http://bioconductor.org/biocLite.R")
## biocLite("limma")
## 
## # Every time you need to use the limma package
## library(limma)


## ----, message=FALSE-----------------------------------------------------
# Bioconductor packages.
# Use the installation instructions at http://www.bioconductor.org/install/
# if you haven't already installed these (install once, load every time)
library(Biobase)
library(GEOquery)
library(limma)
library(arrayQualityMetrics)


## ----, eval=FALSE--------------------------------------------------------
## gset <- getGEO("GSE25724", GSEMatrix =TRUE, destdir=".", AnnotGPL=TRUE)
## gset <- gset[[1]]
## class(gset)


## ----, eval=FALSE--------------------------------------------------------
## # Removes spaces, colons, other weird symbols, inserts "." instead.
## fvarLabels(gset) <- make.names(fvarLabels(gset))


## ----, eval=FALSE--------------------------------------------------------
## # group names for all samples
## View(pData(gset))
## rep("ctrl", 7)
## rep("diab", 6)
## condition <- c(rep("ctrl", 7), rep("diab", 6))
## condition <- as.factor(condition)
## gset$description <- condition


## ----, eval=FALSE--------------------------------------------------------
## arrayQualityMetrics(gset, outdir="aqm-report", intgroup="description")


## ----, eval=FALSE--------------------------------------------------------
## # Set up and view the design matrix
## design <- model.matrix(~ description + 0, data=gset)
## colnames(design) <- levels(condition)
## design
## 
## # Set up and view the contrast matrix
## cont.matrix <- makeContrasts(diab-ctrl, levels=design)
## cont.matrix
## 
## # Fit a linear model, compute contrast coefficients and std errs, moderated t-tests and p-values:
## fit <- lmFit(gset, design)
## fit <- contrasts.fit(fit, cont.matrix)
## fit <- eBayes(fit, proportion=0.01)
## 
## # Extract the top 250 most differentially expressed genes, with an FDR correction
## tt  <- topTable(fit, adjust="fdr", sort.by="p", number=250)
## 
## # Look at the results in RStudio
## View(tt)
## 
## # Export the results to file
## write.table(tt, file='top-250-genes.csv', row.names=F, sep=",")


## ------------------------------------------------------------------------
sessionInfo()


