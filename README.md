# Introduction to R for Life Scientists

This workshop is directed toward life scientists with little to no experience with statistical computing or bioinformatics. This interactive workshop will introduce the R statistical computing environment, including basic instruction in data types, variables, array manipulation, functions, data frames, data import/export, visualization, and using packages. At the end of the workshop, participants will see a live demonstration of a real biomedical application - analysis of publicly available gene expression microarray data. This will demo (1) how to search for and acquire publicly accessible data from NCBI Gene Expression Omnibus, and (2) how to use Bioconductor packages to import, process, QC, analyze, and visualize the results of the analysis. By the end of the workshop, participants will be able to use R for basic data manipulation and visualization, and will know where to look for further help and instruction. Participants will also be exposed to downloading and analyzing publicly available gene expression data. An advanced follow-on course will go through a gene expression data analysis in detail. 

Link to slides: *coming soon.*

## Before coming

Prior to the workshop, please download the software as below and take the pre-workshop survey.

### Software setup

0. Download and install R for [Mac OS X](http://cran.rstudio.com/bin/macosx/R-3.0.3.pkg) or [Windows](http://cran.rstudio.com/bin/windows/base/R-3.0.3-win.exe)
0. Download and install RStudio Desktop: <http://www.rstudio.com/ide/download/desktop>.
0. Run RStudio, and enter the following commands into the "Console" window (usually the lower-right, by default). This will download and install necessary add-on packages we will use in class.


```coffee
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biobase")
biocLite("GEOquery")
biocLite("limma")
biocLite("arrayQualityMetrics")
```


### Take the pre-workshop survey

Please complete [this very short survey](https://docs.google.com/forms/d/1Ef4r-5yTOZO-rMGyjZ5M-wP3Q_j2WPtkp1HM_ksApnw/viewform) (it should take you less than 60 seconds). 

## R basics

The first part of this workshop will demonstrate very basic functionality in r, including functions, functions, vectors, creating variables, getting help, subsetting, data frames, plotting, and reading/writing files.

### Basic operations

R can be used as a glorified calculator. Try typing this in directly into the console. Then start typing this into the editor, and save your script. Use the run button, or press `CMD`+`Enter` (`Ctrl`+`Enter` on Windows).


```coffee
2 + 2
5 * 4
2^3
```


R Knows order of operations.


```coffee
2 + 3 * 4/(5 + 3) * 15/2^2 + 3 * 4^2
```


### Functions

R has built-in functions.


```coffee
# Notice that this is a comment.  Anything behind a # is 'commented out' and
# is not run.
sqrt(144)
log(1000)
```


Get help by typing a question mark in front of the function's name, or `help(functionname)`:

```
help(log)
?log
```

Note syntax highlighting when typing this into the editor. Also note how we pass *arguments* to functions. Finally, see how you can *next* one function inside of another (here taking the square root of the log-base-10 of 1000).


```coffee
log(1000)
log(1000, base = 10)
sqrt(log(1000, base = 10))
```


### Vectors

Let's create some numeric vectors. Vectors (aka "arrays" in Perl, "lists" in Python) are single *objects* containing an ordered collection of *elements*. A simple vector is a numeric vector, a single *object* containing several numbers. Here let's display a few vectors. We can also do vector arithmetic. When printing vectors to the screen that have lots of elements, notice that the bracketed number in the gutter of the output is just a counter indexing the number of elements in the vector.


```coffee
# Some simple numeric vectors:
1:5
6:10
1:5 + 6:10
1:100

# Get some help with the seq() function, then create a vector from 2 to 200
# by 2s.  Notice how the seq() function works -- the `to` argument will
# never be exceeded.
help(seq)
seq(from = 2, to = 200, by = 4)
```


### Variables / objects

You can store values in a variable or object. Use the `<-` operator for assignment. `=` also will work, but `<-` is conventional and preferred. Objects should start with a letter and can include numbers and underscores. Named variables are objects containing whatever is assigned to them. Remember, *everything is an object*.


```coffee
x <- 5
x

y <- 42
y

y - x
z <- y - x
z

x <- 1:5
y <- 6:10
x
y
x + y
x * y
x^y
```


You can see what objects (variables) are stored by viewing the Environment tab in Rstudio. You can also use the `ls()` function. You can remove objects (variables) with the `rm()` function. You can do this one at a time or remove several objects at once.


```coffee
ls()
rm(x)
ls()
x  # oops! you should get an error because x no longer exists!
```

```
## Error: object 'x' not found
```

```coffee
rm(y, z)
```


### Classes: everything is an object

Use the `class()` function to see what *kind* of object a something is. You can run `class()` on constants, built-in objects, or objects you create. Let's create a character object and then get it's class. 


```coffee
class(42)
class(log)
name <- "Stephen"
name
class(name)
```


Certain *functions* operate only on certain *classes* of object. Here, `name` is a `character` class, assigned to `"Stephen"`. The built-in `toupper()` function will operate on character objects, but not others. 


```coffee
toupper(name)  # name is an object of class character. methods or functions are associated with certain classes.
```

```
## [1] "STEPHEN"
```

```coffee
toupper(log)  # can't run a function that expects character on an object of class function
```

```
## Error: cannot coerce type 'special' to vector of type 'character'
```


We can combine values into a vector with the built-in `c()` function.


```coffee
# Get some help with ?c
x <- c(1, 3, 5)
x
class(x)
length(x)
```


Let's create and manipulate a character vector:


```coffee
y <- c("My", "name", "is", "Stephen")
y
class(y)
length(y)
y <- c(y, "Turner")
y
length(y)
```


Try running the built-in `sum()` function on a numeric vector:


```coffee
sum(x)
```


Combining characters with numerics results in coercing everything to be a character class. 


```coffee
y
z <- c(x, y)
z
class(z)
```


Certain functions only operate on certain classes. You can't compute the `sum()` of a character vector!


```coffee
z
```

```
## [1] "1"       "3"       "5"       "My"      "name"    "is"      "Stephen"
## [8] "Turner"
```

```coffee
sum(z)
```

```
## Error: invalid 'type' (character) of argument
```


### Slicing/indexing vectors

Let's create a vector of 50 integers going from 101 to 150. We can access certain elements of that vector by putting the element's *index(es)* in square brackets. E.g., `x[1]` will return the first element in vector `x`. Calling `x[c(3,5)]` will access the third and fifth elements. Calling `x[1:10]` will return the first ten elements of `x`. 

*Special note: R indexes vectors starting at 1. This is different from many other languages, including Perl and Python, which index starting from zero.*


```coffee
# Create the vector.
x <- 101:150

# Get the first element.
x[1]

# Get the 42nd element.
x[42]

# Get the 20th through the 25th elements.
x[20:25]

# If you try to access elements that don't exist, you'll return missing
# values.  Missing values are represented as NA
x[45:55]  #NA is missing value!
```


### Data Frames

Data frames are a standard way to store heterogeneous tabular data in R: tabular, meaning that individuals or observations are typically represented in rows, while variables or features are represented as columns; heterogeneous, meaning that columns/features/variables can be different classes (on variable, e.g. age, can be numeric, while another, e.g., cause of death, can be text).

Later on we'll go over how we load our own data, but for now, let's use a built-in data frame called `mtcars`. This data was extracted from the 1974 Motor Trend US magazine, and comprises fuel consumption and 10 aspects of automobile design and performance for 32 automobiles (1973–74 models). We can load this built-in data with `data(mtcars)`. By the way, running `data()` without any arguments will list all the available built-in datasets included with R.

Let's load the data first. Type the name of the object itself (`mtcars`) to view the entire data frame. *Note: doing this with large data frames can cause you trouble.*


```coffee
data(mtcars)
class(mtcars)
mtcars
```


There are several built-in functions that are useful for working with data frames.
* `head()` prints the first few lines of a large data frame.
* `length()` tells you the number of features (variables, columns) in a data frame.
* `dim()` returns a two-element vector containing the number of rows and the number of columns in a data frame.
* `str()` displays the structure of a data frame, printing out details and a preview of every column.
* `summary()` works differently depending on what kind of object you pass to it. Passing a data frame to the `summary()` function prints out some summary statistics about each column (min, max, median, mean, etc.)


```coffee
head(mtcars)
length(mtcars)
dim(mtcars)
dim(mtcars)[1]  # number of rows (individual cars in the survey)
dim(mtcars)[2]  # number of columns (number of variables measured)
str(mtcars)
```


We can access individual variables within a data frame using the `$` operator, e.g., `mydataframe$specificVariable`. Let's print out the number of cylinders for every car, and calculate the average miles per gallon for ever car in the dataset (using the built-in `mean()` function).


```coffee
# display the number of cylinders for each car.
mtcars$cyl
# first display MPG for all vehicles, then calculate the average.
mtcars$mpg
mean(mtcars$mpg)
```


We can also access certain rows or columns of a dataset by providing multiple indices using the syntax `mydataframe[rows, columns]`. Let's get the first 4 rows and the first two rows (MPG and # cylinders) from the dataset:


```coffee
head(mtcars)
mtcars[1:4, 1:2]
```


We can also use the `subset()` function to return a subset of the data frame that meets a specific condition. The first argument is the data frame you want to subset. The second argument is a condition you must satisfy. If you want to satisfy *all* of multiple conditions, you can use the "and" operator, `&`. The "or" operator `|` (the pipe character, usually shift-backslash) will return a subset that meet *any* of the conditions. 

The commands below will:

0. Return only cars with 6 cylinder engines.
0. Return only cars with greater than 6 cylinders.
0. Return only the cars that get at least 20 miles per gallon or have a displacement volume of less than 100cc.
0. Return cars with 6 cylinder engines, but using the `select=` argument, only the MPG and displacement columns. Note the syntax there -- we're passing a vector of variables created with the `c()` function to the `select=` argument, which only returns certain columns. 
0. Return cars that have greater than or equal to 6 cylinders *and* get at least 15 miles per gallon, but display only the MPG, cylinders, and qsec columns (qsec is the 1/4 mile time).

Try some subsetting on your own.


```coffee
subset(mtcars, cyl == 6)
subset(mtcars, cyl > 6)
subset(mtcars, mpg >= 20 | disp < 100)
subset(mtcars, cyl == 6, select = c(mpg, disp))
subset(mtcars, cyl >= 6 & mpg >= 15, select = c(mpg, cyl, qsec))
```


The `with()` function is particularly helpful. Let's say you wanted to compute some (senseless) value by computing the MPG times the number of cylinders divided by the car's displacement. You could access the dataset's variables using the `$` notation, or you could use `with()` to temporarily *attach* the data frame, and call the variables directly. The first argument to `with()` is the name of the data frame, and the second argument is all the stuff you'd like to do with the particular features in that data frame.

Try typing the following commands:


```coffee
# Display the number of cylinders.
mtcars$cyl
with(mtcars, cyl)

# Compute the senseless value described above. Both return the same results.
mtcars$mpg * mtcars$cyl/mtcars$disp
with(mtcars, mpg * cyl/disp)
```


### Plotting

Plotting a single numeric variable goes down the rows and plots a value on the y-axis for each observation (index) in the data frame. 


```coffee
plot(mtcars$mpg)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 


This isn't a very useful figure. More appropriate might be a histogram. We can try to let R decide how many breaks to insert in the histogram, or we can set that manually. We can also set the color of the bars. 



```coffee
hist(mtcars$mpg)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-241.png) 

```coffee
hist(mtcars$mpg, breaks = 10)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-242.png) 

```coffee
hist(mtcars$mpg, breaks = 10, col = "black")
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-243.png) 


We can create a scatterplot between two variables with `plot(varX, varY)`.


```coffee
# This would also work, but let's use with().  plot(mtcars$disp, mtcars$mpg)
with(mtcars, plot(disp, mpg))
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25.png) 


There are hundreds of plotting parameters you can use to make your plot look exactly like you want. Let's use a solid-filled point instead of an open circle with the `pch=` argument, color the points red with the `col=` argument, give it a title by passing a character object to the `main=` argument, and change the x and y axis titles with the `xlab=` and `ylab=` arguments, respectively. Let's go through this one step at a time. 


```coffee
with(mtcars, plot(disp, mpg, pch = 16))
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-261.png) 

```coffee
with(mtcars, plot(disp, mpg, pch = 16, col = "red"))
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-262.png) 

```coffee
with(mtcars, plot(disp, mpg, pch = 16, col = "red", main = "MPG vs Displacement"))
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-263.png) 

```coffee
with(mtcars, plot(disp, mpg, pch = 16, col = "red", main = "MPG vs Displacement", 
    ylab = "Fuel Economy (MPG)", xlab = "Displacement (cu. in.)"))
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-264.png) 


Notice how on that last line I broke the command up into two lines for better readability. I broke the command at the comma separating arguments, and indented the following line for readability.

On your own, try plotting horsepower vs displacement for vehicles with more than 4 cylinders. Give the graph a title and label the axes. Make the points solid (hint, `pch=16`) blue (hint, `col="blue"`) circles. Your plot should look something like this:

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27.png) 


### Reading in / writing out data

First, lets create a small dataset consisting of only 8 cylinder cars. 


```coffee
mtcars_8cyl <- subset(mtcars, cyl == 8)
mtcars_8cyl
```


Next, check what your working directory is with `getwd()` with no arguments, and look up some help for `write.table()` and `write.csv()`. 


```coffee
getwd()
help(write.table)
help(write.csv)
```


Using RStudio, go to the Session menu, and select the directory (folder) you want to work from under the "Set Working Directory" menu. You can also do this manually with the `setwd()` command. 


```coffee
getwd()
setwd("/Users/sdt5z/Desktop/R")
```


Once you've set your working directory either using RStudio or on the command line, save the new reduced data frame to a comma-separated file called `cars8.csv` using the `write.csv()` function. 


```coffee
write.csv(mtcars_8cyl, file = "cars8.csv")
```


Data can be loaded using the Tools -- Import Dataset -- From text file menu in R studio. Or you can also load a dataset manually using `read.table()` or `read.csv()`. First, read the help on these functions:


```coffee
help(read.table)
help(read.csv)
```


Here let's remove the dataset, and re-import it into an object called cars8 from the file we just saved.


```coffee
rm(mtcars_8cyl)
mtcars_8cyl
cars8 <- read.table(file = "cars8.csv", header = TRUE, sep = ",", row.names = 1)
cars8
rm(cars8)
cars8 <- read.csv(file = "cars8.csv", header = TRUE, row.names = 1)
cars8
```


## Using R and Bioconductor for gene expression analysis

In this section we'll analyze some publicly available gene expression data using R and bioinformatics-focused R packages in [Bioconductor](http://bioconductor.org/). But first, a bit about R *packages*. 

### Packages

Most generic R packages are hosted on the Comprehensive R Archive Network (CRAN, <http://cran.us.r-project.org/>). To install one of these packages, you would use `install.packages("packagename")`. You only need to install a package once, then load it each time using `library(packagename)`. Let's install the `ggplot2` package, and load it.


```coffee
# Install only once.
install.packages("ggplot2")

# Load the package every time you want to use it.
library(ggplot2)
```


### Bioconductor

Bioconductor packages work a bit different, and are not hosted on CRAN. Go to <http://bioconductor.org/> to learn more about the Bioconductor project. To use any Bioconductor package, you'll need a few "core" Bioconductor packages. Run the following commands to (1) download the installer script, and (2) install some core Bioconductor packages. You'll need internet connectivity to do this, and it'll take a few minutes, but it only needs to be done once.


```coffee
# Download the installer script
source("http://bioconductor.org/biocLite.R")

# biocLite() is the bioconductor installer function.  Run it without any
# arguments to install the core packages or update any installed packages.
# This requires internet connectivity and will take some time!
biocLite()
```


To install specific packages, first download the installer script if you haven't done so, and use `biocLite("packagename")`. This only needs to be done once then you can load the package like any other package. Let's download the [limma package](http://www.bioconductor.org/packages/release/bioc/html/limma.html): 


```coffee
# Do only once
source("http://bioconductor.org/biocLite.R")
biocLite("limma")

# Every time you need to use the limma package
library(limma)
```


Bioconductor packages usually have great documentation in the form of *vignettes*. For a great example, take a look at the [limma user's guide](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf). 

### Analyzing publicly available microarray data

Now, let's analyze some publicly available gene expression microarray data. NCBI Gene Expression Omnibus (<http://www.ncbi.nlm.nih.gov/geo/>) is an international public repository that archives and freely distributes microarray, next-generation sequencing, and other forms of high-throughput functional genomics data submitted by the research community. Many publishers require gene expression data be submitted to GEO and made publicly available before publication. You can learn a lot more about GEO by reading their [overview](http://www.ncbi.nlm.nih.gov/geo/info/overview.html) and [FAQ](http://www.ncbi.nlm.nih.gov/geo/info/faq.html) pages. At the time of this writing, GEO hosts over 45,000 studies comprising over 1,000,000 samples on over 10,000 different technology platforms.

In this demonstration, we're going to be using data from GEO Series accession number GSE25724. You can enter this number in the search box on the GEO homepage, or use [this direct link](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25724). In this study, the authors performed microarray analysis to evaluate differences in the transcriptome of type 2 diabetic human islets compared to non-diabetic islet samples. Human pancreatic islets were isolated from the pancreas of organ donors by collagenase digestion followed by density gradient purification, then hand-picked and cultured two days in M199 culture medium. Seven non-diabetic islet samples and six diabetic islet samples were hybridized to the Affymetrix GeneChip Human Genome U133A microarray. The results were published in: Dominguez V, et al. Class II phosphoinositide 3-kinase regulates exocytosis of insulin granules in pancreatic beta cells. *J Biol Chem 2011* Feb 11;286(6):4216-25. PMID: [21127054](http://www.ncbi.nlm.nih.gov/pubmed/21127054). 

You can go to the [GEO accession page](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25724) and download the raw .CEL file data directly (at the very bottom, under "Supplementary file"). However, we are going to use a Bioconductor package called [GEOquery](http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html) that allows us to *programmatically* access all data in GEO from within R. Take a look at the [GEOquery vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf) or Google around for "GEOquery tutorial" to learn more.

#### Load packages

First, we'll need to load the Bioconductor packages we'll be using:


```coffee
# Bioconductor packages.  Use the installation instructions at
# http://www.bioconductor.org/install/ if you haven't already installed
# these (install once, load every time)
library(Biobase)
library(GEOquery)
library(limma)
library(arrayQualityMetrics)
```


#### Acquire and pre-process the data

Next, let's get the data and verify that it's an [ExpressionSet](http://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf)-class object.


```coffee
gset <- getGEO("GSE25724", GSEMatrix = TRUE, destdir = ".", AnnotGPL = TRUE)
gset <- gset[[1]]
class(gset)
```


Next we'll fix the names of the variables describing the features (probes) on the array. The variable names are non-standard, having things like spaces (e.g. `"Gene title"`) and colons (e.g. `"GO:Function"`), which can cause problems downstream. The `make.names()` function makes syntactically valid names out of character vectors.


```coffee
# Removes spaces, colons, other weird symbols, inserts '.' instead.
fvarLabels(gset) <- make.names(fvarLabels(gset))
```


Next, let's take a look at the sample metadata. We can see that the first seven samples are non-diabetic islets and the last six are diabetic. Let's create a `condition` variable that indicates the disease state of each sample, and annotate the dataset's description with this new condition variable.


```coffee
# group names for all samples
View(pData(gset))
rep("ctrl", 7)
rep("diab", 6)
condition <- c(rep("ctrl", 7), rep("diab", 6))
condition <- as.factor(condition)
gset$description <- condition
```


#### Quality assessment 

Before we do any analysis, let's do some quality assessment using the [arrayQualityMetrics](http://www.bioconductor.org/packages/release/bioc/html/arrayQualityMetrics.html) Bioconductor package. Note that this can take a few minutes and a good deal of RAM for multiple samples. The `arrayQualityMetrics()` function produces a number of useful outputs in a tidy HTML report. Open up the index.html file produced in the output directory and look at the metadata overview, between-array comparisons (distance heatmap/dendrograms, PCA), array density distrubutions (boxplots and density plots), variance-mean dependence, and MA plots. The text below each figure in the QA report gives a description of each figure and what you should be seeing.


```coffee
arrayQualityMetrics(gset, outdir = "aqm-report", intgroup = "description")
```


#### Differential expression analysis

Finally, we'll do some data analysis. Below I'm using the [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html) (linear models for microarray) Bioconductor package for the analysis. Here we have a very simple two-group comparison, but limma can handle very complex multi-factorial designs including time-series, paired samples, and nested interaction models, among others. 

First we specify a *design matrix* which indicates the experimental condition of each sample. Then we specify a *contrast* matrix that indicates the exact comparisons we want to do between samples. For an experiment that's this simple we don't really need a contrast matrix since the contrast that we want to perform (diabetic vs non-diabetic) is implicit in the design, but we di it here anyway for clarity.


```coffee
# Set up and view the design matrix
design <- model.matrix(~description + 0, data = gset)
colnames(design) <- levels(condition)
design

# Set up and view the contrast matrix
cont.matrix <- makeContrasts(diab - ctrl, levels = design)
cont.matrix

# Fit a linear model, compute contrast coefficients and std errs, moderated
# t-tests and p-values:
fit <- lmFit(gset, design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit, proportion = 0.01)

# Extract the top 250 most differentially expressed genes, with an FDR
# correction
tt <- topTable(fit, adjust = "fdr", sort.by = "p", number = 250)

# Look at the results in RStudio
View(tt)

# Export the results to file
write.table(tt, file = "top-250-genes.csv", row.names = F, sep = ",")
```


#### Record package and version info with `sessionInfo()`

The `sessionInfo()` prints version information about R and any attached packages. It's a good practice to always run this command at the end of your R session and record it for the sake of reproducibility in the future.


```coffee
sessionInfo()
```

```
## R version 3.0.2 (2013-09-25)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  methods   stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
## [1] arrayQualityMetrics_3.18.0 limma_3.18.10             
## [3] GEOquery_2.28.0            Biobase_2.22.0            
## [5] BiocGenerics_0.8.0         BiocInstaller_1.12.0      
## 
## loaded via a namespace (and not attached):
##  [1] affy_1.40.0           affyio_1.30.0         affyPLM_1.38.0       
##  [4] annotate_1.40.0       AnnotationDbi_1.24.0  beadarray_2.12.0     
##  [7] BeadDataPackR_1.14.0  Biostrings_2.30.1     Cairo_1.5-5          
## [10] cluster_1.14.4        colorspace_1.2-4      DBI_0.2-7            
## [13] evaluate_0.5.1        formatR_0.10          Formula_1.1-1        
## [16] gcrma_2.34.0          genefilter_1.44.0     grid_3.0.2           
## [19] Hmisc_3.14-1          hwriter_1.3           IRanges_1.20.6       
## [22] knitr_1.5             lattice_0.20-24       latticeExtra_0.6-26  
## [25] plyr_1.8              preprocessCore_1.24.0 RColorBrewer_1.0-5   
## [28] RCurl_1.95-4.1        reshape2_1.2.2        RSQLite_0.11.4       
## [31] setRNG_2011.11-2      splines_3.0.2         stats4_3.0.2         
## [34] stringr_0.6.2         survival_2.37-7       SVGAnnotation_0.93-1 
## [37] tools_3.0.2           vsn_3.30.0            XML_3.95-0.2         
## [40] xtable_1.7-1          XVector_0.2.0         zlibbioc_1.8.0
```


