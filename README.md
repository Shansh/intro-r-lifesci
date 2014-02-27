# Introduction to R and Bioconductor for Cell Biology 8401

## R basics

### Basic operations

R as a glorified calculator. Try typing this in directly into the console. Then start typing this into the editor, and save your script. Use the run button, or press `CMD`+`Enter` (`Ctrl`+`Enter` on Windows).


```coffee
2 + 2
5 * 4
2^3
```


R Knows order of operations.


```coffee
2 + 3 * 4
```


### Functions

R has built-in functions.


```coffee
# Notice that this is a comment.  Anything behind a # is 'commented out' and
# is not run.
sqrt(144)
log(1000)
```


Get help by typing a question mark in front of the function's name:

```
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
`?`(seq)
seq(from = 2, to = 200, by = 4)
```


### Variables

You can store values in a variable. Use the `<-` operator for assignment. `=` also will work, but `<-` is conventional and preferred. Named variables are objects containing whatever is assigned to them. 


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


### Subsetting

Lorem ipsum.




