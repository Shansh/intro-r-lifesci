OUTS = README.html README.md

all: $(OUTS)

clean:
	rm -f $(OUTS)

%.html: %.Rmd
	Rscript -e "knitr::knit2html('$*.Rmd')"
	gsed -i 's/```r/```coffee/g' README.md
