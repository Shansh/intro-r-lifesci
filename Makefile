OUTS := $(patsubst %.Rmd,%.html,$(wildcard *.Rmd)) $(patsubst %.Rmd,%.md,$(wildcard *.Rmd))

all: $(OUTS)

clean:
	rm $(OUTS)

%.html: %.Rmd
	Rscript -e "knitr::knit2html('$*.Rmd')"
	Rscript -e "knitr::purl('$*.Rmd')"
	gsed -i 's/```r/```coffee/g' $*.md
