SRC_CPP=$(filter-out src/RcppExports.cpp, $(wildcard src/*.cpp))
TARGET_RCPP_ATTRIBUTE=R/RcppExports.R src/RcppExports.cpp

.PHONY:	prebuild
prebuild: $(TARGET_RCPP_ATTRIBUTE) Rd
	@

$(TARGET_RCPP_ATTRIBUTE) : $(SRC_CPP)
	R --slave -e "library(Rcpp); Rcpp::compileAttributes();" && \
	touch $(TARGET_RCPP_ATTRIBUTE)
## `compileAttributes`, called by `document` may not touch the exports
## files, if no update is needed. Use touch to update the modification
## time.

Rd: $(filter-out R/RcppExports.R, $(wildcard R/*.R))
	R -e 'roxygen2::roxygenise(); devtools::document()'
# build package documentation and vignettes

vig: vignettes/qtlpvl.html
vignettes/qtlpvl.html: vignettes/qtlpvl.Rmd Rd
	cd vignettes; R -e 'devtools::install("../"); knitr::knit2html("qtlpvl.Rmd", "qtlpvl.html")'
# generate vignette html

clean:
	rm -rf $(TARGET_RCPP_ATTRIBUTE) man/*.Rd NAMESPACE

cleanall:
	rm -rf $(TARGET_RCPP_ATTRIBUTE) man/*.Rd NAMESPACE \
           vignettes/figure vignettes/*.md vignettes/*.html 
