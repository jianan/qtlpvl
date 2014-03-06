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

Rd:
	R -e 'roxygen2::roxygenise(); devtools::document()'
# build package documentation and vignettes


vig:
	cd vignettes; R -e 'library(devtools);install("../");library(knitr);knit2html("qtlpvl.Rmd")'
# generate vignette html

clean:
	rm -rf $(TARGET_RCPP_ATTRIBUTE) man/*.Rd


