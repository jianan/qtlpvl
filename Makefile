# build package documentation and vignettes
doc:
	R -e 'library(devtools);document(roclets=c("namespace", "rd"))'
	cd vignettes; R -e 'library(knitr);knit2html("qtlpvl.Rmd")'
