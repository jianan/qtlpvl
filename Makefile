
# build package documentation and vignettes
doc:
	R -e 'library(devtools);document(roclets=c("namespace", "rd"))'

vig:
	cd vignettes; R -e 'library(devtools);install("../");library(knitr);knit2html("qtlpvl.Rmd")'
