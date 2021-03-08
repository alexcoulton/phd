options(encoding = 'utf-8')
options(encoding = 'UTF-8')
options(encoding = '')
getOption('encoding')



render('~/project.phd.main/rotation1scripts_v4/rmarkdown/axc.analysis.markdown.Rmd', output_format = 'html_document')
render('~/project.phd.main/rotation1scripts_v4/rmarkdown/axc.analysis.markdown.Rmd', output_format = 'pdf_document')

#clear cache if needed
system('rm -rf ~/project.phd.main/rotation1scripts_v4/rmarkdown/axp.v2_cache')


render('~/project.phd.main/rotation1scripts_v4/rmarkdown/axp.v2.Rmd', output_format = 'html_document')
system('cpmarkdown.sh')
render('~/project.phd.main/rotation1scripts_v4/rmarkdown/axp.v2.Rmd', output_format = 'pdf_document')


render('~/project.phd.main/rotation1scripts_v4/rmarkdown/recombination_distribution.Rmd', output_format = 'html_document')
render('~/project.phd.main/rotation1scripts_v4/rmarkdown/recombination_distribution.Rmd', output_format = 'pdf_document')

render('~/project.phd.main/rotation1scripts_v4/rmarkdown/03-seg.dist.chapter.Rmd', output_format = 'html_document')



render('~/project.phd.main/rotation1scripts_v4/rmarkdown/01-introduction.Rmd', output_format = 'pdf_document')
render('~/project.phd.main/rotation1scripts_v4/rmarkdown/02-exome.Rmd', output_format = 'pdf_document')
render('~/project.phd.main/rotation1scripts_v4/rmarkdown/03-axp.v2.Rmd', output_format = 'pdf_document')
render('~/project.phd.main/rotation1scripts_v4/rmarkdown/04-seg.dist.chapter.Rmd', output_format = 'pdf_document')
render('~/project.phd.main/rotation1scripts_v4/rmarkdown/05-autocloner.Rmd', output_format = 'pdf_document')




library(ggplot2)
library(dplyr)
library(bookdown)
library(grid)
library(rmarkdown)
library(gridExtra)
library(pander)
library(ape)
library(kableExtra)
source('~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R')

#load('~/project.phd.main/rotation1scripts_v4/scripts/r/recombination/workspace22052020.RData')
#load('~/project.phd.main/24022020.FIPS.workspace.RData')
#load('~/project.phd.main/recombination.workspace.updated.thesis.RData') #make thesis and paper recomb data consistent

load('~/project.phd.main/thesis.workspace.25022021.RData')

load('~/project.phd.main/rotation1scripts_v4/saved.objects/main_plots')
load('~/project.phd.main/rotation1scripts_v4/saved.objects/supp_plots')
load('~/project.phd.main/rotation1scripts_v4/saved.objects/main_plots_exome_comp')


setwd('~/project.phd.main/rotation1scripts_v4/rmarkdown/')
file.remove('~/project.phd.main/rotation1scripts_v4/rmarkdown/_main.Rmd')



render_book('index.Rmd')   #note - to render to latex pdf, remove other output formats from index.Rmd and don't specify output_format in render_book()

#also note that i've temporarily disabled reference rendering in each individual Rmd file as it slows rendering of the overall pdf. This should be turned on at some point
							#otherwise latex_engine argument for xelatex is not recognised, and an error is generated on
#unicode characters

render_book('index.Rmd', output_format = 'word_document2')

