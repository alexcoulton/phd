#!/home/ac14037/bin/Rscript
bibfile = readLines('full_zotero_library.bib')
bibfile2 = bibfile[grep('@', bibfile)]
bf3 = strsplit(bibfile2, '\\{')
bf4 = lapply(bf3, function(x) tryCatch(x[[2]], error = function(e) NA))
bf5 = unlist(lapply(bf4, function(x) paste0('@', x)))
bf6 = unlist(lapply(bf5, function(x){
		x = strsplit(x, split = '')
		paste0(x[[1]][1:(length(x[[1]]) - 1)], collapse = "")

}))

writeLines(bf6, 'full_zotero_library_tags.bib')
