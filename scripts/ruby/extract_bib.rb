#!/usr/bin/env ruby

bibfile = "~/phd/rotation1scripts_v4/rmarkdown/full_zotero_library.bib"

matches = ARGF.read.scan(/@(.*?)[\.,:;\] ]/)
reg = "\\(" + matches.join("\\)\\|\\(") + "\\)"

system 'bibtool -X \'' + reg + '\' ' + bibfile
