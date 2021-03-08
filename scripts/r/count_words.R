#!/home/ac14037/bin/Rscript
suppressMessages(library(dplyr))

args = commandArgs(trailingOnly = T)

options(warn = -1)
suppressMessages(library(academicWriteR))
#count_words = function(file){
    #wc <- readr::read_lines(file) %>% tibble::as_data_frame(.) %>%
        #remove_front_matter(.) %>% remove_code_chunks(.) %>%
        #remove_inline_code(.) %>% remove_html_comment(.) %>%
        #tidytext::unnest_tokens(., output = .data$words, input = value) %>%
        #nrow(.)
    #return(wc)
#}

paste0('Words: ', count_words(args[1]))
