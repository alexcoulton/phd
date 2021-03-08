install.packages("gmailr")
suppressPackageStartupMessages(library(gmailr))
test_email <- mime(
    To = "alexcoulton@gmail.com",
    From = "apodforce@gmail.com",
    Subject = "this is just a gmailr test",
    body = "Can you hear me now?")
send_message(test_email)