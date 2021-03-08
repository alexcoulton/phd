source("rotation1scripts_v4/scripts/r/functions.R")

sacha.old = read_csv("rotation1scripts_v4/original_data/sacha.introgression.data/OneDrive_2018-10-18/Introgression.detection/sacha.barley.order.cons.csv")

sacha.new = convert.aas.to.rqtl("rotation1scripts_v4/original_data/sacha.introgression.data/Introgression analysis input files/Analysis_2_genotype_calls.txt", call.code.numeric = T)

sacha.samp.type = read_delim("rotation1scripts_v4/original_data/sacha.introgression.data/Introgression analysis input files/Analysis_2_sample type.txt", delim = "\t")

#correct column names
colnames(sacha.samp.type) = c("sample", "sample_type")


sacha.samp.type = sacha.samp.type[c(which(sacha.samp.type$sample_type == "relative"), which(sacha.samp.type$sample_type == "wheat")), ]

sacha.samp.type = sacha.samp.type[c(which(sacha.samp.type$sample_type == "relative"), which(sacha.samp.type$sample_type == "wheat")), ]

r1 = min(which(sacha.samp.type$sample_type == "relative"))
r2 = min(which(sacha.samp.type$sample_type == "wheat"))

data.frame(c("relative start coord", "relative end coord", "wheat start coord", "wheat end coord"), c(r1, r2-1, r2, nrow(sacha.samp.type)))

sacha.new2 = as.data.frame(t(sacha.new))
sacha.new3 = cbind(rownames(sacha.new2), sacha.new2)
sacha.new3 = convert.to.character.data.frame(sacha.new3)

sacha.new4 = sacha.new3[, c(1, na.omit(match(sacha.samp.type$sample, sacha.new3[2, ])))]
sacha.new4 = sacha.new4[-1, ]


sacha.new5 = sacha.new4[c(1, match(sacha.old$probe_ID, sacha.new4[, 1])), ]
colnames(sacha.new5) = sacha.new5[1, ]
sacha.new5 = sacha.new5[-1, ]
sacha.new6 = cbind(sacha.old[, 1:3], sacha.new5)
sacha.new6 = sacha.new6[, -4]

#add a new row
sacha.new6 = rbind(sacha.new6[1, ], sacha.new6)
#clear new row
sacha.new6[1, ] = ""
#add sample type information to dataframe
sacha.new6[1, 4:ncol(sacha.new6)] = sacha.samp.type$sample_type

# write.csv(sacha.new6, "rotation1scripts_v4/original_data/sacha.introgression.data/new15112018/sacha.data.barley.cons.csv", row.names = F)

write.csv(sacha.new6, "rotation1scripts_v4/original_data/sacha.introgression.data/Introgression analysis input files/data_analysis2.csv", row.names = F)


