#spreadsheet for axiom samples
#need create 24x16 dataframe
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
b=data.frame(matrix(NA, nrow=16, ncol=24))

q1c=seq(1,24,2)
q1r=seq(1,16,2)

for(i in 1:length(q1c)){
    for(p in 1:length(q1r)){
        b[q1r[p],q1c[i]]=paste("Q1_10_", LETTERS[1:8][p], i, sep="")
    }
}



q2c=seq(2,24,2)
q2r=seq(1,16,2)



for(i in 1:length(q2c)){
    for(p in 1:length(q2r)){
        b[q2r[p],q2c[i]]=paste("Q2_26_", LETTERS[1:8][p], i, sep="")
    }
}

q3c=seq(1,24,2)
q3r=seq(2,16,2)

for(i in 1:length(q3c)){
    for(p in 1:length(q3r)){
        b[q3r[p],q3c[i]]=paste("Q3_14_", LETTERS[1:8][p], i, sep="")
    }
}

q4c=seq(2,24,2)
q4r=seq(2,16,2)

for(i in 1:length(q4c)){
    for(p in 1:length(q4r)){
        b[q4r[p],q4c[i]]=paste("Q4_28_", LETTERS[1:8][p], i, sep="")
    }
}

b[b=="Q1_10_B12"]="Paragon_2_samp1"
b[b=="Q1_10_A12"]="Apogee_2_samp1"

b[b=="Q4_28_A12"]="Apogee_1_samp1"
b[b=="Q4_28_H12"]="Paragon_1_samp1"

b[b=="Q3_14_A12"]="Apogee_2_samp2"
b[b=="Q3_14_C12"]="Paragon_2_samp2"

b[b=="Q2_26_A5"]="Apogee_1_samp2"
b[b=="Q2_26_B5"]="Apogee_2_samp3"
b[b=="Q2_26_C11"]="Paragon_1_samp_2"
b[b=="Q2_26_A12"]="Paragon_2_samp3"
b[b=="Q2_26_B9"]="Paragon_3_samp1"


rownames(b)=LETTERS[1:16]
rownames(b)
colnames(b)=1:24
colnames(b)
write.csv(b, "alex_axiom_384_spreadsheet.csv")
