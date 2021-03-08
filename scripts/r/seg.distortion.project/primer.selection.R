#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####
conv.mult.align.dnastringset = function(mult.align){
    #converts a multiplealignment object to dnastringset
    mult.mat = as.matrix(mult.align)
    
    mult.mat3 = apply(mult.mat, 1, function(x) paste(x, collapse = ""))
    
    #convert each object in the list to a DNAStringSet object
    mult.mat4 = lapply(mult.mat3, DNAStringSet)
    
    #make new DNAStringSet object
    mult.mat5 = DNAStringSet()
    
    #add our sequences to the new DNAStringSet object one by one
    for(i in 1:length(mult.mat4)){
        mult.mat5 = c(mult.mat5, mult.mat4[[i]])
    }
    
    names(mult.mat5) = rownames(mult.align)
    
    mult.mat5
}

trim.multiple.alignment = function(input.path, main.seq.id, start.cut, end.cut){
    #takes a multiple alignment fasta file and trims it so that all sequences are the same length.
    #output is written to the same directory with .trimmed.fa appeneded to the original file name
    #
    #args:
    # input.path - character string of the path to the multiple alignment fasta file
    # main.seq.id - integer; number of the reference sequence to make cuts from in the fasta file
    # start.cut - integer; number of bases to cut upstream of the start of the reference 
    # end.cut - integer; number of bases to cut downstream of the end of the reference
    
    mult.align = readDNAMultipleAlignment(input.path)
    
    #convert the multiple alignment object to a matrix for easier
    #processing
    mult.mat = as.matrix(mult.align)
    
    #get the start coordinate of the chinese spring gene sequence
    #in the multiple alignment
    start.main.seq = min(which(mult.mat[main.seq.id, ] == "A"))
    #get the end coordinate
    end.main.seq = max(which(mult.mat[main.seq.id, ] == "A" | mult.mat[main.seq.id, ] == "T" | mult.mat[main.seq.id, ] == "G" | mult.mat[main.seq.id, ] == "C"))
    
    #subset the matrix
    mult.mat2 = mult.mat[, (start.main.seq - start.cut):(end.main.seq + end.cut)]
    
    #now begin conversion back to a multiple alignment object
    #first grab the sequences from the matrix and put them into a list
    mult.mat3 = apply(mult.mat2, 1, function(x) paste(x, collapse = ""))
    
    #convert each object in the list to a DNAStringSet object
    mult.mat4 = lapply(mult.mat3, DNAStringSet)
    
    #make new DNAStringSet object
    mult.mat5 = DNAStringSet()
    
    #add our sequences to the new DNAStringSet object one by one
    for(i in 1:length(mult.mat4)){
        mult.mat5 = c(mult.mat5, mult.mat4[[i]])
    }
    
    names(mult.mat5) = rownames(mult.align)
    
    #parse input path to make an output path - the original file with .trimmed.fa appended to the end
    output.path = strsplit(input.path, "\\/")
    output.path = output.path[[1]]
    file.name = output.path[length(output.path)]
    file.name = paste(file.name, ".trimmed.fa", sep = "")
    output.path = p(paste(output.path[1:(length(output.path)-1)], collapse = "/"), "/", file.name)
    
    
    writeXStringSet(mult.mat5, output.path)
    
}

remove.inserts = function(dna.string){
    #removes inserted bases ("-") from a DNAString object
    raw.seq = as.character(dna.string)
    
    #following code removes "-" from raw.seq
    raw.seq2 = strsplit(raw.seq, "")
    raw.seq2 = raw.seq2[[1]]
    raw.seq3 = raw.seq2[!raw.seq2 == "-"]
    raw.seq4 = paste(raw.seq3, collapse = "")
    
    DNAString(raw.seq4)
}

calculate.start.and.end.ranges = function(sequence.to.cut, product.size){
    #calculate coordinates for cutting of sequence into product sizes of ~ 700 bp
    #args:
    # sequence.to.cut - a DNAString object
    # product.size - integer specifying the desired length of the PCR product in bases
    
    if(missing(product.size)) product.size = 700
    
    if(product.size > length(sequence.to.cut)){
        start.range = 1
        end.range = length(sequence.to.cut)
    } else {
        num.divisions = ceiling(length(sequence.to.cut) / product.size)
        iterations = ceiling(length(sequence.to.cut) / num.divisions)
        
        start.range = seq(1, length(sequence.to.cut) + 100, iterations)
        end.range = seq(1, length(sequence.to.cut) + 100, iterations)
        
        start.range = start.range[-length(start.range)]
        end.range = end.range[-1] - 1
        
        end.range[length(end.range)] = length(sequence.to.cut)
        
        start.range[2:length(start.range)] = start.range[2:length(start.range)] - 100
        end.range[1:(length(end.range)-1)] = end.range[1:(length(end.range)-1)] + 100
    }
    
    ranges1 = list(start.range, end.range)
    names(ranges1) = c("start.range", "end.range")
    ranges1    
}

examine.primer.homologue.concordance = function(left.primer, right.primer, template.sequence){
    #given two primer sequences, checks how well the primers complement a given template sequence. 
    #returns the percentage sequence match against the template for each primer
    #args:
    # left.primer - character string of left primer sequence
    # right.primer - character string of right primer sequence
    left.align = pairwiseAlignment(DNAString(left.primer), template.sequence)
    right.align = pairwiseAlignment(reverseComplement(DNAString(right.primer)), template.sequence)
    
    alignment.pids = list(pid(left.align), pid(right.align))
    names(alignment.pids) = c("left.align.pid", "right.align.pid")
    alignment.pids
}

grab.primer.seq = function(primer.id, left.or.right, primer3.output){
    #extracts primer sequence from primer3 output file
    #args:
    # primer.id - integer (starting at 0)
    # left.or.right - character string; either "LEFT" or "RIGHT"
    primer1 = primer3.output[grep(p("PRIMER_", left.or.right, "_", primer.id, "_SEQUENCE"), primer3.output)]
    primer1 = strsplit(primer1, "=")
    primer1 = primer1[[1]][2]
    primer1
}


#     ____________________________________________________________________________
#     CHECK IF HOMOLOGUES ARE DIFF LENGTHS                                                                        ####


test1 = readDNAMultipleAlignment("rotation1scripts_v4/auto.primer.picker/gene4/seq/extended/alignments/alignment.w.primers.set11primer3.130-4834.txt.output.txt.fa")

test2 = conv.mult.align.dnastringset(test1)

test3 = lapply(test2, function(x) x[153:6231])

test4 = lapply(test3, remove.inserts)

test1 = readDNAMultipleAlignment("rotation1scripts_v4/auto.primer.picker/Traes700.full/seq/extended/alignments/all.align.rev.fa.trimmed.fa")

test2 = conv.mult.align.dnastringset(test1)

test3 = lapply(test2, function(x) x[1192:5992])

test4 = lapply(test3, remove.inserts)

lapply(test4, length)

#     ____________________________________________________________________________
#     PRIMER3 MULTIPLE ALIGNMENT                                                                                            ####

trim.multiple.alignment("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/seq/all.align.fa",
                                                7, 600, 600)

mult.align1 = readDNAMultipleAlignment("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/extended/all.align.rev.fa")


#     ____________________________________________________________________________
#     DETECT SNPS                                                                                                                         ####


count1 = make.counter()

grab.homeologous.snps = function(variety.rows, homologue.rows, multiple.alignment, allow.insertion){
    #takes a DNAMultipleAlignment class and returns the columns containing homeologous SNPs
    #args:
    # variety.rows - numeric vector containing the row coordinates of the varietal sequences from the same locus
    # homologue.rows - numeric vector containing the row coordinates of the homologous sequences (either paralogous or homeologous)
    # multiple.alignment - DNAMultipleAlignment class 
    # allow.insertion - boolean, indicates whether or not insertions "-" into the homologous sequences are valid SNPs
    
    mult.align.mat1 = convert.to.character.data.frame(as.data.frame(as.matrix(multiple.alignment)))    
    g = lapply(mult.align.mat1, function(x){
        if(length(unique(x[variety.rows])) == 1){ #are all of the varieties the same base at this locus?
            snp.check = unique(x[homologue.rows]) %in% unique(x[variety.rows]) #are any of the bases in any of the homologous sequences the same as the 
            #varietal sequences?
            if(T %in% snp.check){
                snp = 0
            } else {
                if(allow.insertion == T){
                    snp = 1    
                } else {
                    if("-" %in% x[homologue.rows] | "N" %in% x[homologue.rows] | "-" %in% x[variety.rows]){
                        snp = 0
                    } else {
                        snp = 1
                    }
                }
            }
        } else {
            snp = 0
        }
        snp    
    })
    
    unlist(g)
    
    
}

homologous.snps = grab.homeologous.snps(3:6, 1:2, mult.align1, F)

mult.align.mat1 = convert.to.character.data.frame(as.data.frame(as.matrix(mult.align1)))
mult.align.mat1 = add_row(mult.align.mat1)
mult.align.mat1[nrow(mult.align.mat1), ] = homologous.snps

mult.align.mat.start.base = min(grep("A|C|T|G", mult.align.mat1[7, ]))
mult.align.mat.end.base = max(grep("A|C|T|G", mult.align.mat1[7, ]))


write.csv(mult.align.mat1, "rotation1scripts_v4/processed_data/TILLING/mult.align.mat1.csv")

mult.align2 = conv.mult.align.dnastringset(mult.align1)

#need to transform coordinates of bases after removing hyphens
seq.mat = as.data.frame(as.matrix(mult.align2[[3]]))
seq.mat$snps = homologous.snps
seq.mat$start.base = ""
seq.mat$start.base[mult.align.mat.start.base] = 1
seq.mat$end.base = ""
seq.mat$end.base[mult.align.mat.end.base] = 1
seq.mat2 = seq.mat[-which(seq.mat[, 1] == "-"), ]
seq.mat2$new.coords = 1:nrow(seq.mat2)
snp.coords.after.filter = which(seq.mat2$snps == 1) #these are the coordinates of the homologous snps in this particular sequence after removing hyphens 
start.coord.after.filter = which(seq.mat2$start.base == 1) #coordinate of the start codon ATG after removing hyphens
end.coord.after.filter = which(seq.mat2$end.base == 1)

#primer3 uses coordinates starting from 0. Need to update the SNP coordinates to reflect this
snp.coords.after.filter = snp.coords.after.filter - 1
start.coord.after.filter = start.coord.after.filter - 1
end.coord.after.filter = end.coord.after.filter - 1

#     ____________________________________________________________________________
#     BEGIN MAIN LOOP                                                                                                                 ####

maximum.snp.coord = start.coord.after.filter
while(best.primer.file.end.coord < end.coord.after.filter){
    
    f.primer.candidates = snp.coords.after.filter[which(snp.coords.after.filter < maximum.snp.coord & snp.coords.after.filter > (maximum.snp.coord - 200))]
    
    #get all valid (ending on a homologous SNP) combinations of forward and reverse primer coordinates
    primer.combinations = lapply(f.primer.candidates, function(x){
        r.primer.candidates = snp.coords.after.filter[which(snp.coords.after.filter > (x + (start.coord.after.filter + 100)) & snp.coords.after.filter < (x + 700))]
        
        primer.combinations = expand.grid(x, r.primer.candidates)
        
        primer.combinations
    })
    
    primer.combinations = combine.list.of.data.frames(primer.combinations)
    
    # browser()
    
    #select one of the 5A varietal sequences for our template to make primers against
    template.sequence = mult.align2[[3]]
    
    template.sequence2 = remove.inserts(template.sequence)
    
    generate.primer3.input.files = function(template.sequence2, p3.seqid, product.size.min, product.size.max, left.end.coord,
                                                                                    right.end.coord){
        #args:
        # template.sequence2 - a DNAString object without inserts ("-"s)
        # p3.seqid - character string indicating name of sequence (used both inside the primer3 input file and in the title of the primer3 input file)
        
        
            
        #primer3 variables:
        p3.template = as.character(template.sequence2)
        p3.product.size.range = "100-3000"
        
        #note here that line breaks "\n" have to be added in manually as 
        #writeLines automatically adds a line break to the end of every line,
        #whilst primer3_core will not accept a file in which the last line 
        #has a line break on it
        primer3.input = c(p("SEQUENCE_ID=", p3.seqid, "\n"), 
                                            p("SEQUENCE_TEMPLATE=", p3.template, "\n"),
                                            p("PRIMER_PRODUCT_SIZE_RANGE=", p3.product.size.range, "\n"),
                                            p("SEQUENCE_FORCE_LEFT_END=", left.end.coord, "\n"),
                                            p("SEQUENCE_FORCE_RIGHT_END=", right.end.coord, "\n"),
                                            "=")
        
        output.filepath = file(p("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/input/primer3.", p3.seqid, ".txt"), "wb")
        writeLines(primer3.input, output.filepath, sep = "")
        close(output.filepath)
        "Done"
    }
    
    #generate all primer input files using every valid combination of forward primer coordinates and reverse primer coordinates
    Map(function(f.primer1, r.primer1){
        generate.primer3.input.files(template.sequence2, p((f.primer1 + 1), "-", (r.primer1 + 1)), 100, 750, f.primer1, r.primer1)
    }, primer.combinations[, 1], primer.combinations[, 2])
    
    print("Sync with Wilkins and run run.primer3.sh, then sync again.")
    browser()


#     ____________________________________________________________________________
#     PARSE OUTPUT FILES                                                                                                            ####

    output.files = list.files("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/output/")
    
    #parse primer3 output files and return the best paired primer penalty score in each one
    output.penalties = unlist(lapply(output.files, function(x){
        current.file = paste("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/output/", x, sep = "")
        current.file.lines = readLines(current.file)
        if(length(grep("PRIMER_PAIR_0_PENALTY", current.file.lines)) > 0){
            pen1 = current.file.lines[grep("PRIMER_PAIR_0_PENALTY", current.file.lines)]
            pen1 = strsplit(pen1, "=")
            return(as.numeric(pen1[[1]][2]))
        } else {
            return(1000) #if no PRIMER_PAIR_0_PENALTY in output file, return arbitrarily large penalty
        }
    }))
    
    browser()
    

    
    best.primer.coord = which(output.penalties == min(output.penalties))
    output.penalties[best.primer.coord]
    output.penalties.sorted = sort(output.penalties)
    
    best.primer.file = output.files[best.primer.coord]
    best.primer.file1 = strsplit(best.primer.file, "\\.")
    best.primer.file1 = best.primer.file1[[1]][2]
    best.primer.file1 = strsplit(best.primer.file1, "-")
    best.primer.file.start.coord = as.numeric(best.primer.file1[[1]])[1]
    best.primer.file.end.coord = as.numeric(best.primer.file1[[1]])[2]
    
    maximum.snp.coord = best.primer.file.end.coord #The name of this variable originates from the first iteration. Of cause this is not truly the start.coord.after.filter on subsequent iterations
    
    list.of.best.primer.files = list()
    list.of.best.primer.files = c(list.of.best.primer.files, best.primer.file)
    
    file.copy(p("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/output/", best.primer.file), p("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/best.primers/", best.primer.file))
    
    #cleanup intermediate input and output files before next iteration
    input.files = list.files("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/input/")
    input.files = paste("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/input/", input.files, sep = "")
    lapply(input.files, function(x) file.remove(x))
    
    output.files.full.path = paste("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/output/", output.files, sep = "")
    
    lapply(output.files.full.path, function(x) file.remove(x))
    
    print("Found best file:")
    print(best.primer.file)
}

#     ____________________________________________________________________________
#     PARSE PRIMER3 OUTPUT                                                                                                        ####

get.number.of.primers = function(primer3.output.list){
    #grab the number of primers from the primer3 output file and make 
    #vector for each primer id
    #args:
    # primer3.output.list - list containing primer 3 output obtained using readLines()
    
    num.primers = strsplit(primer3.output.list[max(grep("PRIMER_PAIR_", primer3.output.list))], "_")
    num.primers = num.primers[[1]]
    num.primers = as.numeric(num.primers[3])
    num.primers = 0:num.primers
    num.primers
}

make.primer.concordance.df = function(dna.stringset1, primer3.output, return.avg){
    #takes a multiple alignment DNAStringSet object (produced with conv.mult.align.dnastringset) as well as primer3 output file
    #and tests all of the primers for percentage id against all sequences in the set
    #args:
    # dna.stringset1 - DNAStringSet object produced with conv.mult.align.dnastringset()
    # primer3.output - primer3 output read with readLines()
    # return.avg - boolean, if true, return average PID for left and right primers
    
    if(missing(return.avg)) return.avg = F
    
    #make vectors for column names - both left and right primers
    col.names1 = paste("p.l.1", names(dna.stringset1), sep = "")
    col.names2 = paste("p.r.1", names(dna.stringset1), sep = "")
    #interleave the previous two vectors
    col.names3 = c(rbind(col.names1, col.names2))
    
    primer.pid.df = newdf(col.names3, no.rows = T)
    primer.pid.df.avg = newdf(names(dna.stringset1), no.rows = T)
    col.vector = seq(1, ncol(primer.pid.df), 2)
    
    
    for(i in num.primers){
        count = 1
        primer.pid.df = add_row(primer.pid.df)
        
        for(i2 in 1:length(dna.stringset1)){
            left.primer1 = grab.primer.seq(i, "LEFT", primer3.output)
            
            right.primer1 = grab.primer.seq(i, "RIGHT", primer3.output)
            
            primer.conc = examine.primer.homologue.concordance(left.primer1, right.primer1, remove.inserts(dna.stringset1[[i2]]))
            
            #populate dataframe with pid values 
            primer.pid.df[(i + 1), col.vector[count]] = primer.conc[[1]]
            primer.pid.df[(i + 1), (col.vector[count] + 1)] = primer.conc[[2]]
            
            primer.pid.df.avg[(i + 1), i2] = (primer.conc[[1]] + primer.conc[[2]]) / 2
            
            
            
            count = count + 1
        }
    }
    
    if(return.avg == T){
        return(primer.pid.df.avg)
    } else {
        return(primer.pid.df)
    }
}

examine.primer.df = function(primer.df, homologue.cols){
    #examine a primer.df made with make.primer.concordance.df() and grab the best primers based on their PIDs
    #args:
    # primer.df - a dataframe made with make.primer.concordance.df (average argument = T)
    # homologue.cols - numeric vector indicating the columns containing homologues (paralogues / homeologues)
    
    max.homologue.pid = as.numeric()
    for(i in 1:nrow(primer.df)){
        max.homologue.pid = c(max.homologue.pid, max(as.numeric(unlist(lapply(homologue.cols, function(x) primer.df[i, x])))))
    }
    max.homologue.pid
}

primer.files = list.files("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/")
primer.files2 = primer.files[grep("primer3.input.cut", primer.files)]

#get the number of primer3 files
num.primer.files = max(as.numeric(substr(primer.files2, 19, 19)))

primer3.output = readLines("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/primer3.input.cut.5.txt.output.txt")

num.primers = get.number.of.primers(primer3.output)

primer1.df = make.primer.concordance.df(mult.align2, primer3.output, T)

best.primer.pids = examine.primer.df(primer1.df, 1:2)

b.p.pid.coords = unique(unlist(lapply(sort(best.primer.pids)[1:2], function(x){
    which(best.primer.pids == x)
})))







