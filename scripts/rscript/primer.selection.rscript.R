#!/home/ac14037/bin/Rscript
setwd("/home/ac14037/project.phd.main/")
library(Biostrings)
library(dplyr)
library(tibble)
source("rotation1scripts_v4/scripts/r/functions.R")

system("/home/ac14037/project.phd.main/rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/clean.folders.sh")

#   ____________________________________________________________________________
#   DEFINE FUNCTIONS                                                        ####
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

#   ____________________________________________________________________________
#   PRIMER3 MULTIPLE ALIGNMENT                                              ####

trim.multiple.alignment("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/extended/alignments/all.align.rev.fa", 9, 2000, 2000)

mult.align1 = readDNAMultipleAlignment("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/extended/alignments/all.align.rev.fa.trimmed.fa")


#   ____________________________________________________________________________
#   DETECT SNPS                                                             ####

homologue.rows = 1:3
variety.rows = 4:8
main.gene.row = 9
template.row = 8


count1 = make.counter()

grab.homeologous.snps = function(variety.rows, homologue.rows, multiple.alignment, allow.insertion, allow.N){
  #takes a DNAMultipleAlignment object and returns a numeric vector of the column coordinates containing homeologous SNPs
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
        if(allow.insertion == T & allow.N == T){
          snp = 1  
        } else {
          if(allow.insertion == F & allow.N == T){
            if("-" %in% x[homologue.rows] | "-" %in% x[variety.rows]){
              snp = 0
            } else {
              snp = 1
            }
          }

          if(allow.insertion == T & allow.N == F){
            if("N" %in% x[homologue.rows]){
              snp = 0
            } else {
              snp = 1
            }
          }

          if(allow.insertion == F & allow.N == F){
            if("-" %in% x[homologue.rows] | "-" %in% x[variety.rows] | "N" %in% x[homologue.rows]){
              snp = 0
            } else {
              snp = 1
            }
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

homologous.snps = grab.homeologous.snps(variety.rows, homologue.rows, mult.align1, F, T)

# print("homologous snps:")
# print(which(homologous.snps == 1))

get.start.and.end.base = function(multiple.alignment1, main.sequence.row.number){
  #Returns the positions of the start base and the end base of a specified sequence in a multiple alignment
  #args:
  # multiple.alignment1 - a DNAMultipleAlignment object
  # main.sequence.row.number - Integer; row number of the desired sequence to return positions for
  multiple.alignment1 = convert.to.character.data.frame(as.data.frame(as.matrix(multiple.alignment1)))
  mult.align.mat.start.base = min(grep("A|C|T|G", multiple.alignment1[main.sequence.row.number, ]))
  mult.align.mat.end.base = max(grep("A|C|T|G", multiple.alignment1[main.sequence.row.number, ]))  
  g = c(mult.align.mat.start.base, mult.align.mat.end.base)
  names(g) = c("start.base", "end.base")
  g
}

main.start.end = get.start.and.end.base(mult.align1, main.gene.row)

get.coordinates.after.removing.hyphens = function(mult.align1, sequence.row.number, homologous.snps1, start.base1, end.base1){
  #Returns the coordinates of homologous SNPs, start base and end base for a specified sequence in a multiple alignment file after removing hyphens
  #args:
  #mult.align1 - a DNAMultipleAlignment file
  #sequence.row.number - Integer, the row coordinate of the desired sequence to transform coordinate for in the multiple alignment
  #homologous.snps1 - Numeric vector, obtained using grab.homoelogous.snps()
  #start.base1 - Integer, obtained using get.start.and.end.base()
  #end.base1 - Integer, obtained using get.start.and.end.base()
  mult.align2 = conv.mult.align.dnastringset(mult.align1)

  #need to transform coordinates of bases after removing hyphens
  seq.mat = as.data.frame(as.matrix(mult.align2[[sequence.row.number]]))

  seq.mat$snps = homologous.snps1
  seq.mat$start.base = ""
  seq.mat$start.base[start.base1] = 1
  seq.mat$end.base = ""
  seq.mat$end.base[end.base1] = 1

  seq.mat2 = seq.mat[-which(seq.mat[, 1] == "-"), ] #remove hyphens
  seq.mat2$new.coords = 1:nrow(seq.mat2)
  snp.coords.after.filter = which(seq.mat2$snps == 1) #these are the coordinates of the homologous snps in this particular sequence after removing hyphens 
  start.coord.after.filter = which(seq.mat2$start.base == 1) #coordinate of the start codon ATG after removing hyphens
  end.coord.after.filter = which(seq.mat2$end.base == 1)

  #primer3 uses coordinates starting from 0. Need to update the SNP coordinates to reflect this
  snp.coords.after.filter = snp.coords.after.filter - 1
  start.coord.after.filter = start.coord.after.filter - 1
  end.coord.after.filter = end.coord.after.filter - 1

  g = list(snp.coords.after.filter, start.coord.after.filter, end.coord.after.filter)
  names(g) = c("snp.coords", "start.coord", "end.coord")
  g
}

coords = get.coordinates.after.removing.hyphens(mult.align1, template.row, homologous.snps, main.start.end[1], main.start.end[2])

print(coords)

#   ____________________________________________________________________________
#   BEGIN MAIN LOOP                                                         ####

find.best.primers = function(multiple.alignment, template.sequence.row.number, snp.coords.after.filter, start.coord.after.filter, end.coord.after.filter, product.size.range, span.whole.gene){
  #Automatically obtains primer sequences
  #args:
  # multiple.alignment - a DNAMultipleAlignment object
  # template.sequence.row.number - Integer; the multiple alignment row of the sequence to use as a template in primer3
  # snp.coords.after.filter - Numeric vector; obtained using grab.homeologous.snps() and then get.coordinates.after.removing.hyphens()
  # start.coord.after.filter - Integer; position of the first base of the start codon after removing hyphens
  # end.coord.after.filter - Integer; position of the final base in the coding sequence after removing hyphens
  # product.size.range - a numeric vector with two elements, the first being the minimum product size, the second the maximum

  if(missing(span.whole.gene)) span.whole.gene = F

  list.best.primer.start.coords = as.numeric()
  list.best.primer.end.coords = as.numeric()
  mult.align2 = conv.mult.align.dnastringset(multiple.alignment)

  for(i in 1:2){
    # Initialize some variables outside of loop  
    best.primer.file.end.coord = 0
    minimum.snp.coord = 0
    maximum.snp.coord = start.coord.after.filter
    count = 1

    while(best.primer.file.end.coord < end.coord.after.filter){
    
      f.primer.candidates = snp.coords.after.filter[which(snp.coords.after.filter < (maximum.snp.coord - 1) & snp.coords.after.filter > minimum.snp.coord)]      
      f.matches = c(na.omit(match(list.best.primer.end.coords, f.primer.candidates)), na.omit(match(list.best.primer.start.coords, f.primer.candidates)))
      
      print("list.best.primer.start.coords")
      print(list.best.primer.start.coords)

      print("f.matches:")
      print(f.matches)

      if(length(f.matches) > 0) f.primer.candidates = f.primer.candidates[-f.matches]
      
      
      #get all valid (ending on a homologous SNP) combinations of forward and reverse primer coordinates
      primer.combinations = lapply(f.primer.candidates, function(x){

        if(span.whole.gene == T){
          r.primer.candidates = snp.coords.after.filter[which(snp.coords.after.filter > (x + product.size.range[1]) & snp.coords.after.filter > end.coord.after.filter & snp.coords.after.filter < (x + product.size.range[2]))]
        } else {
          r.primer.candidates = snp.coords.after.filter[which(snp.coords.after.filter > (x + product.size.range[1]) & snp.coords.after.filter > best.primer.file.end.coord & snp.coords.after.filter < (x + product.size.range[2]))]  
        }
        
        r.matches = c(na.omit(match(list.best.primer.start.coords, r.primer.candidates)), na.omit(match(list.best.primer.end.coords, r.primer.candidates)))
        if(length(r.matches) > 0) r.primer.candidates = r.primer.candidates[-r.matches]
        
        
        primer.combinations = expand.grid(x, r.primer.candidates)
        
        primer.combinations
      })
      
      primer.combinations = combine.list.of.data.frames(primer.combinations)
      
      # browser()
      
      #select one of the 5A varietal sequences for our template to make primers against
      template.sequence = mult.align2[[template.sequence.row.number]]
      
      template.sequence2 = remove.inserts(template.sequence)
      
      generate.primer3.input.files = function(template.sequence2, p3.seqid, product.size.min, product.size.max, left.end.coord,
                                              right.end.coord){
        #args:
        # template.sequence2 - a DNAString object without inserts ("-"s)
        # p3.seqid - character string indicating name of sequence (used both inside the primer3 input file and in the title of the primer3 input file)
        
        
          
        #primer3 variables:
        p3.template = as.character(template.sequence2)
        p3.product.size.range = "100-10000"
        
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
      
      # print("Sync with Wilkins and run run.primer3.sh, then sync again.")
      # browser()
      system("~/project.phd.main/rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/run.primer3.sh")

    #   ____________________________________________________________________________
    #   PARSE OUTPUT FILES                                                      ####

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
      
      # browser()
      

      
      best.primer.coord = which(output.penalties == min(output.penalties))
      output.penalties[best.primer.coord]
      output.penalties.sorted = sort(output.penalties)
      if(all(output.penalties.sorted == 1000)) print("NO VALID PRIMERS FOR THIS RANGE")
      print(output.penalties.sorted[1:10])
      
      best.primer.file = output.files[best.primer.coord]
      best.primer.file1 = strsplit(best.primer.file, "\\.")
      best.primer.file1 = best.primer.file1[[1]][2]
      best.primer.file1 = strsplit(best.primer.file1, "-")
      best.primer.file.start.coord = as.numeric(best.primer.file1[[1]])[1]
      best.primer.file.end.coord = as.numeric(best.primer.file1[[1]])[2]
      
      minimum.snp.coord = best.primer.file.start.coord
      maximum.snp.coord = best.primer.file.end.coord #The name of this variable originates from the first iteration. Of cause this is not truly the start.coord.after.filter on subsequent iterations
      
      list.best.primer.start.coords = c(list.best.primer.start.coords, (best.primer.file.start.coord - 1)) #subtract 1 as these coordinates are on a scale beginning with 1
      list.best.primer.end.coords = c(list.best.primer.end.coords, (best.primer.file.end.coord - 1))
      
      if(!dir.exists(p("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/best.primers/set", i))){
        dir.create(p("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/best.primers/set", i))
      }

      file.copy(p("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/output/", best.primer.file), p("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/best.primers/set", i, "/", count, best.primer.file))
      
      #cleanup intermediate input and output files before next iteration
      input.files = list.files("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/input/")
      input.files = paste("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/input/", input.files, sep = "")
      lapply(input.files, function(x) file.remove(x))
      
      output.files.full.path = paste("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/output/", output.files, sep = "")
      
      lapply(output.files.full.path, function(x) file.remove(x))
      
      print("Found best file:")
      print(best.primer.file)

      count = count + 1
    }
  }
}


find.best.primers(mult.align1, template.row, coords[[1]], coords[[2]], coords[[3]], c(3000, 30000), T)
