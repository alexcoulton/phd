dirs1 = list.dirs(recursive = F)

dirs2 = paste0(dirs1, "/bams/a8.10x.cov")
dirs3 = paste0(dirs1, "/bams/a8.5x.cov")
dirs4 = paste0(dirs1, "/bams/a8.3x.cov")

file.sizes1 = sapply(dirs2, function(x) file.info(x)$size)
file.sizes2 = sum(sapply(dirs3, function(x) file.info(x)$size) / 1000000000)
file.sizes3 = sum(sapply(dirs4, function(x) file.info(x)$size) / 1000000000)


# file.sizes1 = sapply(all.coverage.files, function(x) file.info(x)$size)
# if(length(which(file.sizes1 == 0)) > 0){
# 	all.coverage.files2 = all.coverage.files[-which(file.sizes1 == 0)]
# } else {
# 	all.coverage.files2 = all.coverage.files
# }

make.shared.cov = function(dirs.list, coverage){
	coverage1 = read.delim(dirs.list[[1]], sep = "\t", stringsAsFactors = F, header = F)

	for(i in 2:length(dirs.list)){	
		g = read.delim(dirs.list[i], sep = "\t", stringsAsFactors = F, header = F)
		num.shared = length(which(coverage1$V2 %in% g$V2)) / nrow(coverage1)
		coverage1 = coverage1[which(coverage1$V2 %in% g$V2), ]
		
		print(paste0("done comparison with ", dirs.list[i]))
		print(paste0("% shared is ", num.shared))
	}

	write.csv(coverage1, paste0("../csv/all.coverage.", coverage, "x.csv"), row.names = F)
}

make.shared.cov(dirs3, 5) # 5x coverage

make.shared.cov(dirs4, 3) # 3x coverage


#### GENERATE NEW VCFS #### 

all.vcf.files = paste0(dirs1, "/all_w_indels2.vcf")


coverage1.5 = split(coverage1, coverage1$V1) 


library(dplyr)

#grab only SNP locations that have at least 20x coverage in all samples
for(i in all.vcf.files){	
	
	i2 = gsub("all_w_indels2.vcf", "all_w_indels2.filtered.coverage_10x.vcf", i)        
	g = read.delim(i, sep = "\t", comment.char = "#", stringsAsFactors = F, header = F)
	
	g2 = split(g, g$V1)

	g3 = Map(function(vcf1, cov1){
		vcf1[which(vcf1$V2 %in% cov1$V2), ]
		}, g2, coverage1.5)

	g4 = bind_rows(g3)
	
	write.table(g4, i2, sep = "\t", col.names = F, quote = F, row.names = F)
	print(paste0("done ", i))
}
