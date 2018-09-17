#to extract from an array of values the N, mean, SD, SE, ymin, ymax, and distance between for the T-test and graphs
args = commandArgs(TRUE)

if(length(args)!=3){
        cat("No arguments supplied.","\n")
    ##supply default values
        cat("Usage: Rscript extract_sample_info_tidyverse.R [1st FILE] [Type1] [Type2]","\n")
        cat("Ex: Rscript extract_sample_info_tidyverse.R t_test1.txt Female Male","\n")
        q(save="no")
}else{
	infile1 = args[1]
	infile2 = args[2]
	infile3 = args[3]
}

	options(warn=0)

	data1 <- read.table(infile1,header=T)
	ifelse(colnames(data1)[1]=="Gender",type1<-subset(data1,Gender==infile2),type1<-subset(data1,Type==infile2))
	ifelse(colnames(data1)[1]=="Gender",type2<-subset(data1,Gender==infile3),type2<-subset(data1,Type==infile3))
	mean_type1 <- mean(type1[,2])
	mean_type2 <- mean(type2[,2])
	cat("mean\t",mean_type1,"\t",mean_type2)

	
	