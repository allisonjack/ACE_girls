#to extract from an array of values the N, mean, SD, SE, ymin, ymax, and distance between for the T-test and graphs
args = commandArgs(TRUE)

if(length(args)!=2){
        cat("No arguments supplied.","\n")
    ##supply default values
        cat("Usage: Rscript extract_sample_values_info.R [1st FILE] [2nd FILE]","\n")
        cat("Ex: Rscript extract_sample_values_info.R t_test1.txt outfile_info.txt","\n")
        q(save="no")
}else{
	infile1 = args[1]
	infile2 = args[2]
}

	options(warn=0)
	library(stringi)
	library(plyr)

	data1 <- read.table(infile1,header=T)
	ifelse(colnames(data1)[1]=="Gender",tsize<-ddply(data1,c("Gender"),summarize,N=length(Gender),m=mean(Size),sds=sd(Size),se=sds/sqrt(N)),tsize<-ddply(data1,c("Type"),summarize,N=length(Type),m=mean(Size),sds=sd(Size),se=sds/sqrt(N)))
	tsize$ymin=tsize$m-tsize$se
	tsize$ymax=tsize$m+tsize$se
	ifelse(tsize$m[1]>tsize$m[2],(tsize$distance=tsize$ymin[1]-tsize$ymax[2]),(tsize$distance=tsize$ymin[2]-tsize$ymax[1]))
	tsize$mean_difference=tsize$m[1]-tsize$m[2]
	#print (tsize,quote=FALSE,sep="\t")
	#toprint[,1]=tsize$Type
	#toprint$mean=tsize$m
	#toprint
	
	output <- tsize[1:2,1:9]
	print (output,quote=FALSE,sep="\t",dec=".",append=FALSE,row.names=FALSE)
	write.table(tsize,file=infile2,append=TRUE,quote=FALSE,sep="\t",dec=".",col.names=TRUE,row.names=FALSE)
	
	