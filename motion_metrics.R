#---Load files---#
#Here define your file paths. These will be different for everyone depending
#on how you are organizing things on your local machine. I'm specifying 
#that my local files associated with the study are kept on an external hard-drive
#called "HIPPOCAMPUS" in a directory called "ACE_girls"...then I am specifying
#that I keep files with individual differences information in a sub-directory.
#So, change the following two lines to be appropriate for your setup.

STUDYDIR <- "/Volumes/HIPPOCAMPUS/ACE_girls"
DIFFDIR <- "/ANALYSIS/INDIV_DIFF"

#Read in the rms.txt file. Respect the header. Fill in empty cells with "NA."
rms <- read.table(sprintf("%s/%s/rms.txt", STUDYDIR, DIFFDIR), 
                  header     = TRUE,
                  fill       = TRUE,
                  na.strings = "")

#---Load packages---#
#Make sure you have installed the package first
library(plyr)

#---Summarize RMS data by participant
rms.abs.summ <- ddply(rms, .(ID), summarise,
                      rms.abs.M   = mean(ABS, na.rm = TRUE),
                      rms.abs.SD  = sd(ABS, na.rm = TRUE),
                      rms.abs.max = max(ABS, na.rm = TRUE))
rms.rel.summ <- ddply(rms, .(ID), summarise,
                      rms.rel.M   = mean(REL, na.rm = TRUE),
                      rms.rel.SD  = sd(REL, na.rm = TRUE),
                      rms.rel.max = max(REL, na.rm = TRUE),
                      rms.rel.min = min(REL, na.rm = TRUE))
#Create a new dataframe with the summary info
rms.df <- merge(rms.abs.summ, rms.rel.summ, by = c("ID"))

#Write this dataframe to a .csv file
write.csv(rms.df, sprintf("%s/%s/rms_df.csv", STUDYDIR, DIFFDIR), row.names=FALSE)

#Read in a motion outliers count file.
mo.df <- read.table(sprintf("%s/%s/motion_outliers_count.txt", STUDYDIR, DIFFDIR),
                    header = TRUE)

df <- merge(rms.df, mo.df, by = "ID")

#---Specify inclusion cutoffs---#
#These can be whatever you think is appropriate. Here I specify that 
#a participant cannot have more than 20% of total volumes be identified
#as motion outliers, and that the RMS motion cutoff is 4 mm. 

Total.vols <- 172
MO.cutoff <- .2*Total.vols
RMS.cutoff <- 4

#---Check RMS values---#

df$RMS.abs.include[rms.df$rms.abs.max>=RMS.cutoff] <- "no"
df$RMS.abs.include[rms.df$rms.abs.max<RMS.cutoff] <- "yes"
df$RMS.rel.include[rms.df$rms.rel.max>=RMS.cutoff] <- "no"
df$RMS.rel.include[rms.df$rms.rel.max<RMS.cutoff] <- "yes"

#---Check number of outliers (according to DVARS)---#

df$MO.include[df$OutlierCount>=MO.cutoff] <- "no"
df$MO.include[df$OutlierCount<MO.cutoff] <- "yes"

#---Decision rule---#

df$Include[df$RMS.abs.include =="no" | df$RMS.rel.include == "no" | df$MO.include == "no"] <- "no"
df$Include[df$RMS.abs.include == "yes" & df$RMS.rel.include == "yes" & df$MO.include == "yes"] <- "yes"

#Create a .csv file with all the information in the current dataframe.
write.csv(df, sprintf("%s/%s/motion.csv", STUDYDIR, DIFFDIR), row.names=FALSE)

#Create a .txt file with just the ID numbers of participants who "passed inspection"
df.include <- subset(df, df$Include == "yes", select = ID)
write.table(df.include, 
            sprintf("%s/%s/IDs_include.txt", STUDYDIR, DIFFDIR),
            row.names = FALSE,
            quote     = FALSE,
            col.names = FALSE)
            
