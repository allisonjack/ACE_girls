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
