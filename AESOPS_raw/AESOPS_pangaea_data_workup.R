library(reshape2)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(lubridate)

############################################################################

#create list of .tab files in folder to loop over
fnames <- list.files(getwd(), pattern = "*.tab")

############################################################################

for (i in fnames) {

    site <- substr(i, 1, 13)
    
    #import raw tab file, as downloaded from PANGAEA
    # raw <- read.csv(fname, header=F)
    raw <- read.csv(i, header=F)
    
    #find row where column headers begin, whatever is after */ in this case
    count_start<-which(raw[,1]=='*/')+1
    #find row that data ends
    file_end <- nrow(raw)
    
    #extract rows and columns that specifically correspond to data, not anything else
    raw_counts <- as.character(raw[count_start:file_end,1])
    
    #split tab-delimited data into columns, coerce into dataframe, then transpose to get it right
    counts <- t(as.data.frame(strsplit(raw_counts, "\t")))
    #remove silly, long row names
    rownames(counts) <- NULL
    
    #convert to df
    counts_df <- as.data.frame(counts[2:nrow(counts),])
    names(counts_df) <- counts[1,]
    names(counts_df)[1:5] <- c("depth", "datetime", "datetime_end", "duration", "valve_flux")
    
    ###adjust where "G. species" shorthand is identical across genera
    #adjust: first "C. criophilum [#]" "T. antarctica [#]" to "Ch. criophilum [#]" "Tm. antarctica [#]"
    names(counts_df)[which(names(counts_df)=="C. criophilum [#]")[1]] <- "Ch. criophilum [#]"
    names(counts_df)[which(names(counts_df)=="T. antarctica [#]")[1]] <- "Tr. antarctica [#]"
    
    #find range of rows that correspond to taxonomic names used in dataset
    genus_names_start <- which(raw[,1]=='*/')-(ncol(counts) +1)
    genus_names_end <- which(raw[,1]=='License:\tCreative Commons Attribution 3.0 Unported (CC-BY-3.0)')-1
    
    #extract names vector based on above range
    genus_names_full <- as.character(raw[genus_names_start:genus_names_end,1])
    
    #truncate to genus only, remove initial \t
    genus_names_cut <- word(genus_names_full, 1)
    genus_names_cut <- substring(genus_names_cut, 2)
    
    #check that length of genus_names_cut == number of taxa in dataset, -1 because of depth
    length(genus_names_cut) == ncol(counts_df) - 1
    
    #make table that relates header names ("T. tumida"), to truncated genus name ("Thalassiosira")
    names_df <- data.frame("genus" = genus_names_cut, "full_unwieldy" = names(counts_df)[2:ncol(counts)])
    names_df$full_unwieldy <- as.character(names_df$full_unwieldy)
    
    #gather counts_df
    counts_df_long <- gather(counts_df, key=full_unwieldy, value=count, -depth, -datetime, -datetime_end, -duration, -valve_flux) 
    # counts_df_long <- gather(counts_df_long, key=full_unwieldy, value=count, full_unwieldy, "C. criophilum [#]", "T. antarctica [#]")
    counts_df_long$full_unwieldy <- as.character(counts_df_long$full_unwieldy)
    names_df$full_unwieldy <- as.character(names_df$full_unwieldy)
    counts_df_long <- left_join(counts_df_long, names_df, by="full_unwieldy")
    
    #aggregate genus by depth
    counts_df_long$count <- as.numeric(counts_df_long$count)
    counts_df_long <- aggregate(count ~ depth + datetime + datetime_end + duration + valve_flux + genus, data=counts_df_long, FUN="sum")
    
    #add site name
    counts_df_long$site <- site
    
    #write.csv()
    write.csv(counts_df_long, paste(site, "by_genus.csv", sep="_"), row.names = FALSE)
}

#combine created csv files
myMergedData <- 
  do.call(rbind,
          lapply(list.files(getwd(), pattern = "*.csv"), read.csv))

#export merged
write.csv(myMergedData, "merged_data_long.csv", row.names = FALSE)


#spread merged data
myMergedData_wide <- spread(myMergedData, key=genus, value=count)

