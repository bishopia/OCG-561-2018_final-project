#BIO-561 Final Poster
#Ian Bishop
#December 6th, 2018


library(here)


#libraries
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(dplyr)
library(scales)
library(gridExtra)


#import data
hasle_data <- read.csv("hasle_data/Hasle_Bretagg_counts.csv", header=TRUE, stringsAsFactors = FALSE)
grig_data <- read.csv("grig_data/datasets/workup/merged_data_long.csv", header=TRUE, stringsAsFactors = FALSE) 
grig_data$project <- "AESOPS"
assmy_data <- read.csv("assmy_data/datasets/raw_out_of_patch/workup/merged_data_long.csv", header=TRUE, stringsAsFactors = FALSE) 
assmy_data$project <- "EISENEX"

#lat long depth
metadata <- read.csv("site_project_metadata.csv", header=TRUE, stringsAsFactors = FALSE)

#latitude diverging gradient palette
mypal <- c("#7fbf7b", "#f7f7f7", "#af8dc3")
mypal2 <- c("#7fbf7b", "#1f78b4", "#af8dc3")
my.palette <- c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')

mypal3 <- c("#05668D", "#CA5310", "#840032", "#5A5353")



#########################
###-------HASLE-------###
#########################


###---initial data transformation---###

#collapse to genus
hasle_gen <- aggregate(cells_per_liter ~ project + site + genus, data=hasle_data, FUN=sum)

hasle_sums <- hasle_gen %>% group_by(project, site) %>% summarise(valve_sum = sum(cells_per_liter))
hasle_sums$site <- as.character(hasle_sums$site)
#add zone
hasle_sums <- left_join(hasle_sums, metadata[,c("site", "apf_side")], by="site")
hasle_sums$site <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q")

psum1 <- ggplot(data=hasle_sums, aes(x=site, y=valve_sum, fill=apf_side)) + 
  geom_bar(width=.5, stat="identity") +
  xlab("Station") +
  ylab("Diatom valve abundance") +
  scale_fill_manual(values=c(mypal3)) +
  theme_linedraw(30)


#create sum df
df1 <- hasle_gen %>%
  group_by(site) %>%
  summarize(sum = sum(cells_per_liter))
#merge sums to hasle_gen
df2 <- left_join(hasle_gen, df1, by="site")
#calculate rel_abundance and add as column
df2$rel_abund <- round(df2$cells_per_liter/df2$sum, 4)
#remove intermediate columns
hasle_gen <- select(df2, c(project, site, genus, rel_abund))

#remove all rel_abundances of 0.0000
hasle_gen <- hasle_gen %>% filter(rel_abund!=0)

#delete intermediate dfs
rm(df1)
rm(df2)

# df_0.001_rel
hasle_gen_0.001 <- hasle_gen %>% filter(rel_abund>0.001)

# #add lat long
hasle_gen$site <- as.character(hasle_gen$site)
hasle_gen <- left_join(hasle_gen, metadata[c("site", "lat", "long")], by="site")


###---calculate metrics, including species richness, diversity, evenness---###

#summarize genus count per site
richness_gen <- hasle_gen %>%
  group_by(site) %>%
  count()
#rename columns
names(richness_gen) <- c("site", "richness")
richness_gen$site <- as.character(richness_gen$site)

#create wide rel_abund df
hasle_gen_wide <- spread(hasle_gen, key=genus, value=rel_abund)
hasle_gen_wide[is.na(hasle_gen_wide)] <- 0
#calculate diversity and evenness
div_hasle <- diversity(as.matrix(hasle_gen_wide[,5:ncol(hasle_gen_wide)]), "shannon")
s_hasle <- specnumber(as.matrix(hasle_gen_wide[,5:ncol(hasle_gen_wide)])) ## rowSums(BCI > 0) does the same...
eve_hasle <- div_hasle/log(s_hasle)

#create df with all metrics and site and project; to rbind with other projects downstream
hasle_metrics <- data_frame(project = hasle_gen_wide$project, 
                            site = as.character(hasle_gen_wide$site),
                            datetime = NA,
                            depth = NA,
                            div = div_hasle, 
                            eve = eve_hasle,
                            rich = richness_gen$richness)
#add latitude from metadata
hasle_metrics <- left_join(hasle_metrics, metadata[,c("site", "lat")], by="site")


# #spread hasle_gen
# hasle_gen_wide <- spread(hasle_gen, key=genus, value=rel_abund)
# hasle_gen_wide[is.na(hasle_gen_wide)] <- 0


###---HASLE_NMDS---###

hasle_matrix <- as.matrix(hasle_gen_wide[,5:ncol(hasle_gen_wide)])
hasle_NMDS=metaMDS(hasle_matrix,k=2,trymax=100)
hasle.scores <- as.data.frame(scores(hasle_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
hasle.scores$site <- hasle_gen_wide$site
hasle.scores$lat <- hasle_gen_wide$lat
hasle.scores$long <- hasle_gen_wide$long


###---HASLE_PLOTS---###

p1 <- ggplot() +
  geom_point(data=hasle.scores,aes(x=NMDS1,y=NMDS2,fill=lat),size=6, shape=21, colour="black") +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=station,group=station),alpha=0.30) + # add the convex hulls 
  coord_equal() +
  scale_fill_gradientn(colours = c(mypal3[1], "white", mypal3[3]),
                       values = scales::rescale(c(-55, -61, -62.5, -63.5, -70))) +
  # scale_x_continuous(limits = c(-1.1,1.2)) +
  # scale_y_continuous(limits = c(-1,1)) + 
  # theme(legend.position="none") +
  theme_linedraw(30) +
  theme(legend.key.size = unit(2, "cm")) +
  # theme(plot.title = element_text(size = 12, face = "bold"),
  #       legend.title=element_text(size=10), 
  #       legend.text=element_text(size=9)) +
  # theme(axis.line = element_line(colour = 'black', size = 1)) +
  annotate("text", x=0.64, y=0.59, label=paste("2d stress: ", round(hasle_NMDS$stress,3), sep=""))
  
png("hasle_NMDS.png", res=300, width=3000, height=2000)
p1
dev.off()


#########################
###-------GRIG-------###
#########################
  
#collapse to genus
grig_gen <- aggregate(count ~ project + site + genus + datetime, data=grig_data, FUN=sum)  

#summarise sums by sample
grig_sums <- grig_gen %>% group_by(project, site, datetime) %>% summarise(valve_sum = sum(count))
#summarise for means and sd
grig_sums2 <- grig_sums %>% group_by(project, site) %>% summarise(mean = mean(valve_sum), sd = sd(valve_sum))
#add zone
grig_sums2 <- left_join(grig_sums2, metadata[,c("site", "apf_side")], by="site")
#convert site names for graph
grig_sums2$site <- c("MS2","MS3","MS4","MS5","MS6","MS7")

psum2 <- ggplot(data=grig_sums2, aes(x=site, y=mean, fill=apf_side)) + 
  geom_bar(width=.5, stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  xlab("Station") +
  ylab("Diatom valve abundance") +
  scale_fill_manual(values=c(mypal3)) +
  theme_linedraw(30)


#remove genus "diatoms" as this is the total of the others
grig_gen <- grig_gen %>% filter(genus!="Diatoms")
  
#create sum df
df1 <- grig_gen %>%
  group_by(datetime,site) %>%
  summarize(sum = sum(count))
#merge sums to hasle_gen
df2 <- left_join(grig_gen, df1, by=c("datetime", "site"))
#calculate rel_abundance and add as column
df2$rel_abund <- round(df2$count/df2$sum, 4)
#remove intermediate columns
grig_gen <- select(df2, c(project, site, datetime, genus, rel_abund))

#convert "Diatoms" to "unknown"

#remove all rel_abundances of 0.0000
grig_gen <- grig_gen %>% filter(rel_abund!=0)

#delete intermediate dfs
rm(df1)
rm(df2)
  

# df_0.001_rel
grig_gen_0.001 <- grig_gen %>% filter(rel_abund>0.001)

#add lat long
grig_gen$site <- as.character(grig_gen$site)
grig_gen <- left_join(grig_gen, metadata[c("site", "lat", "long")], by="site")

#summarize genus count per site
richness_grig <- grig_gen %>%
  group_by(datetime, site) %>%
  count()
names(richness_grig) <- c("datetime", "site", "richness")

richness_grig <- left_join(richness_grig, metadata[,c("site", "lat")], by="site")


#create wide rel_abund df
grig_gen_wide <- spread(grig_gen, key=genus, value=rel_abund)
grig_gen_wide[is.na(grig_gen_wide)] <- 0

#calculate diversity and evenness
div_grig <- diversity(as.matrix(grig_gen_wide[,6:ncol(grig_gen_wide)]), "shannon")
s_grig <- specnumber(as.matrix(grig_gen_wide[,6:ncol(grig_gen_wide)])) ## rowSums(BCI > 0) does the same...
eve_grig <- div_grig/log(s_grig)

#create df with all metrics and site and project; to rbind with other projects downstream
grig_metrics <- data_frame(project = grig_gen_wide$project, 
                            site = as.character(grig_gen_wide$site),
                            datetime = grig_gen_wide$datetime,
                            depth = NA,
                            div = div_grig, 
                            eve = eve_grig,
                            rich = richness_grig$richness)
#add latitude from metadata
grig_metrics <- left_join(grig_metrics, metadata[,c("site", "lat")], by="site")



###---GRIG_NMDS---###

grig_matrix <- as.matrix(grig_gen_wide[,6:ncol(grig_gen_wide)])
grig_NMDS=metaMDS(grig_matrix,k=2,trymax=100)
grig.scores <- as.data.frame(scores(grig_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
grig.scores$site <- grig_gen_wide$site
grig.scores$lat <- grig_gen_wide$lat
grig.scores$long <- grig_gen_wide$long
grig.scores$datetime <- grig_gen_wide$datetime



grp.a <- grig.scores[grig.scores$site == "NBP96-04A_MS2", ][chull(grig.scores[grig.scores$site == 
                                                                             "NBP96-04A_MS2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- grig.scores[grig.scores$site == "NBP96-04A_MS3", ][chull(grig.scores[grig.scores$site == 
                                                                             "NBP96-04A_MS3", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.c <- grig.scores[grig.scores$site == "NBP96-04A_MS4", ][chull(grig.scores[grig.scores$site == 
                                                                             "NBP96-04A_MS4", c("NMDS1", "NMDS2")]), ]  # hull values for grp C
grp.d <- grig.scores[grig.scores$site == "NBP96-04A_MS5", ][chull(grig.scores[grig.scores$site == 
                                                                             "NBP96-04A_MS5", c("NMDS1", "NMDS2")]), ]  # hull values for grp D
grp.e <- grig.scores[grig.scores$site == "NBP96-04A_MS6", ][chull(grig.scores[grig.scores$site == 
                                                                             "NBP96-04A_MS6", c("NMDS1", "NMDS2")]), ]  # hull values for grp D
grp.f <- grig.scores[grig.scores$site == "NBP96-04A_MS7A", ][chull(grig.scores[grig.scores$site == 
                                                                             "NBP96-04A_MS7A", c("NMDS1", "NMDS2")]), ]  # hull values for grp D

hull.data <- rbind(grp.a, grp.b, grp.c, grp.d, grp.e, grp.f)  #combine grp.a and grp.b
hull.data




###---GRIG_PLOTS---###

p2 <- ggplot() +
  geom_point(data=grig.scores,aes(x=NMDS1,y=NMDS2,fill=lat),size=6, shape=21) +
  # [grig.scores$datetime=="1996-11-28",]
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=station,group=station),alpha=0.30) + # add the convex hulls 
  geom_polygon(data=hull.data, aes(x=NMDS1,y=NMDS2,group=site),alpha=0.20) + # add the convex hulls 
  coord_equal() +
  # scale_x_continuous(limits = c(-1.1,1.2)) +
  # scale_y_continuous(limits = c(-1,1)) + 
  # theme(legend.position="none") +
  theme_linedraw(30) +
  theme(legend.key.size = unit(2, "cm")) + 
  scale_fill_gradientn(colours = c(mypal3[1], "white", mypal3[3]),
                     values = scales::rescale(c(-55, -57.5, -60, -70, -76.5))) +
# theme(plot.title = element_text(size = 12, face = "bold")
#       legend.title=element_text(size=10), 
#       legend.text=element_text(size=9)) +
# theme(axis.line = element_line(colour = 'black', size = 1)) +
annotate("text", x=0.31, y=0.65, label=paste("2d stress: ", round(grig_NMDS$stress,3), sep=""))

png("grig_NMDS.png", res=300, width=3000, height=2000)
p2
dev.off()




#########################
###-------ASSMY-------###
#########################  

#collapse to genus
assmy_gen <- aggregate(count ~ project + site + genus + depth, data=assmy_data, FUN=sum)  
head(assmy_gen)

assmy_sums <- assmy_gen %>% group_by(project, site, depth) %>% summarise(valve_sum = sum(count))
head(assmy_sums)
#summarise for means and sd
assmy_sums2 <- assmy_sums %>% group_by(project, site) %>% summarise(mean = mean(valve_sum), sd = sd(valve_sum))
#add zone
assmy_sums2 <- left_join(assmy_sums2, metadata[,c("site", "apf_side")], by="site")
#convert site names for graph
assmy_sums2$site <- c("A12", "B42", "C48" ,"D90", "E108")

psum3 <- ggplot(data=assmy_sums2, aes(x=site, y=mean, fill=apf_side)) + 
  geom_bar(width=.5, stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  xlab("Station") +
  ylab("Diatom valve abundance") +
  scale_fill_manual(values=c(mypal3)) +
  theme_linedraw(30)



#create sum df
df1 <- assmy_gen %>%
  group_by(site, depth) %>%
  summarize(sum = sum(count))
#merge sums to hasle_gen
df2 <- left_join(assmy_gen, df1, by=c("depth", "site"))
#calculate rel_abundance and add as column
df2$rel_abund <- round(df2$count/df2$sum, 4)
#remove intermediate columns
assmy_gen <- select(df2, c(project, site, depth, genus, rel_abund))

#remove all rel_abundances of 0.0000
assmy_gen <- assmy_gen %>% filter(rel_abund!=0)

#delete intermediate dfs
rm(df1)
rm(df2)

# df_0.001_rel
assmy_gen_0.001 <- assmy_gen %>% filter(rel_abund>0.001)

# #add lat long
assmy_gen$site <- as.character(assmy_gen$site)
assmy_gen <- left_join(assmy_gen, metadata[c("site", "lat", "long")], by="site")

#summarize genus count per site
richness_assmy <- assmy_gen %>%
  group_by(site, depth) %>%
  count()
names(richness_assmy) <- c("site", "depth", "richness")

richness_assmy <- left_join(richness_assmy, metadata[,c("site", "lat")], by="site")


#create wide rel_abund df
assmy_gen_wide <- spread(assmy_gen, key=genus, value=rel_abund)
assmy_gen_wide[is.na(assmy_gen_wide)] <- 0

#calculate diversity and evenness
div_assmy <- diversity(as.matrix(assmy_gen_wide[,6:ncol(assmy_gen_wide)]), "shannon")
s_assmy <- specnumber(as.matrix(assmy_gen_wide[,6:ncol(assmy_gen_wide)])) ## rowSums(BCI > 0) does the same...
eve_assmy <- div_assmy/log(s_assmy)

#create df with all metrics and site and project; to rbind with other projects downstream
assmy_metrics <- data_frame(project = assmy_gen_wide$project, 
                           site = as.character(assmy_gen_wide$site),
                           datetime = NA,
                           depth = assmy_gen_wide$depth,
                           div = div_assmy, 
                           eve = eve_assmy,
                           rich = richness_assmy$richness)
#add latitude from metadata
assmy_metrics <- left_join(assmy_metrics, metadata[,c("site", "lat")], by="site")



#combine abundance bar graphs
assmy_sums2
# grig_sums2$apf_side.y <- NULL
# grig_sums2$apf_side.x <- "apf_side"
names(hasle_sums)[3] <- "mean"
hasle_sums$sd <- NA
hasle_sums <- hasle_sums[,c("project", "site", "mean", "sd", "apf_side")]

all_sums <- rbind(assmy_sums2, hasle_sums, grig_sums2)


psum4 <- ggplot(data=all_sums, aes(x=site, y=mean, fill=apf_side)) + 
  geom_bar(width=.5, stat="identity") +
  facet_wrap(~ project, nrow=3, scales="free") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  xlab("Station") +
  ylab("Diatom valve abundance") +
  scale_fill_manual(values=c(mypal3)) +
  theme_linedraw(30) +
  theme(legend.position="none")


####----METRICS---###

metrics <- rbind(hasle_metrics, grig_metrics, assmy_metrics)
#gather metrics into tidy single column
metrics <- gather(metrics, key=metric, value=measurement, div, eve, rich)
metrics$metric <- as.factor(metrics$metric)
metrics$project <- as.factor(metrics$project)


###---metrics stats---###
#lm for all projects combined
lm1 <- lm(measurement ~ lat, data=metrics[metrics$metric=="div",])
summary(lm1)
lm2 <- lm(measurement ~ lat, data=metrics[metrics$metric=="eve",])
summary(lm2)
lm3 <- lm(measurement ~ lat, data=metrics[metrics$metric=="rich",])
summary(lm3)

#lms for all projects split up
lm4a <- lm(measurement ~ lat, data=metrics[metrics$metric=="div" & metrics$project=="Brategg",])
summary(lm4a)
lm4b <- lm(measurement ~ lat, data=metrics[metrics$metric=="div" & metrics$project=="AESOPS",])
summary(lm4b)
lm4c <- lm(measurement ~ lat, data=metrics[metrics$metric=="div" & metrics$project=="EISENEX",])
summary(lm4c)

lm5a <- lm(measurement ~ lat, data=metrics[metrics$metric=="eve" & metrics$project=="Brategg",])
summary(lm5a)
lm5b <- lm(measurement ~ lat, data=metrics[metrics$metric=="eve" & metrics$project=="AESOPS",])
summary(lm5b)
lm5c <- lm(measurement ~ lat, data=metrics[metrics$metric=="eve" & metrics$project=="EISENEX",])
summary(lm5c)

lm6a <- lm(measurement ~ lat, data=metrics[metrics$metric=="rich" & metrics$project=="Brategg",])
summary(lm6a)
lm6b <- lm(measurement ~ lat, data=metrics[metrics$metric=="rich" & metrics$project=="AESOPS",])
summary(lm6b)
lm6c <- lm(measurement ~ lat, data=metrics[metrics$metric=="rich" & metrics$project=="EISENEX",])
summary(lm6c)


###---metrics plots---###

levels(metrics$metric)[levels(metrics$metric)=="div"] <- "Shannon diversity"
levels(metrics$metric)[levels(metrics$metric)=="eve"] <- "Pielou's Evenness"
levels(metrics$metric)[levels(metrics$metric)=="rich"] <- "Species Richness"

m1 <- ggplot(data=metrics, aes(x=lat, y=measurement, colour=project), shape=21, fill="black") +
  geom_point(size=3) +
  scale_colour_manual(values=mypal3) +
  facet_wrap(~ metric, nrow=3, scales="free_y") +
  scale_x_continuous(breaks=seq(-75, -45, 5)) +
  theme_linedraw(30) +
  # theme_bw() +
  xlab("Latitude (dd)") +
  ylab("Index") + 
  theme(legend.position="none") + 
  geom_smooth(method = "lm", se = FALSE)
  


###---combined wide---###
#combine long forms, then covert to combined_wide form
gen_combined <- bind_rows(hasle_gen, grig_gen, assmy_gen)

#add rel_apf factor before spreading
gen_combined <- left_join(gen_combined, metadata[,c("lat", "apf_side")], by="lat")

gen_wide <- spread(gen_combined, key=genus, value=rel_abund)



#fill in datetimes
gen_wide[gen_wide$project=="Brategg",]$datetime <- "1947-12-13"
gen_wide[gen_wide$project=="EISENEX",]$datetime <- "2000-08-29"
gen_wide[gen_wide$project=="Brategg",]$depth <- "unknown"
gen_wide[gen_wide$project=="AESOPS",]$depth <- "unknown"

gen_wide[is.na(gen_wide)] <- 0



#combined NMDS
combo_matrix <- as.matrix(gen_wide[,8:ncol(grig_gen_wide)])
all_NMDS=metaMDS(combo_matrix,k=2,trymax=100)
all.scores <- as.data.frame(scores(all_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
all.scores$site <- gen_wide$site
all.scores$project <- gen_wide$project
all.scores$lat <- gen_wide$lat
all.scores$long <- gen_wide$long
all.scores$datetime <- gen_wide$datetime
all.scores$rel_apf <- as.factor(gen_wide$apf_side)

grp.a <- all.scores[all.scores$project == "Brategg", ][chull(all.scores[all.scores$project == 
                                                                                "Brategg", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- all.scores[all.scores$project == "AESOPS", ][chull(all.scores[all.scores$project == 
                                                                                "AESOPS", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.c <- all.scores[all.scores$project == "EISENEX", ][chull(all.scores[all.scores$project == 
                                                                                "EISENEX", c("NMDS1", "NMDS2")]), ]  # hull values for grp C

hull.data <- rbind(grp.a, grp.b, grp.c)  #combine grps
hull.data


p3 <- ggplot() +
  geom_point(data=all.scores, aes(x=NMDS1,y=NMDS2,shape=rel_apf, colour=project), size=6) +#, show.legend = FALSE) +
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=project),alpha=0.30) + # add the convex hulls
  coord_equal() +
  # scale_x_continuous(limits = c(-1.1,1.2)) +
  # scale_y_continuous(limits = c(-1,1)) + 
  theme_linedraw(30) +
  theme(legend.key.size = unit(2, "cm")) + 
  scale_colour_manual(values=c(mypal3)) +
  # theme(plot.title = element_text(size = 12, face = "bold")
  #       legend.title=element_text(size=10), 
  #       legend.text=element_text(size=9)) +
  # theme(axis.line = element_line(colour = 'black', size = 1)) +
  annotate("text", x=-.65, y=0.9, label=paste("2d stress: ", round(grig_NMDS$stress,3), sep=""))








########-----HYPOTHESIS TESTING-----#######

#test for zone(factor) differences, both within project and for all project together
#permanova, permdisp
#simper, to demonstrate with species are most important

#dissimilarity matrices
dist_gen_wide <- vegdist(gen_wide[,8:ncol(gen_wide)], method="bray")
dist_test <- decostand(gen_wide[,8:ncol(gen_wide)])

#adonis, test by project factor
# perm_by_project <- adonis(formula = wisconsin(gen_wide[,8:ncol(gen_wide)]) ~ as.factor(gen_wide$project))
perm_by_project <- adonis(formula = dist_gen_wide ~ as.factor(gen_wide$project))

#check group dispersion
betadisper(dist_gen_wide, as.factor(gen_wide$project), type = "centroid")


#adonis, test by zone factor
perm_by_zone <- adonis(formula = dist_gen_wide ~ as.factor(gen_wide$apf_side))
#check group dispersion
betadisper(dist_gen_wide, as.factor(gen_wide$apf_side), type = "centroid")


#adonis, test by both factors
perm_by_zone <- adonis(formula = dist_gen_wide ~ as.factor(gen_wide$apf_side) + as.factor(gen_wide$project))


#simper for project factor
simp1 <- simper(wisconsin(gen_wide[,8:ncol(gen_wide)]) ~ as.factor(gen_wide$project))

  
#correlation between latitude and species richness?
lm1 <- lm(Shannon diversity ~ lat, data=metrics)
metrics




