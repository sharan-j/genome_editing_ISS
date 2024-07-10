## Load packages
library(tidyverse)
library(reshape2)
library(tiff)
library(imager)

## ------------------------------------------------------------------------------------------------
## set barcodes used in the experiment
## ------------------------------------------------------------------------------------------------

reference = read.csv("barcodes/reference.csv")

## ------------------------------------------------------------------------------------------------
## Set path to the folder containing the files to be analysed, open the files and extract meta data
## ------------------------------------------------------------------------------------------------

setwd("/Volumes/G_PHARMA_DataGenEdit$/Users/janjuhas/ISS_2020-2022/analysis/images/")
samCond = "mq_liver_used/"

filenames = list.files(samCond, pattern = "_sep_dots.csv$", full.names = TRUE, recursive = TRUE)
rawRep_names = lapply(samCond, grep, x = filenames, value = TRUE)


rawRep = data.frame()
for (i in 1:lengths(rawRep_names)) {
  print(rawRep_names[[1]][i])
  rr.pre = read.csv(rawRep_names[[1]][i], sep = ",", header = T, stringsAsFactors = F)
  print(nrow(rr.pre))
  try(
    {
  rr.pre$subRegionType = str_extract(rawRep_names[[1]][i], "periportal|pericentral|mid")
  rawRep = rbind(rawRep, rr.pre)
  })
}

 rm(rr.pre)


rawRep$injection_type = str_extract(rawRep$PathName_MAXb1, "48h|2w|coinjection|Mq")
rawRep$tissue_location = str_extract(rawRep$PathName_MAXb1, "RL|LML")

## ------------------------------------------------------------------------------------------------

chred = rawRep %>% 
  select(contains("Intensity_UpperQuartileIntensity_enhanced"), Location_Center_X, Location_Center_Y, Number_Object_Number, subRegionType, tissue_location, injection_type, FileName_MAXb1) %>%
  mutate(rep_num = str_extract(as.character(FileName_MAXb1), "rep[0-9]*|bact_[0-9]*|liver[0-9]|Mq[0-9][0-9][0-9]|201|202")) %>%
  mutate(region_num = str_extract(as.character(FileName_MAXb1), "reg[0-9]*")) %>%
  mutate(ID = paste(injection_type, rep_num, region_num,subRegionType, Number_Object_Number, sep = "_")) %>%
  unite(parentID, c(injection_type, rep_num, region_num),sep = "_", remove = F )%>%
  unite(normalisationID, c(injection_type, tissue_location),sep = "_", remove = F )

head(chred)

##-------------------------------------------------------------------------------------------------------------------
chred.norm = as.data.frame(chred)
##-------------------------------------------------------------------------------------------------------------------

chred.norm$s1 = apply(chred.norm[, endsWith( names(chred.norm), "_b1"   )], 1, sum)
chred.norm$s2 = apply(chred.norm[, endsWith( names(chred.norm), "_b2"   )], 1, sum)
chred.norm$s3 = apply(chred.norm[, endsWith( names(chred.norm), "_b3"   )], 1, sum)

chred.norm$m1 = apply(chred.norm[, endsWith( names(chred.norm), "_b1"   )], 1, max)
chred.norm$m2 = apply(chred.norm[, endsWith( names(chred.norm), "_b2"   )], 1, max)
chred.norm$m3 = apply(chred.norm[, endsWith( names(chred.norm), "_b3"   )], 1, max)

chred.norm$q1 = chred.norm$m1/chred.norm$s1
chred.norm$q2 = chred.norm$m2/chred.norm$s2
chred.norm$q3 = chred.norm$m3/chred.norm$s3

chred.norm$b1.max = names(chred.norm[grep("_b1", colnames(chred.norm))])[max.col(replace(chred.norm[grep("_b1", colnames(chred.norm))], is.na(chred.norm[grep("_b1", colnames(chred.norm))]), -999), 'first')*NA^(rowSums(chred.norm[grep("_b1", colnames(chred.norm))])==0)]
chred.norm$b2.max = names(chred.norm[grep("_b2", colnames(chred.norm))])[max.col(replace(chred.norm[grep("_b2", colnames(chred.norm))], is.na(chred.norm[grep("_b2", colnames(chred.norm))]), -999), 'first')*NA^(rowSums(chred.norm[grep("_b2", colnames(chred.norm))])==0)]
chred.norm$b3.max = names(chred.norm[grep("_b3", colnames(chred.norm))])[max.col(replace(chred.norm[grep("_b3", colnames(chred.norm))], is.na(chred.norm[grep("_b3", colnames(chred.norm))]), -999), 'first')*NA^(rowSums(chred.norm[grep("_b3", colnames(chred.norm))])==0)]

chred.norm$b1.max.base = chred.norm$b1.max
chred.norm$b1.max.base = tolower(chred.norm$b1.max.base)
chred.norm$b1.max.base = gsub(".*_red_b1", "T", chred.norm$b1.max.base)
chred.norm$b1.max.base = gsub(".*_farred_b1", "A", chred.norm$b1.max.base)
chred.norm$b1.max.base = gsub(".*_gold_b1", "G", chred.norm$b1.max.base)
chred.norm$b1.max.base = gsub(".*_cy7_b1", "C", chred.norm$b1.max.base)

chred.norm$b2.max.base = chred.norm$b2.max
chred.norm$b2.max.base = tolower(chred.norm$b2.max.base)
chred.norm$b2.max.base = gsub(".*_red_b2", "T", chred.norm$b2.max.base)
chred.norm$b2.max.base = gsub(".*_farred_b2", "A", chred.norm$b2.max.base)
chred.norm$b2.max.base = gsub(".*_gold_b2", "G", chred.norm$b2.max.base)
chred.norm$b2.max.base = gsub(".*_cy7_b2", "C", chred.norm$b2.max.base)

chred.norm$b3.max.base = chred.norm$b3.max
chred.norm$b3.max.base = tolower(chred.norm$b3.max.base)
chred.norm$b3.max.base = gsub(".*_red_b3", "T", chred.norm$b3.max.base)
chred.norm$b3.max.base = gsub(".*_farred_b3", "A", chred.norm$b3.max.base)
chred.norm$b3.max.base = gsub(".*_gold_b3", "G", chred.norm$b3.max.base)
chred.norm$b3.max.base = gsub(".*_cy7_b3", "C", chred.norm$b3.max.base)


chred.norm$str.call = paste0(chred.norm$b1.max.base, chred.norm$b2.max.base, chred.norm$b3.max.base)
plot(table(chred.norm$str.call))
table(chred.norm$str.call)


### Clean df 
chred.norm.clean = chred.norm[complete.cases(chred.norm), ]
chred.norm.clean = chred.norm.clean[chred.norm.clean$q1 >= 0.35,]
chred.norm.clean = chred.norm.clean[chred.norm.clean$q2 >= 0.35,]
chred.norm.clean = chred.norm.clean[chred.norm.clean$q3 >= 0.35,]

table(chred.norm.clean$str.call)
plot(table(chred.norm.clean$str.call))

chred.norm.clean$str.call[which(grepl("TA.", chred.norm.clean$str.call))] = "t149_M"
chred.norm.clean$str.call[which(grepl(".TA", chred.norm.clean$str.call))] = "t149_W"
chred.norm.clean$str.call[which(grepl("AG.", chred.norm.clean$str.call))] = "d25_M"
chred.norm.clean$str.call[which(grepl("GA.", chred.norm.clean$str.call))] = "d25_W"

chred.norm.clean$gene = gsub("_.$", "", chred.norm.clean$str.call)

## ------------------------------------------------------------------------------------------------
## Calculate the editing efficiency 
## ------------------------------------------------------------------------------------------------

chred.clean.overallEditing = chred.norm.clean %>% 
  filter(str.call == "d25_M" | str.call == "d25_W" | str.call == "t149_M" | str.call == "t149_W" ) %>%
  group_by(injection_type, gene, rep_num, region_num, subRegionType) %>%
  count(str.call) %>% mutate(tot = sum(n)) %>% mutate(pct = (n*100)/tot) #%>% filter(tot >=10)

chred.clean.overallEditing$mutation_status = gsub("^.*_", "", chred.clean.overallEditing$str.call)
chred.clean.overallEditing.plot = chred.clean.overallEditing  %>% 
  filter(mutation_status == "M") %>% group_by(injection_type, gene) %>% mutate(means = mean(pct))
chred.clean.overallEditing.plot$strategy = ifelse( grepl("Mq2",chred.clean.overallEditing.plot$rep_num), "strategy2", "strategy1"  ) 
chred.clean.overallEditing.plot$subRegionType = factor(chred.clean.overallEditing.plot$subRegionType, levels = c("periportal", "mid", "pericentral", "NOTperiportal"))

# chred.clean.overallEditing.plot %>% filter(mutation_status == "M") %>%  
#   ggplot(aes(subRegionType, pct, fill = subRegionType)) +
#   stat_summary(fun = "mean", geom = "bar", color = "darkgrey", position = "dodge") +
#   geom_jitter( width = 0.1) +
#   theme_classic() +   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   scale_fill_manual(values = c("indianred1", "indianred3", "indianred4")) +
#   ylab("Editing efficiency (%)") + xlab("") + 
#   ylim(0,100) +
#   facet_grid(cols = vars(gene), rows = vars(strategy))


