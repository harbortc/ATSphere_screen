# Playing around with haruhiko screen data
#13.01.18
#also try subsetting the newick file based on Haruhiko's strain IDs

rm(list = ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(multcompView)
library(phylobase)
library(ape)
library(phylocanvas)

screen <- read.delim("/Volumes/FLASHDRIVE/AT_Sphere/Cologne soil data_171219.txt")
screen <- screen[c(1:200),]
str(screen)


#Plot rescue activity - with stats
#Stats:
RA.aov <- aov(screen$Relative.plant.rescue.activity.compared.to.FeEDTA ~ screen$Yang_Bacteria_classify, data = screen)
RA.tHSD <- TukeyHSD(RA.aov, ordered = FALSE, conf.level = 0.95)
generate_label <- function(HSD, group){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[group]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(screen, group, function (x) max(fivenum(x$Relative.plant.rescue.activity.compared.to.FeEDTA))*1.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = group, sort = FALSE)
  
  return(labels.df)
}
Rescue.box <- ggplot(screen, aes(x=screen$Yang_Bacteria_classify, y=screen$Relative.plant.rescue.activity.compared.to.FeEDTA, fill = screen$Yang_Bacteria_classify)) +
  geom_boxplot() +
  geom_point() +
  geom_text(data = generate_label(RA.tHSD, "screen$Yang_Bacteria_classify"), aes(x = plot.labels, y = max(V1), label = labels, fill = NA)) +
  ggtitle("Rescue Activity by Phylum") +
  theme_tufte() +
  scale_y_continuous(name = "Percent Growth\nCompared to FeEDTA") +
  scale_x_discrete(name = "Phylum") +
  guides(fill = F) +
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title.x = element_text(size=12),
        axis.title.y = element_text(size = 12)) +
  scale_fill_brewer(palette = "Set1")

ggsave(filename = "/Volumes/FLASHDRIVE/AT_Sphere/boxplot.pdf", plot = Rescue.box)



#Histogram
histogram <- ggplot(screen, aes(x = screen$Relative.plant.rescue.activity.compared.to.FeEDTA)) +
  geom_histogram(binwidth = 10, center = 5) +
  ggtitle("Rescue activity - all strains") +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title.x = element_text(size=12),
        axis.title.y = element_text(size = 12)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  scale_x_continuous(name = "Percent Growth Compared to FeEDTA") +
  scale_fill_brewer(palette = "Set1")

ggsave(filename = "/Volumes/FLASHDRIVE/AT_Sphere/histogram.pdf", plot = histogram)



# Upload ATsphere strains tree (all strains from ATSphere website)
ATsphere.new <- readNewick("/Volumes/FLASHDRIVE/AT_Sphere/spp_tree.newick.txt")

#Create vector with names of strains Haruhiko used
Haruhiko_names <- as.character(screen$Strain.name) #strains used in Haru's screen

ATSphere_strains <- as.character(ATsphere.new@label) #character vector with ATSphere names


setdiff <- setdiff(ATSphere_strains, Haruhiko_names) #strains not in both lists(prune these)
### I think this loses some strain names in the screen that weren't in ATSphere strains...

pruned <- prune(ATsphere.new, tips.exclude = setdiff, trim.internal = TRUE)

#convert phylo4 to newick tree format
is.new <- as_tree(pruned)

write.csv(is.new, file = "/Volumes/FLASHDRIVE/AT_Sphere/screened_tree.newick", quote = F, row.names = F)
#used this file to make tree in iTOL
#Have to remove the "x" at the beginning 
