# Creates stacked barplots with clones/rearrangements sorted from largest at
# top to smallest at bottom. Size of each bar segment represents the size of
# the clone/rearrangement. The most abundant amino acids in the specified column
# will be colored in each barplot to show the size and abundance within each
# sample.


stacked <- function(data, minimum = 50, maximum = 102, pull = 1, border = "grey30"){

 
# Minimum and maximum alter the scale of the y axis to focus on a certain range of
# percentages.
# Pull number gives column that the data is sorted by before aa are pulled from.
# 1 for "sum of all templates", 2 for "present in", 3 for first sample, 4 for second
# sample, etc.

 
library(lattice)
library(tidyr)
 

# Sort based on pull value and identify the most abundant aa sequences
temp <- data
order.temp <- temp[order(temp[,(pull+1)], decreasing = TRUE),]

aa1 <- order.temp[1,1]
aa2 <- order.temp[2,1]
aa3 <- order.temp[3,1]
aa4 <- order.temp[4,1]
aa5 <- order.temp[5,1]
aa6 <- order.temp[6,1]
aa7 <- order.temp[7,1]
aa8 <- order.temp[8,1]
aa9 <- order.temp[9,1]
aa10 <- order.temp[10,1]

 
# Isolate data needed by removing the sum and present in columns
df <- data[,c(1,5:ncol(data))]
df$Lineage <- factor(df$Lineage)


# Change the data frame to long form
df_long <- gather(df, sample, sequences, 2:ncol(df))
 

# Identify and change to percent of total templates

sampleNames <- unique(df_long$sample)

for(i in 1:length(sampleNames)) {

                df_long$sequences[df_long$sample == sampleNames[i]] <- df_long$sequences[df_long$sample == sampleNames[i]]/sum(df_long$sequences[df_long$sample == sampleNames[i]]) * 100
}


# Identify aa's to be colored
g <- order(df_long$sequences)
colVec <- ifelse(df_long$Lineage == aa1, "red",
                                ifelse(df_long$Lineage == aa2, "green",
                                ifelse(df_long$Lineage == aa3, "purple",
                                ifelse(df_long$Lineage == aa4, "blue",
                                ifelse(df_long$Lineage == aa5, "yellow",
                                ifelse(df_long$Lineage == aa6, "coral",
                                ifelse(df_long$Lineage == aa7, "pink",
                                ifelse(df_long$Lineage == aa8, "navy",
                                ifelse(df_long$Lineage == aa9, "darkgreen",
                                ifelse(df_long$Lineage == aa10, "maroon", "white"))))))))))[g]

 
# Plot the stacked barplots
barchart(sequences ~ sample, groups = factor(1:nrow(df_long),levels=g), stack = TRUE, data = df_long,  border = border, col=colVec, ylim = c(minimum,maximum), ylab = "Cumulative Percent of Templates")
#can use border = for border color change, box.width or box.ratio for width of bars
}