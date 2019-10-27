# A naive approach to examining differences in DNA composition across genes
# =========================================================================

# Function A: correlationTable
# ----------------------------
# Input: 
# - directory // containing all sequences with FASTA header removed
# - nucleotideWidth // which indicates with length of subsequence used to examine frequencies
# Output:
# - A table of all pearson correlation coefficients across each gene pair
# =========================================================================

# Function B: heatmapPlot
# Input:
# - corrTable // input the saved correlationTable output variable into this call 
# Output:
# - A simple heatmap of the correlations

correlationTable <- function(directory, nucleotideWidth = 3) {
        
        setwd(directory)
        
        library(dplyr)
        library(data.table)
        library(stringr)
        library(arrangements)
        library(pheatmap)
        library(tibble)
        
        # 1. Build all possible nucleotide sequences of length k (nucleotideWidth); standard is 3 
        permute <- c("A", "G", "C", "T")
        permute <- data.table(permutations(permute, wordcount, replace = TRUE))
        permute$strings <- apply(permute[,1:wordcount] , 1, paste, collapse = "")
        
        # 2. Build a frequency table for each of the genes analyzed
        countTable <- data.table()
        for (i in list.files()) {
            data <- readr::read_file(i)
            data <- gsub("\n", "", data)
            geneLength <- nchar(data)
            for (x in permute$string) {
                countString <- str_count(data, x)
                tableString <- data.table(gene = i, 
                                          string = x, 
                                          count = countString, 
                                          length = geneLength,
                                          frequency = countString/geneLength)
                countTable <- rbind(countTable, tableString)
            }
        }

        # 3. Build a correlation table for each of the genes analyzed
        corrCombos <- t(combn(list.files(), 2))
        corrTable <- data.table()
        for (i in 1:dim(corrCombos)[1]) {
          freqA <- filter(countTable, gene == corrCombos[[i, 1]])$frequency
          freqB <- filter(countTable, gene == corrCombos[[i, 2]])$frequency
          correlationValue <- cor.test(freqA, freqB, method = "pearson")
          correlationValue <- correlationValue$estimate[[1]]
          outputValue <- data.table(geneA = corrCombos[i,1],
                                     geneB = corrCombos[i,2],
                                     pearson = correlationValue)
          corrTable <- rbind(corrTable, outputValue)
        }
        
        addReverse <- data.table(geneA = corrTable$geneB, 
                                 geneB = corrTable$geneA, 
                                 pearson = corrTable$pearson)
        
        corrTable <- rbind(corrTable, addReverse)
        
        return(corrTable)
}

correlationHeatmap <- function(corrTable) {
        heatmapMatrix <- dcast(corrTable, geneA ~ geneB, value.var = "pearson", fill = 1)
        heatmapMatrix <- column_to_rownames(heatmapMatrix, "geneA")
        heatmapMatrix <- as.matrix(heatmapMatrix)
        
        heatmapPlot <- pheatmap(heatmapMatrix, 
                                cluster_rows = FALSE, 
                                cluster_cols = FALSE)
        return(heatmapPlot)
}

# Test

testTable <- correlationTable(directory = "~/Documents/github/compbio/horizontalGeneTransferDetection/sampleGenes", nucleotideWidth = 3)
correlationHeatmap(testTable)

# Supplementary
# =========================================================================
# Potential k-means and heirarchical clustering solutions
library(tidyverse)
library(cluster)
library(factoextra)

# K means clustering
distance <- get_dist(heatmapMatrix, method = "pearson")
fviz_dist(distance)
k2 <- kmeans(heatmapMatrix, centers = 3, nstart = 25)
fviz_cluster(k2, data = heatmapMatrix)

# Heirarchical clustering
build <- dist(heatmapMatrix, method = "euclidean")
buildClust <- hclust(build, method = "complete")
plot(buildClust, cex = 0.6, hang = -1)








