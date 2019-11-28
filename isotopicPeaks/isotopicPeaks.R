# This program uses three functions that can be piped to generate theoretic MS2 spectra


# Functions
# ============================================================
# ============================================================
# ============================================================

# Load libraries
library(data.table)
library(tidyverse)

# Read in sequence
sequence <- commandArgs(trailingOnly = TRUE)

# Generate list of all y-ions and b-ions
getIons <- function(sequence = sequence) {
      b_ions <- list()
      y_ions <- list()
      split <- strsplit(sequence, "")
      for (i in 1:(length(split[[1]])-1)) {
        b_ion <- split[[1]][1:i]
        y_ion <- split[[1]][(i+1):length(split[[1]])]
        b_ions <- append(b_ions, paste(b_ion, collapse = ""))
        y_ions <- append(y_ions, paste(y_ion, collapse = ""))
      }
      ions <- data.table(b_ions = unlist(b_ions), y_ions = unlist(y_ions))
      return(ions)
}

# Get masses for each peptide in both ion series
getMasses <- function(ion_list = getIonsOutput, PTM = "none") {
      mass_table <- data.table::fread("~/Documents/github/compbio/isotopicPeaks/aminoAcidMasses.csv")
      b_ions <- ion_list$b_ions
      y_ions <- ion_list$y_ions
      b_ion_masses <- c()
      y_ion_masses <- c()
      for (i in b_ions) {
        mass <- sum(subset(mass_table, 
                           oneLetter %in% strsplit(i, "")[[1]])$monoisotopic)
        b_ion_masses[i] <- mass
      }
      for (i in y_ions) {
        mass <- sum(subset(mass_table, 
                           oneLetter %in% strsplit(i, "")[[1]])$monoisotopic)
        y_ion_masses[i] <- mass
      }
      ion_list$b_ion_masses <- b_ion_masses
      ion_list$y_ion_masses <- y_ion_masses
      
      y_ion_list <- data.table(series = "y-ions",
                               ions = y_ions,
                               ion_masses = y_ion_masses)
      b_ion_list <- data.table(series = "b-ions",
                               ions = b_ions,
                               ion_masses = b_ion_masses)
      
      ion_list <- rbind(y_ion_list, b_ion_list)
      ion_list$`m/z` <- ion_list$ion_masses/2
      ion_list$PTM <- "None"
      
      if (PTM == "crotonylation") {
        crotonylation <- data.table()
        for (i in 1:dim(ion_list)[1]) {
          if ("K" %in% strsplit(ion_list$ions[i],"")[[1]]) {
            kcr <- ion_list[i,]
            kcr$ion_masses <- as.numeric(kcr$ion_masses) + 68/2
            kcr$PTM <- "Kcr"
            crotonylation <- rbind(crotonylation, kcr)
          } 
         }
        ion_list <- rbind(ion_list, crotonylation)
        }
      return(ion_list)
}

plotSpectra <- function(mass_table = masses) {
  library(ggplot2)
  library(ggrepel)
  ggplot(mass_table, aes(x = ion_masses, y = 1)) + 
    geom_col(aes(fill = PTM), width = 5) +
    geom_label_repel(aes(label = ions, color = series), size = 6, alpha = 0.9) +
    scale_y_continuous(expand = c(0.01,0)) +
    scale_fill_manual(values = c("red", "black")) +
    scale_color_manual(values = c("blue", "green3")) +
    theme_light() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18)) +
    ylab("") +
    xlab("m/z")
}

# Run the functions in a pipe
# ============================================================
# ============================================================
# ============================================================

plotThis <- getIons(sequence = "VALKGR") %>%
             getMasses() %>%
              plotSpectra()
plotThisKcr <- getIons(sequence = "VALKGR") %>%
             getMasses(PTM = "crotonylation") %>%
             plotSpectra()
