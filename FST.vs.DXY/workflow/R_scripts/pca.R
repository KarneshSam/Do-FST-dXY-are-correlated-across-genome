# pca.R

# load tidyverse package
library(tidyverse)
library(ggplot2)

# Get Snakemake inputs and outputs
pca_vec_file <- snakemake@input[["pca_vec"]]
pca_val_file <- snakemake@input[["pca_val"]]
output_png <- snakemake@output[["png"]]

# read in data
pca <- read_table(pca_vec_file, col_names = FALSE)
eigenval <- scan(pca_val_file)

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
# set the names from PC1..
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual 
# create a pop vector
pop <- rep(NA, length(pca$ind))
pop[grep("^8N", pca$ind)] <- "Pop1"
pop[grep("^Lesina", pca$ind)] <- "Pop2"
pop[grep("^K0", pca$ind)] <- "Pop3"
pop[grep("^Naxos", pca$ind)] <- "Outgroup"  

# location
# condider each pop from one geometric coordinate
loc <- rep(NA, length(pca$ind))
loc[grep("^Lesina", pca$ind)] <- "Lesina"       
loc[grep("^K0", pca$ind)] <- "K0"           
loc[grep("^8N", pca$ind)] <- "8N"       
loc[grep("^Naxos", pca$ind)] <- "Naxos" 

# remake data.frame
pca <- as.tibble(data.frame(pca, pop, loc))
# first convert to percentage variance explained
pve <- data.frame(PC = 1:16, pve = eigenval/sum(eigenval)*100)

# make plot
# for distribution of PCs
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a <- a + ylab("Percentage variance explained") + 
     ggtitle("Proportion of Variance per Principal Component") +
     theme_light()

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = pop, shape = loc)) + 
  geom_point(size = 4)
b <- b + scale_colour_manual(values = c("red", "blue", "green", "yellow"))
b <- b + coord_equal() + theme_light()
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
     ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
     ggtitle("Population Structure Visualization via PCA") + 
     theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# make sure the directory present
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)

# save the plot
ggsave(output_png, b, width=12, height=8, limitsize = FALSE)
