# quality_assess.R

# load tidyverse package 
# and ggplot2
# patchwork
library(tidyverse)
library(ggplot2)
library(patchwork)

# Snakemake inputs/outputs
Var_Freq <- snakemake@input[["Var_Freq"]]
Var_Depth <- snakemake@input[["Var_Depth"]]
Var_Quality <- snakemake@input[["Var_Quality"]]
Var_s_miss <- snakemake@input[["Var_s_miss"]]
Ind_Depth <- snakemake@input[["Ind_Depth"]]

out_png <- snakemake@output[["png"]]

# MAF - minor allele freq analysis
# load the file
var_freq <- read_delim(Var_Freq, delim = "\t",
                       col_names = c("chr", "pos", "nalleles", 
                                     "nchr", "a1", "a2"), skip = 1)
# find the minor allele freq
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

# plot 
a <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a <- a + labs(title = "MAF Density") + theme_light() + theme(plot.title.position = "plot")

# summary of maf
summary(var_freq$maf)

# Per site mean depth
# load the file
var_depth <- read_delim(Var_Depth, delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", 
                                      "var_depth"), skip = 1)
# plot
b <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b <- b + labs(title = "Per-site Mean Depth") + theme_light() + xlim(0, 100) +
     theme(plot.title.position = "plot")

# summary of the mean depth per site
summary(var_depth$mean_depth) 

# Quality per site
# load the file
var_qual <- read_delim(Var_Quality, delim = "\t",
                       col_names = c("chr", "pos", 
                                     "qual"), skip = 1)

# plot
c <- ggplot(var_qual, aes(qual)) + 
  geom_density(fill = "dodgerblue1", colour = "black", 
               alpha = 0.3)
c <- c + labs(title = "Quality per Site")  + theme_light() + xlim(0, 6.556e+03) + 
     theme(plot.title.position = "plot")

sum(var_qual$qual<=30, na.rm = TRUE)

# summary
summary(var_qual$qual)

# misssing data per site
# load the file
var_miss <- read_delim(Var_s_miss, delim = "\t",
                       col_names = c("chr", "pos", 
                                     "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

# plot
d <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
d <- d + labs(title = "Missing Data per Site") + theme_light() + 
     theme(plot.title.position = "plot")

# summary
summary(var_miss$fmiss)

# individual depth
# load the file
ind_depth <- read_delim(Ind_Depth, delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# plot
e <- ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
e <- e + labs(title = "Individual Depth") + theme_light() + 
     theme(plot.title.position = "plot")

# summary
summary(ind_depth$depth)

# for combined visualization
combine <- (a | b | c)/(d | e)

# make the directory present
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

# save the png
ggsave(
  filename = out_png,       
  plot = combine + plot_annotation(title = "Overall Quality Summary"),
  width = 16,               
  height = 12,               
  dpi = 300)







