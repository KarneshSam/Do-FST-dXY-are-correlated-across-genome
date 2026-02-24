# fst_dxy.R

# load the neccessary libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggridges)
library(ggsci)

# Snakemake inputs/outputs
fst_file <- snakemake@input[["fst_file"]]
dxy_file <- snakemake@input[["dxy_file"]]

heatmap_png <- snakemake@output[["heatmap_png"]]
fst_gw_png <- snakemake@output[["fst_gw_png"]]
dxy_gw_png <- snakemake@output[["dxy_gw_png"]]
bar_png <- snakemake@output[["bar_png"]]

# fst
pixy_fst <- read.table(fst_file, header = T) %>% na.omit()
# dxy
pixy_dxy <- read.table(dxy_file, header = T) %>% na.omit()

# add the label as fst
# change the -ve to 0
# create the pop pairs
pixy_fst$data_type <- "fst"
pixy_fst <- mutate(pixy_fst,
                   avg_wc_fst = ifelse(avg_wc_fst < 0, 0, avg_wc_fst),
                   comparison = paste(pop1, pop2, sep = '_v_'))

# add the label as dxy
# change the -ve to 0
# create the pop pairs
pixy_dxy$data_type <- "dxy"
pixy_dxy <- mutate(pixy_dxy,
                   avg_dxy = ifelse(avg_dxy < 0, 0, avg_dxy),
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge data frames
pixy_dxy_sub <- pixy_dxy %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(pixy_dxy_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

pixy_fst_sub <- pixy_fst %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(pixy_fst_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# fix to get matching rows from both datasets
tmp_dxy <- semi_join(pixy_dxy_sub, pixy_fst_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(pixy_fst_sub, pixy_dxy_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))

# join the datasets together to create a single dataframe
data <- full_join(tmp_dxy, tmp_fst)

# add numbering to help with plotting
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = row_number())

# get position for vertical lines used in plot to delineate the linkage groups
data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

# summarise median data for Fst and Dxy
# based on the pop pairs and chr
data %>%
  group_by(comparison, chromosome, data_type) %>%
  summarise(median = median(value, na.rm = TRUE))

# plot heatmap
# for each population pairs
a <- ggplot(data, aes(x = comparison, y = window_pos_1, fill = value)) +
     geom_tile() +
     facet_grid(chromosome~data_type, scales = "free") +
     scale_fill_gradient2(low = "blue", high = "red", 
                       midpoint = 0.2) +
     xlab("Pop pairs") +
     ylab("Genomic Position (window)") +
     theme_bw() + ggtitle("Heatmap of FST and dXY Across Genome")

# genome wide analysis
plot_gw_fst <- function(pair, type){
  # extract for each pop pair
  tmp_data <- data %>% filter(data_type==type & comparison==pair )
  # plot the scatter plot
  # for each chromosome
  plot_gw <- ggplot(tmp_data, aes(window_pos_1, value, col = chromosome)) +
    geom_point(size=0.25) +
    facet_wrap(~chromosome, scales = "free_x") +
    scale_colour_manual(values = c("#1f78b4", "#33a02c")) +
    ylim(0,1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x="Genomic Position (window)", y=type)
  
  # plot the density
  # to know the median and how the distribution is
  plot_density <- ggplot(tmp_data, aes(value, chromosome, 
                                       fill=chromosome)) +
    geom_density_ridges(quantile_lines=TRUE, 
                        quantile_fun=function(x,...) median(x)) +
    theme_bw() + 
    theme(legend.position = "none", 
          axis.text.y=element_blank()) +
    facet_grid(comparison~.) +
    xlim(0,1) +
    scale_fill_npg() +
    labs(x=type, y="Density")
  # to combine the plot
  plot_gw + plot_density + plot_layout(widths = c(5, 1))
}

# make fst plots and assemble them
plot_8N_v_K0_fst <- plot_gw_fst("8N_v_K0", "fst")
plot_8N_v_Lesina_fst <- plot_gw_fst("8N_v_Lesina", "fst")
plot_K0_v_Lesina_fst <- plot_gw_fst("K0_v_Lesina", "fst")

b <- (plot_8N_v_K0_fst / plot_8N_v_Lesina_fst / plot_K0_v_Lesina_fst)

# make the dxy plot and assemble them
plot_8N_v_K0_dxy <- plot_gw_fst("8N_v_K0", "dxy")
plot_8N_v_Lesina_dxy <- plot_gw_fst("8N_v_Lesina", "dxy")
plot_K0_v_Lesina_dxy <- plot_gw_fst("K0_v_Lesina", "dxy")

c <- (plot_8N_v_K0_dxy / plot_8N_v_Lesina_dxy / plot_K0_v_Lesina_dxy)

# for correlation study
cor_data <- data %>%
  pivot_wider(
    id_cols = c(comparison, chromosome, window_pos_1, window_pos_2, position),
    names_from = data_type,
    values_from = value)

# using the linear model
# the correlation is analysed
# between the population
# across the genome window
model_fit <- ggplot(cor_data, aes(x=dxy, y=fst)) +
             geom_point(alpha=0.3) +
             geom_smooth(method="lm", color="red") +
             facet_wrap(~comparison, scales="free_x") +
             labs(x="dXY", y="FST") +
             theme_bw() + ggtitle("Fitting Linear model (dxy vs Fst)")

# to see the different range in the fst and dxy
# how it helps in correlation?
cor_classified <- cor_data %>%
  group_by(comparison) %>%
  mutate(
    fst_high = fst > quantile(fst, 0.90, na.rm = TRUE),
    dxy_high = dxy > quantile(dxy, 0.90, na.rm = TRUE))

# classified into: high - fst and dxy
# low - fst and dxy 
# under different combination
cor_classified <- cor_classified %>%
  mutate(category = case_when(
    fst_high & dxy_high ~ "HighFST_HighDXY",
    fst_high & !dxy_high ~ "HighFST_LowDXY",
    !fst_high & dxy_high ~ "LowFST_HighDXY",
    !fst_high & !dxy_high ~ "LowFST_LowDXY"))

# plot the different range of dxy and fst
# to know how it varies between population across the genome
Variation <- ggplot(cor_classified, aes(dxy, fst, color = category)) +
             geom_point(alpha = 0.4) +
             facet_wrap(~comparison) +
             theme_bw() + ggtitle("Variation of dxy vs Fst Across Populations") +
             guides(color = guide_legend(override.aes = list(size = 5)))

# correlation test
# create a data frame
# perform correlation based on the spearman correlation
corr_table <- cor_data %>%
  group_by(comparison, chromosome) %>%
  summarise(
    Spearman_rho = cor(fst, dxy, method = "spearman"),
    p_value = cor.test(fst, dxy, method = "spearman")$p.value,
    .groups = "drop")

# barplot for analysing the correlation 
# between the fst and dxy
# for each pop pair 
# also to know the support of fst and dxy corr on chromosomes
d <- ggplot(corr_table, aes(x=chromosome, y=Spearman_rho, fill=Spearman_rho)) +
     geom_bar(stat="identity") +
     scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
     facet_wrap(~ comparison) +
     labs(x="Chromosomes",
       y="Spearman correlation (FST vs dXY)",
       title="Correlation per chromosome within each population pair") +
     theme_bw()

# Save heatmap
dir.create(dirname(heatmap_png), recursive = TRUE, showWarnings = FALSE)
ggsave(heatmap_png, a, width=16, height=12, dpi=300)

# Save genome-wide FST plot
ggsave(fst_gw_png, b + plot_annotation(title = "Genome-wide Analysis Based on Fst"),
 width=16, height=12, dpi=300)

# Save genome-wide dXY plot
ggsave(dxy_gw_png, c + plot_annotation(title = "Genome-wide Analysis Based on dxy"),
 width=16, height=12, dpi=300)

# Save FST vs dXY correlation plot
ggsave(bar_png, d, width=16, height=12, dpi=300)


