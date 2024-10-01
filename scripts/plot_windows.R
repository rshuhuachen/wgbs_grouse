### load packages ###
pacman::p_load(tidyverse, cowplot)

### load theme ###
source("/home/nioo/rebeccash/PhD_grouse/epi_lekking/scripts/plotting_theme.R")

### load data ####
cpg_windows_min2_sum <- read.csv(file = "output/summary_cpg_meth_min2.csv")

### plot ###
# assign colours

clr_4 <- c(clrs[1], clrs[2], clrs[8], clrs[13])

# mean TSS methylation (for y value label)

mean_tss_meth <- mean(cpg_windows_min2$methperc_mean[which(cpg_windows_min2$region == "TSS")], na.rm=T)      

ggplot(cpg_windows_min2_sum, aes(x = window_total, y = mean_meth)) + 
  geom_ribbon(aes(ymin = mean_meth - se_meth, ymax =mean_meth + se_meth), fill = clrs[3], alpha = 0.7)+
  geom_line(col = "black", linewidth = 0.8) +
  geom_point() + labs(y = "Mean CpG methylation %", title = "CpG methylation across gene regions")+
  annotate("text", label="10 kb upstream", x = 20, y = 66, col = clr_4[1], size = 6)+
  annotate("text", label="Gene body", x = 62, y = 66, col = clr_4[3], size = 6)+
  annotate("text", label="10 kb downstream", x = 104, y = 66, col = clr_4[4], size = 6)+
  annotate("text", label="TSS", x = 51, y = mean_tss_meth, col = clr_4[2], size = 6)+ 
  geom_segment(aes(xend = 43, y = mean_tss_meth, x = 47, yend = mean_tss_meth), arrow = arrow(length = unit(0.2, "cm")))+
  geom_vline(xintercept = 41.5, linetype = "dotted", col = clr_4[1])+ 
  geom_vline(xintercept = 42.5, linetype = "dotted", col = clr_4[2])+ 
  geom_vline(xintercept = 41.5, linetype = "dotted", col = clr_4[3])+ 
  geom_vline(xintercept = 83.5, linetype = "dotted", col = clr_4[4])+ 
  scale_color_manual(values=clr_4)+
  ylim(55, 67)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") -> across_genes_min2

ggsave(across_genes_min2, file = "plots/test.png", width=12, height=12)

tss_min2 <- subset(cpg_windows_min2, region == "TSS")
tss_min2$methperc_n <- as.numeric(tss_min2$methperc_n)

ggplot(tss_min2, aes(x = methperc_mean)) + geom_histogram(fill=) +
  labs(x = "Mean CpG methylation %", title = "Histogram TSS CpG methylation % ", subtitle = "(min sites = 2)", y = "Count") -> tss_hist_min2

ggplot(tss_min2, aes(x = methperc_n, y = methperc_mean)) + geom_point() +
  labs(x = "Number of CpGs", title = "CpG sites vs mean methylation in TSS", subtitle = "(min sites = 2)", y = "Mean CpG methylation %") +
  ylim(0, 100)+
  geom_smooth(method="lm")-> tss_n_vs_meth_min2

cowplot::plot_grid(tss_hist_min2, tss_n_vs_meth_min2, ncol = 2, align = "hv", axis = "bl") -> tss_a_min2
cowplot::plot_grid(tss_a_min2, across_genes_min2, ncol = 1, align = "h", axis = "l") -> tss_sum_min2

ggsave(tss_sum_min2, file = "plots/summary_wgbs_TSS_min2.png", width=16, height=16)

### or just the histogram and across the genes

cowplot::plot_grid(tss_hist_min2, across_genes_min2, ncol = 1, align = "h", axis = "l") -> hist_across_min2

ggsave(hist_across_min2, file = "plots/summary_wgbs_TSS_min2_hist_across.png", width=16, height=16)

## min 3

tss_min3 <- subset(cpg_windows_min3, region == "TSS")
tss_min3$methperc_n <- as.numeric(tss_min3$methperc_n)

ggplot(tss_min3, aes(x = methperc_mean)) + geom_histogram() +
  labs(x = "Mean CpG methylation %", title = "Histogram TSS CpG methylation %" , subtitle = "(min sites = 3)", y = "Count")-> tss_hist_min3

ggplot(tss_min3, aes(x = methperc_n, y = methperc_mean)) + geom_point() +
  labs(x = "Number of CpGs", title = "#CpG sites vs mean methylation in TSS", subtitle = "(min sites = 3)",y = "Mean CpG methylation %") +
  ylim(0, 100)+
  geom_smooth(method="lm")-> tss_n_vs_meth_min3

mean_tss_meth <- mean(cpg_windows_min3$methperc_mean[which(cpg_windows_min3$region == "TSS")], na.rm=T)      

ggplot(cpg_windows_min3_sum, aes(x = window_total, y = mean_meth)) + 
  geom_ribbon(aes(ymin = mean_meth - se_meth, ymax =mean_meth + se_meth), fill = "lightblue", alpha = 0.7)+
  geom_line(col = "black", linewidth = 0.8) +
  geom_point() + labs(y = "Mean CpG methylation %", title = "CpG methylation across gene regions")+
  geom_text(label="10 kb upstream", x = 20, y = 67, col = clr_4[1], family = "Arial", size = 6)+
  geom_text(label="Gene body", x = 62, y = 67, col = clr_4[3], family = "Arial", size = 6)+
  geom_text(label="10 kb downstream", x = 104, y = 67, col = clr_4[4], family = "Arial", size = 6)+
  geom_text(label="TSS", x = 51, y = mean_tss_meth, col = clr_4[2], family = "Arial", size = 6)+ 
  geom_segment(aes(xend = 43, y = mean_tss_meth, x = 47, yend = mean_tss_meth), arrow = arrow(length = unit(0.2, "cm")))+
  geom_vline(xintercept = 41.5, linetype = "dotted", col = clr_4[1])+ 
  geom_vline(xintercept = 42.5, linetype = "dotted", col = clr_4[2])+ 
  geom_vline(xintercept = 41.5, linetype = "dotted", col = clr_4[3])+ 
  geom_vline(xintercept = 83.5, linetype = "dotted", col = clr_4[4])+ 
  scale_color_manual(values=clr_4)+
  ylim(52, 68)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") -> across_genes_min3

cowplot::plot_grid(tss_hist_min3, tss_n_vs_meth_min3, ncol = 2, align = "hv", axis = "bl") -> tss_a_min3
cowplot::plot_grid(tss_a_min3, across_genes_min3, ncol = 1, align = "h", axis = "l") -> tss_sum_min3

ggsave(tss_sum_min3, file = "plots/summary_wgbs_TSS_min3.png", width=16, height=16)

#without coverage

cowplot::plot_grid(tss_hist_min3, across_genes_min3, ncol = 1, align = "h", axis = "l") -> tss_sum_min3_v2

ggsave(tss_sum_min3_v2, file = "plots/summary_wgbs_TSS_min3_v2.png", width=10, height=10)
