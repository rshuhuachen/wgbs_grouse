#### Here we will annotate the methylation call files, calculate mean methylation (CpG and non-CpG) per window, and plot the result ####

#### packages ####
pacman::p_load(data.table, genomation, GenomicFeatures, rtracklayer, GenomicRanges, windowscanr, tidyverse)

#### load data ####
meth <- fread("data/processed/methcall/methcall_report.CX_report.txt.gz")

meth$V1 <- gsub(";", "__", meth$V1)
meth$V1 <- gsub("=", "_", meth$V1)

names(meth) <- c("chr", "start", "strand", "nC", "nT", "dicontext", "tricontext")
meth$end <- meth$start

cpg <- subset(meth, dicontext == "CG")
ch <- subset(meth, dicontext != "CG")
rm(meth)

##### Annotation ####
### function to annotate dataframe

annotate_methfile <- function(long_df, annotation_dir, regions){
  pacman::p_load(genomation, GenomicFeatures, rtracklayer, GenomicRanges, tidyverse)
  
  ### annotate df accordingly and add region name
  long_df$end <- long_df$start
  df_gr <- as(long_df, "GRanges")
  
  #more complicated than just subsetting, want to include gene information!
  annotated <- data.frame()
  
  for (i in regions){
    #load in gff file
    region <- unique(gffToGRanges(paste0(annotation_dir, i, ".gff3")))
    #find overlaps in df and gene region, because can have multiple hits include all
    hits <- findOverlaps(df_gr, region)
    idx_c <- queryHits(hits)
    idx_region <- subjectHits(hits)
    
    #get values from overlaps for cpg and region separately
    values_region <- data.frame(chr=region@seqnames[idx_region],
                                start = region@ranges@start[idx_region], 
                                end= region@ranges@start[idx_region] + region@ranges@width[idx_region], 
                                gene_id = mcols(region)$gene_id[idx_region],
                                region = i)
    
    values_cytosines <- data.frame(pos = df_gr@ranges@start[idx_c], 
                                   nC = mcols(df_gr)$nC[idx_c], 
                                   nT = mcols(df_gr)$nT[idx_c])
    
    annotated_region <- cbind(values_cytosines, values_region) #combine
    annotated <- rbind(annotated, annotated_region) #add to full df
  }
  annotated$region <- as.factor(annotated$region)
  
  return(annotated)
}

## Annotate CpG and CH file
cpg_annotate <- annotate_methfile(long_df = cpg, 
                                  annotation_dir = "/vol/cluster-data/rchen/git/grouse-annotation/output/", 
                                  regions = c("TSS", "genes", "upstream", "downstream"))

ch_annotate <- annotate_methfile(long_df = ch, 
                                  annotation_dir = "/vol/cluster-data/rchen/git/grouse-annotation/output/", 
                                  regions = c("TSS", "genes", "upstream", "downstream"))

cpg_annotate$methperc <- cpg_annotate$nC / (cpg_annotate$nC + cpg_annotate$nT) * 100

save(cpg_annotate, file = "data/processed/annotated_cpg.RData")
save(ch_annotate, file = "data/processed/annotated_chh.RData")
#### Per gene, calculate meth% per bin
#cpg_annotate <- cpg_annotate[,c(1:5,7:9,11)]
cpg_list <- cpg_annotate %>% group_split(gene_id, region)

save(cpg_list, file = "data/processed/annotated_cpg_list.RData")

meth_per_window <- function(gene_region){tryCatch({
  #determine start and end of each region
 
    #determine window length for each region
  if (gene_region$region[1] == "TSS"){
    start_region <- gene_region$start[1]
    end_region <- gene_region$end[1]
    length = end_region - start_region
    overlap = 0
    } else if (gene_region$region[1] == "genes"){
    start_region <- gene_region$start[1]
    end_region <- gene_region$end[1]
    length = floor((end_region - start_region) / 40)
    overlap = 250 #maybe change this to 0.5 * length?
    } else if (gene_region$region[1] == "downstream"){
    start_region <- gene_region$start[1]
    end_region <- gene_region$end[1]
    length = 250
    overlap = 125 
    } else if (gene_region$region[1] == "upstream"){
    start_region <- gene_region$start[1]
    end_region <- gene_region$end[1] #exlude TSS
    length = 250
    overlap = 125}
  
  ### function window slider adapted from windowscanr
  
  source("scripts/function_sliding_window_custom.R")
  window <- win_scan(x = gene_region,
                     position = "pos",
                     values = "methperc",
                     win_size = length,
                     win_step = length,
                     overlap = overlap,
                     funs = c("mean", "sd"),
                     region = as.character(gene_region$region[1]))
  window[window == "NaN"] <- NA
  sum_window <- data.frame(window)
  sum_window$chr<- gene_region$chr[1]
  sum_window$gene_id<- gene_region$gene_id[1]
  sum_window$region<- gene_region$region[1]
  sum_window$start_bp<- start_region
  sum_window$end_bp<- end_region
  sum_window$length_bin<- length
  sum_window$win_nr <- row.names(sum_window)
  
  return(sum_window)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(paste0(gene_region$chr[1], "_", gene_region$pos[1]))})}

cpg_windows <- parallel::mclapply(cpg_list, meth_per_window, mc.cores=12)
cpg_windows <- do.call(rbind.data.frame, cpg_windows)

save(cpg_windows, file = "data/processed/cpg_meth_per_window_per_region.RData")


#### plotting
#extract for one gene
cpg_test <- subset(cpg_windows, gene_id == "ANN00003" | gene_id == "ANN00004")

cpg_test$region <- factor(cpg_test$region, levels = c("upstream", "TSS", "genes", "downstream"))
cpg_test$win_nr <- as.numeric(cpg_test$win_nr)
cpg_test$methperc_mean <- as.numeric(cpg_test$methperc_mean)
cpg_test <- cpg_test %>% group_by(gene_id) %>% arrange(region, win_nr) %>% mutate(window_total = row_number())
#problem when gene is too short and don't have 41 windows

cpg_test_n <- cpg_test %>% group_by(window_total) %>% 
  summarise(n_meth = sum(!is.na(methperc_mean)))

cpg_test_mean <- cpg_test %>% group_by(window_total) %>% 
  summarise(mean_meth = mean(methperc_mean, na.rm=T))

cpg_test_sd <- cpg_test %>% group_by(window_total) %>% 
  summarise(sd_meth = sd(methperc_mean, na.rm=T))

cpg_test_sum <- left_join(cpg_test_n, cpg_test_mean, by = "window_total") #leave out region for now, so wrong - talk to kees
cpg_test_sum <- left_join(cpg_test_sum, cpg_test_sd, by = "window_total")
cpg_test_sum$se_meth <- cpg_test_sum$sd_meth/sqrt(cpg_test_sum$n_meth)
cpg_test_sum[cpg_test_sum == "NaN"] <- NA

# real data
cpg_windows$region <- factor(cpg_windows$region, levels = c("upstream", "TSS", "genes", "downstream"))
cpg_windows$win_nr <- as.numeric(cpg_windows$win_nr)
cpg_windows$methperc_mean <- as.numeric(cpg_windows$methperc_mean)
cpg_windows <- cpg_windows %>% group_by(gene_id) %>% arrange(region, win_nr) %>% mutate(window_total = row_number())

cpg_windows_n <- cpg_windows %>% group_by(window_total) %>% 
  summarise(n_meth = sum(!is.na(methperc_mean)))

cpg_windows_mean <- cpg_windows %>% group_by(window_total) %>% 
  summarise(mean_meth = mean(methperc_mean, na.rm=T))

cpg_windows_sd <- cpg_windows %>% group_by(window_total) %>% 
  summarise(sd_meth = sd(methperc_mean, na.rm=T))

cpg_windows_sum <- left_join(cpg_windows_n, cpg_windows_mean, by = "window_total") #leave out region for now, so wrong - talk to kees
cpg_windows_sum <- left_join(cpg_windows_sum, cpg_windows_sd, by = "window_total")
cpg_windows_sum$se_meth <- cpg_windows_sum$sd_meth/sqrt(cpg_windows_sum$n_meth)
cpg_windows_sum[cpg_windows_sum == "NaN"] <- NA

# plot 
# colours
pacman::p_load(MoMAColors, prismatic)
clr <- MoMAColors::moma.colors("Fritsch", 8) %>% color()
clr_4 <- clr[c(1,7,3,2)]

theme_set(theme_classic() + theme(title = element_text(size=16),
                                  plot.subtitle = element_text(size=14),
                                  axis.title = element_text(size = 18, family = "Arial"),
                                  axis.text = element_text(size = 16, family = "Arial"),
                                  text=element_text(size=14, family = "Arial"),
                                  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                                              color = "black"),
                                  plot.margin = margin(1,1,1,1, "cm")))
                                  
ggplot(cpg_windows_sum, aes(x = window_total, y = mean_meth)) + 
  geom_line(col = "black") +
  geom_ribbon(aes(ymin = mean_meth - se_meth, ymax =mean_meth + se_meth), fill = "lightblue", alpha = 0.8)+
  geom_point() + labs(y = "Mean CpG methylation %", title = "CpG methylation across gene regions")+
  geom_vline(xintercept = 41.5, linetype = "dotted", col = clr_4[1])+ 
  geom_vline(xintercept = 42.5, linetype = "dotted", col = clr_4[2])+ 
  geom_vline(xintercept = 41.5, linetype = "dotted", col = clr_4[3])+ 
  geom_vline(xintercept = 83.5, linetype = "dotted", col = clr_4[4])+ 
  scale_color_manual(values=clr_4)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") -> across_genes

tss <- subset(cpg_windows, region == "TSS")
tss$methperc_n <- as.numeric(tss$methperc_n)

ggplot(tss, aes(x = methperc_mean)) + geom_histogram() +
  labs(x = "Mean CpG methylation %", title = "Histogram TSS CpG methylation %", y = "Count") -> tss_hist

ggplot(tss, aes(x = methperc_n, y = methperc_mean)) + geom_point() +
  labs(x = "Number of CpGs", title = "Number of CpG sites vs mean methylation", y = "Mean CpG methylation %") +
  ylim(0, 100)+
  geom_smooth(method="lm")-> tss_n_vs_meth

cowplot::plot_grid(tss_hist, tss_n_vs_meth, ncol = 2, align = "hv", axis = "bl") -> tss_a
cowplot::plot_grid(tss_a, across_genes, ncol = 1, align = "hv", axis = "bl") -> tss_sum

ggsave(tss_sum, file = "plots/summary_wgbs_TSS.png", width=12, height=12)


