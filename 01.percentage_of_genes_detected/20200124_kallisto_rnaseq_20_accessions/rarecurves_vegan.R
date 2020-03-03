# library vegan
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))


######### Rare curves ########

# counts matrix
#df = read.table("RNAcountsmatrix.txt", header = T, sep = "\t")
#dim(df)

# par(mfrow = c(2,2))
#rarecurve(df[31:33, -(1:5)], step = 100, label = F, main = "All replicates of the wt.rpi_day5 (step = 100)", ylab = "Counts", col = c("#020275", "#009E00", "#C40000"))
#legend("bottomright", c("rep1", "rep2", "rep3"), lty = 1, col = c("#11A7CC", "#E79211", "#8011CC"), y.intersp = 0.75)



# counts matrix
df = fread(input = "01.percentage_of_genes_detected/20200124_kallisto_rnaseq_20_accessions/20_accessions_mRNA_raw_count_wide.tsv", 
            header = T, 
            sep = "\t", 
            integer64 = "character")

df_filtered = filter(df, grepl(pattern = "Moneymaker*", sample)) # only moneymaher

col_sums = 
  df_filtered %>%
  select( - sample) %>% 
  colSums() 

genes_non_null = as.data.frame(col_sums) %>% mutate(target_id = row.names(.)) %>% filter(col_sums > 0)


df_filtered = dplyr::select(df_filtered, genes_non_null$target_id) %>% as.matrix()
df_filtered <- round(df_filtered, digits = 0)

par(mfrow = c(1,2))

# curves for all genes
rarecurve(df_filtered, 
          step = 100000, 
          label = F, 
          col = c("#020275", "#009E00", "#C40000","#C40000"),
          main = "Rarefaction curve", 
          
          ylab = "Counts")
legend("bottomright", c("Moneymaker-01", "Moneymaker-02", "Moneymaker-03","Moneymaker-04"), 
       lty = 1, 
       col = c("#11A7CC", "#E79211", "#8011CC","#11A7CC"), 
       y.intersp = 0.75)

# curves for TFs
# Reference: https://doi.org/10.1016/j.aggene.2017.08.002
transcription_factors <- read.delim("01.percentage_of_genes_detected/20200124_kallisto_rnaseq_20_accessions/transcription_factors.txt",header = T)
tfs <- as.vector(transcription_factors[,c("Gene.name")])
#tfs <- sapply(tfs, function(x) {substr(x,1,14)}, simplify = T) %>%  as.vector() # extract gene locus number 

df_filtered_tfs = df_filtered[, colnames(df_filtered) %in% tfs]

rarecurve(df_filtered_tfs, 
          step = 1000, 
          label = F, 
          main = "Rarefaction curve", 
          col = c("#020275", "#009E00", "#C40000","#C40000"),
          ylab = "Counts")
legend("bottomright", c("Moneymaker-01", "Moneymaker-02", "Moneymaker-03","Moneymaker-04"), 
       lty = 1, 
       col = c("#11A7CC", "#E79211", "#8011CC","#11A7CC"), 
       y.intersp = 0.75)


#########################################
# Different thresholds for gene detection
#########################################
