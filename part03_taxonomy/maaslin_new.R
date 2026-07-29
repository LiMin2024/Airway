setwd("/share/data1/limin/Airway/02.taxo/")
getwd()
rm(list=ls())

library(Maaslin2)


myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

###############
load("taxonomy.known.RData")

dt <- dt_norm
dt2 <- rownames_to_column(dt,var="name")
dt2$id <- paste0("V", rownames(dt2))
dt3 <- dt2 %>%
  select(name,id)
#write.table(dt3,file = "species_id.txt", sep = "\t", row.names = T, quote = FALSE)

dt <- dt2 %>% select(-c(name, id))
rownames(dt) <- paste0("V", rownames(dt))

sample_map = read.table("../00.data/sample.info.txt", header = T,  sep = "\t")
names(sample_map)[1] <- "Sample"
names(sample_map)[3] <- "Site"
names(sample_map)[5] <- "Country"
rownames(sample_map) <- sample_map$Sample

sample_map <- sample_map[sample_map$Country %in% c("China", "USA", "Germany"), ]

result <- align_dt_sample(dt, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

sample_map <- sample_map[sample_map$Site == "Throat",]#"Sputum", "BALF", "Nasopharynx", "Oropharynx", 
#                                                                 "Throat", "Nose"
names(sample_map)[1] <- "Sample"

intersect_id = intersect(colnames(aligned_dt), sample_map$Sample)
sampf = unique(subset(sample_map, Sample %in% intersect_id, c("Country","Sample")))
rownames(sampf) = sampf$Sample
dtf = t(aligned_dt[,sampf$Sample])
dtf <- as.data.frame(dtf)
rownames(dtf) = rownames(sampf)

#unique(sampf$Site) "Sputum"   "BALF" "Oropharynx" "Nasopharynx"    "Throat" "Nose"
#unique(sampf$Country)  "China" "USA"  "Germany"    
#sampf2 <- sampf %>%
#  mutate(Site = ifelse(Site == "Sputum", Site, "Others"))

sampf2 <- sampf %>%
  mutate(Country = ifelse(Country == "USA", Country, "Others"))

library(Maaslin2)
output_dir <- "Site/Throat/USA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fit_data <- Maaslin2(
  input_data = dtf,
  input_metadata = sampf2,
  output = output_dir,
  min_prevalence = 0.1,
  normalization = "NONE",
  fixed_effects = c("Country"),       
  random_effects = c(),           
  reference = c("Country,Others"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

##############Country
China <- read.table("Site/Throat/China/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)
USA <- read.table("Site/Throat/USA/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)
Germany <- read.table("Site/Throat/Germany/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)

all_long <- bind_rows(
  China %>% select(feature, coef, qval) %>% mutate(type = "China"),
  USA   %>% select(feature, coef, qval) %>% mutate(type = "USA"),
  Germany %>% select(feature, coef, qval) %>% mutate(type = "Germany"))

all_long <- bind_rows(
  China %>% select(feature, coef, qval) %>% mutate(type = "China"),
  USA   %>% select(feature, coef, qval) %>% mutate(type = "USA"))

all_long <- bind_rows(
  China %>% select(feature, coef, qval) %>% mutate(type = "China"),
  Germany %>% select(feature, coef, qval) %>% mutate(type = "Germany"))

all_long <- all_long[all_long$qval < 0.05,]
#write.table(all_long,file = "Nasopharynx_maaslin.result", sep = "\t", row.names = F, quote = FALSE)

##################################
load("taxonomy.known.RData")

dt <- dt_norm
dt2 <- rownames_to_column(dt,var="name")
dt2$id <- paste0("V", rownames(dt2))
dt3 <- dt2 %>%
  select(name,id)
#write.table(dt3,file = "species_id.txt", sep = "\t", row.names = T, quote = FALSE)

dt <- dt2 %>% select(-c(name, id))
rownames(dt) <- paste0("V", rownames(dt))

sample_map = read.table("../00.data/sample.info.txt", header = T,  sep = "\t")
names(sample_map)[1] <- "Sample"
names(sample_map)[3] <- "Site"
names(sample_map)[5] <- "Country"
rownames(sample_map) <- sample_map$Sample

sample_map <- sample_map[sample_map$Site %in% c("Sputum", "BALF", "Nasopharynx", "Oropharynx", "Throat", "Nose"), ]

result <- align_dt_sample(dt, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

sample_map <- sample_map[sample_map$Country == "Germany",]
names(sample_map)[1] <- "Sample"

intersect_id = intersect(colnames(aligned_dt), sample_map$Sample)
sampf = unique(subset(sample_map, Sample %in% intersect_id, c("Site","Sample")))
rownames(sampf) = sampf$Sample
dtf = t(aligned_dt[,sampf$Sample])
dtf <- as.data.frame(dtf)
rownames(dtf) = rownames(sampf)

################
#unique(sampf$Site) "Sputum"   "BALF" "Oropharynx" "Nasopharynx"    "Throat" "Nose"
#unique(sampf$Country)  "China" "USA"  "Germany"    
sampf2 <- sampf %>%
  mutate(Site = ifelse(Site == "Nose", Site, "Others"))

#sampf2 <- sampf %>%
#  mutate(Country = ifelse(Country == "Switzerland", Country, "Others"))

library(Maaslin2)
output_dir <- "Nose"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 运行 MaAsLin2，Day 为固定效应，PRJ 作为协变量
fit_data <- Maaslin2(
  input_data = dtf,
  input_metadata = sampf2,
  output = output_dir,
  min_prevalence = 0.1,
  normalization = "NONE",
  fixed_effects = c("Site"),       
  random_effects = c(),           
  reference = c("Site,Others"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

##############
# 获取所有的采样部位
sites <- unique(sampf$Site)

# 遍历每一个部位
for (current_site in sites) {
  
  # 1. 动态修改元数据：将当前部位保留，其余全部归为 "Others"
  sampf2 <- sampf %>%
    mutate(Site = ifelse(Site == current_site, current_site, "Others"))
  
  # 2. 为每个部位创建独立的输出目录
  output_dir <- paste0("Maaslin2_", current_site)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 3. 打印进度，方便在控制台查看运行状态
  message(paste("正在分析部位:", current_site))
  
  # 4. 运行 MaAsLin2
  fit_data <- Maaslin2(
    input_data = dtf,
    input_metadata = sampf2,
    output = output_dir,
    min_prevalence = 0.1,
    normalization = "NONE",
    fixed_effects = c("Site"),       
    random_effects = c(),           
    reference = c("Site", "Others"), # 注意：这里建议用 c("Site", "Others") 而不是逗号连接的字符串
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
}

################test
Sputum <- read.table("Germany/Maaslin2_Sputum/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)
BALF <- read.table("Germany/Maaslin2_BALF/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)
Nasopharynx <- read.table("Germany/Maaslin2_Nasopharynx/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)
Oropharynx <- read.table("Germany/Maaslin2_Oropharynx/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)
Throat <- read.table("Germany/Maaslin2_Throat/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)
Nose <- read.table("Germany/Maaslin2_Nose/all_results.tsv", header = T,  sep = "\t",check.names = FALSE)

all_long <- bind_rows(
  Sputum %>% select(feature, coef, qval) %>% mutate(type = "Sputum"),
  BALF   %>% select(feature, coef, qval) %>% mutate(type = "BALF"),
  Nasopharynx %>% select(feature, coef, qval) %>% mutate(type = "Nasopharynx"),
  Oropharynx   %>% select(feature, coef, qval) %>% mutate(type = "Oropharynx"),
  Throat %>% select(feature, coef, qval) %>% mutate(type = "Throat"),
  Nose %>% select(feature, coef, qval) %>% mutate(type = "Nose"))

all_long <- bind_rows(
  Sputum %>% select(feature, coef, qval) %>% mutate(type = "Sputum"),
  BALF   %>% select(feature, coef, qval) %>% mutate(type = "BALF"),
  Nasopharynx %>% select(feature, coef, qval) %>% mutate(type = "Nasopharynx"),
  Throat %>% select(feature, coef, qval) %>% mutate(type = "Throat"),
  Nose %>% select(feature, coef, qval) %>% mutate(type = "Nose"))

all_long <- bind_rows(
  Sputum %>% select(feature, coef, qval) %>% mutate(type = "Sputum"),
  Nasopharynx %>% select(feature, coef, qval) %>% mutate(type = "Nasopharynx"),
  Oropharynx   %>% select(feature, coef, qval) %>% mutate(type = "Oropharynx"),
  Throat %>% select(feature, coef, qval) %>% mutate(type = "Throat"))
all_long <- all_long[all_long$qval < 0.05,]
#write.table(all_long,file = "Germany_maaslin.result", sep = "\t", row.names = F, quote = FALSE)

Sputum <- read.table("Sputum_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
Sputum$Type <-  "Sputum"
BALF <- read.table("BALF_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
BALF$Type <-  "BALF"
Nasopharynx <- read.table("Nasopharynx_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
Nasopharynx$Type <-  "Nasopharynx"
Oropharynx <- read.table("Oropharynx_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
Oropharynx$Type <-  "Oropharynx"
Throat <- read.table("Throat_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
Throat$Type <-  "Throat"
Nose <- read.table("Nose_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
Nose$Type <-  "Nose"
China <- read.table("China_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
China$Type <-  "China"
USA <- read.table("USA_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
USA$Type <-  "USA"
Germany <- read.table("Germany_maaslin.result", header = T,  sep = "\t",check.names = FALSE)
Germany$Type <-  "Germany"
final_data1 <- rbind(Sputum,BALF,Nasopharynx,Oropharynx,Throat,Nose)
final_data2 <- rbind(China,USA,Germany)
final_data3 <- rbind(Sputum,Nasopharynx,Throat)

taxo <- read.table("merged_data", header = F,  sep = "\t",check.names = FALSE)
nrow(taxo[taxo$V3 == "Bacteria",])
nrow(taxo[taxo$V3 == "Eukaryota",])
nrow(taxo[taxo$V3 == "Virus",])

length(unique(final_data1$feature))
length(unique(final_data2$feature))
length(unique(Sputum$feature))
length(unique(BALF$feature))
length(unique(Nasopharynx$feature))
length(unique(Oropharynx$feature))
length(unique(Throat$feature))
length(unique(Nose$feature))
length(unique(China$feature))
length(unique(USA$feature))
length(unique(Germany$feature))

##############
load("xx/Sputum/Sputum_masslin2.RData")
Sputum_xx <- res[res$qval < 0.05,]
load("xx/BALF/BALF_masslin2.RData")
BALF_xx <- res[res$qval < 0.05,]
load("xx/BALF/BALF_masslin2.RData")
BALF_xx <- res[res$qval < 0.05,]
load("xx/Nasopharynx/Nasopharynx_masslin2.RData")
Nasopharynx_xx <- res[res$qval < 0.05,]
load("xx/Oropharynx/Oropharynx_masslin2.RData")
Oropharynx_xx <- res[res$qval < 0.05,]
load("xx/Throat/Throat_masslin2.RData")
Throat_xx <- res[res$qval < 0.05,]
load("xx/Nose/Nose_masslin2.RData")
Nose_xx <- res[res$qval < 0.05,]

final_dataxx <- rbind(Sputum_xx,BALF_xx,Nasopharynx_xx,Oropharynx_xx,Throat_xx,Nose_xx)
final_dataxx <- final_dataxx[!grepl("\\.", final_dataxx$feature), ]
length(unique(final_dataxx$feature))

####################
dt <- aligned_dt
dt2 <- rownames_to_column(dt,var="name")
high_abundance_features <- dt2 %>%
  # 将除了第1列（物种名）之外的所有列转换为数值型（防止文本干扰）
  mutate(across(-1, as.numeric)) %>%
  # 计算行平均值
  mutate(Mean_Abundance = rowMeans(across(-1), na.rm = TRUE)) %>%
  # 只保留物种名和平均丰度
  select(1, Mean_Abundance)
high_abundance_features2 <- high_abundance_features[high_abundance_features$Mean_Abundance > 0.05,]

############
sig.features <- final_dataxx
sig.features2 <- sig.features[sig.features$feature %in% high_abundance_features2$name,]
#write.table(sig.features2,file = "Site_maaslin.343", sep = "\t", row.names = T, quote = FALSE)

sig.features3 <- sig.features2[!duplicated(sig.features2$feature),]

id <- read.table("species_id.txt", sep="\t", header=T, check.names=F)

data <- left_join(sig.features3,id,by=c("feature"="id"))

tax <- read.table("geneset.taxonomy.select", sep="\t", header=F, check.names=F)
tax <- tax[,-c(1,4)]
names(tax) <- c("id","domain","phylum","class","order","family","genus","species")

merge_dt <- left_join(data,tax,by=c("name"="species"))
#write.table(merge_dt,file = "Site_maaslin.343_tax.txt", sep = "\t", row.names = T, quote = FALSE)

#################
library(tidyr)
library(dplyr)
dt_long <- dt2 %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Abundance")

dt_grouped <- dt_long %>%
  left_join(aligned_sample, by = "Sample")

mean_abundance <- dt_grouped %>%
  group_by(name, Site) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

mean_abundance_wide <- mean_abundance %>%
  pivot_wider(names_from = Site, values_from = MeanAbundance) %>%
  replace(is.na(.), 0) 
mean_abundance_wide2 <- mean_abundance_wide[mean_abundance_wide$name %in% merge_dt$feature,]
#write.table(mean_abundance_wide2,file = "Site.343_mean_abun.txt", sep = "\t", row.names = T, quote = FALSE)

##############
cols <- c('Sputum', 'BALF', 'Nasopharynx', 'Oropharynx', 'Throat', 'Nose')

result2 <- mean_abundance_wide2 %>%
  rowwise() %>%
  mutate(
    high = {
      vals <- c_across(all_of(cols))
      if (all(is.na(vals))) NA_character_ else cols[which.max(vals)]
    },
    low = {
      vals <- c_across(all_of(cols))
      if (all(is.na(vals))) NA_character_ else cols[which.min(vals)]
    }
  ) %>%
  ungroup()

result3 <- result2[result2$name %in% sig.features3$feature,]
#write.table(result3,file = "Site.343_high_low.txt", sep = "\t", row.names = T, quote = FALSE)

############high
Sputum <- result3[result3$high == "Sputum",]
BALF <- result3[result3$high == "BALF",]
Nasopharynx <- result3[result3$high == "Nasopharynx",]
Oropharynx <- result3[result3$high == "Oropharynx",]
Throat <- result3[result3$high == "Throat",]
Nose <- result3[result3$high == "Nose",]

############low
Sputum <- result3[result3$low == "Sputum",]
BALF <- result3[result3$low == "BALF",]
Nasopharynx <- result3[result3$low == "Nasopharynx",]
Oropharynx <- result3[result3$low == "Oropharynx",]
Throat <- result3[result3$low == "Throat",]
Nose <- result3[result3$low == "Nose",]

aligned_dt2 <- rownames_to_column(aligned_dt,var="name")
aligned_dt3 <- aligned_dt2[aligned_dt2$name %in% Sputum$name,]#"Sputum"   "BALF" "Oropharynx" "Nasopharynx"    "Throat" "Nose"
rownames(aligned_dt3) <- NULL
aligned_dt3 <- column_to_rownames(aligned_dt3,var="name")

result <- align_dt_sample(aligned_dt3, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

intersect_id = intersect(colnames(aligned_dt), sample_map$Sample)
sampf = unique(subset(sample_map, Sample %in% intersect_id, c("Site","Sample")))
rownames(sampf) = sampf$Sample
dtf = t(aligned_dt[,sampf$Sample])
dtf <- as.data.frame(dtf)
rownames(dtf) = rownames(sampf)

#unique(sampf$SampleType) "Sputum"  "Oropharynx" "Nasopharynx"   "BALF" "Anterior nares"
#unique(sampf$Country)  "America"  "China" "European.countries" 
sampf2 <- sampf %>%
  mutate(Site = ifelse(Site == "Sputum", Site, "Others"))

library(Maaslin2)
output_dir <- "Abun/Sputum"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 运行 MaAsLin2，Day 为固定效应，PRJ 作为协变量
fit_data <- Maaslin2(
  input_data = dtf,
  input_metadata = sampf2,
  output = output_dir,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = c("Site"),       
  random_effects = c(),           
  reference = c("Site,Others"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

Sputum_high <- read.table("Abun/Sputum/all_results.tsv", sep="\t", header=T, check.names=F)
Sputum_high$high <- "Sputum"

BALF_high <- read.table("BALF/all_results.tsv", sep="\t", header=T, check.names=F)
BALF_high$high <- "BALF"
Nasopharynx_high <- read.table("Nasopharynx/all_results.tsv", sep="\t", header=T, check.names=F)
Nasopharynx_high$high <- "Nasopharynx"
Oropharynx_high <- read.table("Oropharynx/all_results.tsv", sep="\t", header=T, check.names=F)
Oropharynx_high$high <- "Oropharynx"
Sputum_high <- read.table("Sputum/all_results.tsv", sep="\t", header=T, check.names=F)
Sputum_high$high <- "Sputum"

high_data <- rbind(Anterior_high,BALF_high,Nasopharynx_high,Oropharynx_high,Sputum_high)
high_data$enriched <- ifelse(high_data$qval < 0.05, high_data$high, NA_character_)

Anterior_low <- read.table("Anterior2/all_results.tsv", sep="\t", header=T, check.names=F)
Anterior_low$low <- "Anterior nares"
BALF_low <- read.table("BALF2/all_results.tsv", sep="\t", header=T, check.names=F)
BALF_low$low <- "BALF"
Nasopharynx_low <- read.table("Nasopharynx2/all_results.tsv", sep="\t", header=T, check.names=F)
Nasopharynx_low$low <- "Nasopharynx"
Oropharynx_low <- read.table("Oropharynx2/all_results.tsv", sep="\t", header=T, check.names=F)
Oropharynx_low$low <- "Oropharynx"
Sputum_low <- read.table("Sputum2/all_results.tsv", sep="\t", header=T, check.names=F)
Sputum_low$low <- "Sputum"

low_data <- rbind(Anterior_low,BALF_low,Nasopharynx_low,Oropharynx_low,Sputum_low)
low_data$deficient <- ifelse(low_data$qval < 0.05, low_data$low, NA_character_)

merge_dt <- merge(high_data[,c(1,ncol(high_data))],low_data[,c(1,ncol(low_data))],by="feature")
#write.table(merge_dt,file = "SampleType_enriched_deficient.txt", sep = "\t", row.names = F, quote = FALSE)
