setwd("/share/data1/limin/airway_geneset/Analysis/rgi/01.com/")
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(ggpubr)
library(tibble)
library(dplyr) 

align_dt_sample <- function(dt, sample_map, ID=NA){
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dt = dt[,sample_map[,ID]] %>% filter(rowSums(.) !=0)
  list(dt=dt, sample_map=sample_map)
}

load("vfdb.known.RData")
vfdb_dt <- dt_norm
load("rgi.known.RData")
rgi_dt <- dt_norm

sample_map = read.table("sample.info.txt", header = T,  sep = "\t")
names(sample_map)[1] <- "Sample"
sample_map <- sample_map%>%
  select(Sample,site,country)
names(sample_map) <- c("Sample","Site","Country")
rownames(sample_map) <- sample_map$Sample

data_map <- sample_map[sample_map$Country == "China",]#"China", "USA", "Germany"
#"Sputum", "BALF", "Nasopharynx", "Oropharynx", "Throat", "Nose"

result <- align_dt_sample(vfdb_dt, data_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

data <- as.data.frame(t(aligned_dt))
data <- rownames_to_column(data,var="Sample")
merge_data <- merge(data,aligned_sample[,c(1,2)],by="Sample")

merge_data[,-c(1,ncol(merge_data))] <- (merge_data[,-1] > 0) * 1

merge_data1 <- merge_data[merge_data$Site == "Sputum",]
merge_data1 <- merge_data1[,-ncol(merge_data1)]
merge_data2 <- merge_data[merge_data$Site == "BALF",]
merge_data2 <- merge_data2[,-ncol(merge_data2)]
merge_data3 <- merge_data[merge_data$Site == "Nasopharynx",]
merge_data3 <- merge_data3[,-ncol(merge_data3)]
merge_data4 <- merge_data[merge_data$Site == "Oropharynx",]
merge_data4 <- merge_data4[,-ncol(merge_data4)]
merge_data5 <- merge_data[merge_data$Site == "Throat",]
merge_data5 <- merge_data5[,-ncol(merge_data5)]
merge_data6 <- merge_data[merge_data$Site == "Nose",]
merge_data6 <- merge_data6[,-ncol(merge_data6)]

get_expressed_genes <- function(data, threshold = 0) {
  # 提取基因列（去掉 Sample 列）
  gene_cols <- data[, -1]  # 假设第一列是 Sample
  # 找出在任一样本中丰度 > 阈值的基因
  expressed <- apply(gene_cols, 2, function(x) any(x > threshold, na.rm = TRUE))
  return(names(gene_cols)[expressed])
}

genes_list <- list(
  "Sputum" = get_expressed_genes(merge_data1),
  "BALF"           = get_expressed_genes(merge_data2),
  "Nasopharynx"    = get_expressed_genes(merge_data3),
  "Oropharynx"     = get_expressed_genes(merge_data4),
  "Throat"         = get_expressed_genes(merge_data5),
  "Nose"         = get_expressed_genes(merge_data6)
)
lapply(genes_list, length)
all_genes <- unique(unlist(genes_list))
#save(genes_list,file="China.upset.RData")

# install.packages("UpSetR") 
library(UpSetR)
library(gridExtra)
upset(fromList(genes_list), 
      nsets = 5,           # 显示全部5个集合
      nintersects = 30,    # 显示前30个最显著的交集（可根据需要调整或设为NA显示全部）
      order.by = "freq",   # 按照交集的频率（元素数量）降序排列
      mb.ratio = c(0.5, 0.5),
      sets.bar.color = "#56B4E9", # 左侧集合大小条形图的颜色
      point.size = 3,      # 矩阵中点的大小
      line.size = 1,       # 连线的粗细
      text.scale = 1.5     # 整体文字标签的缩放比例
)

###############
data_map <- sample_map[sample_map$Site == "Nose",]#"China", "USA", "Germany"
#"Sputum", "BALF", "Nasopharynx", "Oropharynx", "Throat", "Nose"

result <- align_dt_sample(vfdb_dt, data_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

data <- as.data.frame(t(aligned_dt))
data <- rownames_to_column(data,var="Sample")
merge_data <- merge(data,aligned_sample[,c(1,3)],by="Sample")

merge_data[,-c(1,ncol(merge_data))] <- (merge_data[,-1] > 0) * 1

merge_data1 <- merge_data[merge_data$Country == "China",]
merge_data1 <- merge_data1[,-ncol(merge_data1)]
merge_data2 <- merge_data[merge_data$Country == "USA",]
merge_data2 <- merge_data2[,-ncol(merge_data2)]
merge_data3 <- merge_data[merge_data$Country == "Germany",]
merge_data3 <- merge_data3[,-ncol(merge_data3)]

get_expressed_genes <- function(data, threshold = 0) {
  # 提取基因列（去掉 Sample 列）
  gene_cols <- data[, -1]  # 假设第一列是 Sample
  # 找出在任一样本中丰度 > 阈值的基因
  expressed <- apply(gene_cols, 2, function(x) any(x > threshold, na.rm = TRUE))
  return(names(gene_cols)[expressed])
}

genes_list <- list(
  "China" = get_expressed_genes(merge_data1),
  "USA"           = get_expressed_genes(merge_data2),
  "Germany"    = get_expressed_genes(merge_data3)
)

lapply(genes_list, length)
all_genes <- unique(unlist(genes_list))
#save(genes_list,file="Nose.upset.RData")

# install.packages("UpSetR") 
library(UpSetR)
library(gridExtra)
upset(fromList(genes_list), 
      nsets = 5,           # 显示全部5个集合
      nintersects = 30,    # 显示前30个最显著的交集（可根据需要调整或设为NA显示全部）
      order.by = "freq",   # 按照交集的频率（元素数量）降序排列
      mb.ratio = c(0.5, 0.5),
      sets.bar.color = "#56B4E9", # 左侧集合大小条形图的颜色
      point.size = 3,      # 矩阵中点的大小
      line.size = 1,       # 连线的粗细
      text.scale = 1.5     # 整体文字标签的缩放比例
)
