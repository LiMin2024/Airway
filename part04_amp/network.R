setwd("/share/data1/limin/Airway/03.AMP/")
taxo <- myread("../02.taxo/all.taxonomy.profile.0.0005")

sample_map = read.table("../00.data/sample.info.txt", header = T,  sep = "\t")
names(sample_map)[1] <- "Sample"
sample_map <- sample_map%>%
  select(Sample,site,country)
names(sample_map) <- c("Sample","Site","Country")
rownames(sample_map) <- sample_map$Sample

sample_map <- sample_map[sample_map$Site %in% c("Sputum", "Nasopharynx",  "Throat"), ]
sample_map <- sample_map[sample_map$Country %in% c("China", "USA", "Germany"), ]

result <- align_dt_sample(taxo, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

taxo_dt <- aligned_dt

#############
data <- myread("profs.amp")
amp_profile <- data
amp <- rownames_to_column(data,var="gene")

result <- align_dt_sample(amp_profile, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map
amp_t <- t(aligned_dt)

prevalence <- colMeans(amp_t > 0, na.rm = TRUE)
prevalence_df <- data.frame(
  AMP = names(prevalence),
  Prevalence = prevalence
)
prevalence_df2 <- prevalence_df[prevalence_df$Prevalence > 0.01,]
#save(prevalence_df2,file="amp_profile_0.1.RData")

amp_dt <- aligned_dt[rownames(aligned_dt)%in%prevalence_df2$AMP, ]

sample_map <- aligned_sample

sample_map1 <- sample_map[sample_map$Site == "Sputum",]
sample_map2 <- sample_map[sample_map$Site == "Nasopharynx",]
sample_map3 <- sample_map[sample_map$Site == "Throat",]

sample_map4 <- sample_map[sample_map$Country == "China",]
sample_map5 <- sample_map[sample_map$Country == "USA",]
sample_map6 <- sample_map[sample_map$Country == "Germany",]

result1 <- align_dt_sample(amp_dt, sample_map1, ID = "Sample")
aligned_dt1 <- result1$dt
aligned_sample1 <- result1$sample_map

result2 <- align_dt_sample(taxo_dt, sample_map1, ID = "Sample")
aligned_dt2 <- result2$dt
aligned_sample2 <- result2$sample_map
aligned_dt2 <- aligned_dt2[, colnames(aligned_dt1) ]  # 对齐样本顺序
write.table(aligned_dt1,file = "amp_Sputum", sep = "\t", row.names = T, quote = FALSE)
write.table(aligned_dt2,file = "taxo_Sputum", sep = "\t", row.names = T, quote = FALSE)

result1 <- align_dt_sample(amp_dt, sample_map2, ID = "Sample")
aligned_dt1 <- result1$dt
aligned_sample1 <- result1$sample_map

result2 <- align_dt_sample(taxo_dt, sample_map2, ID = "Sample")
aligned_dt2 <- result2$dt
aligned_sample2 <- result2$sample_map
aligned_dt2 <- aligned_dt2[, colnames(aligned_dt1) ]  # 对齐样本顺序
write.table(aligned_dt1,file = "amp_Nasopharynx", sep = "\t", row.names = T, quote = FALSE)
write.table(aligned_dt2,file = "taxo_Nasopharynx", sep = "\t", row.names = T, quote = FALSE)

result1 <- align_dt_sample(amp_dt, sample_map3, ID = "Sample")
aligned_dt1 <- result1$dt
aligned_sample1 <- result1$sample_map

result2 <- align_dt_sample(taxo_dt, sample_map3, ID = "Sample")
aligned_dt2 <- result2$dt
aligned_sample2 <- result2$sample_map
aligned_dt2 <- aligned_dt2[, colnames(aligned_dt1) ]  # 对齐样本顺序
write.table(aligned_dt1,file = "amp_Throat", sep = "\t", row.names = T, quote = FALSE)
write.table(aligned_dt2,file = "taxo_Throat", sep = "\t", row.names = T, quote = FALSE)

result1 <- align_dt_sample(amp_dt, sample_map4, ID = "Sample")
aligned_dt1 <- result1$dt
aligned_sample1 <- result1$sample_map

result2 <- align_dt_sample(taxo_dt, sample_map4, ID = "Sample")
aligned_dt2 <- result2$dt
aligned_sample2 <- result2$sample_map
aligned_dt2 <- aligned_dt2[, colnames(aligned_dt1) ]  # 对齐样本顺序
write.table(aligned_dt1,file = "amp_China", sep = "\t", row.names = T, quote = FALSE)
write.table(aligned_dt2,file = "taxo_China", sep = "\t", row.names = T, quote = FALSE)

result1 <- align_dt_sample(amp_dt, sample_map5, ID = "Sample")
aligned_dt1 <- result1$dt
aligned_sample1 <- result1$sample_map

result2 <- align_dt_sample(taxo_dt, sample_map5, ID = "Sample")
aligned_dt2 <- result2$dt
aligned_sample2 <- result2$sample_map
aligned_dt2 <- aligned_dt2[, colnames(aligned_dt1) ]  # 对齐样本顺序
write.table(aligned_dt1,file = "amp_USA", sep = "\t", row.names = T, quote = FALSE)
write.table(aligned_dt2,file = "taxo_USA", sep = "\t", row.names = T, quote = FALSE)

result1 <- align_dt_sample(amp_dt, sample_map6, ID = "Sample")
aligned_dt1 <- result1$dt
aligned_sample1 <- result1$sample_map

result2 <- align_dt_sample(taxo_dt, sample_map6, ID = "Sample")
aligned_dt2 <- result2$dt
aligned_sample2 <- result2$sample_map
aligned_dt2 <- aligned_dt2[, colnames(aligned_dt1) ]  # 对齐样本顺序
write.table(aligned_dt1,file = "amp_Germany", sep = "\t", row.names = T, quote = FALSE)
write.table(aligned_dt2,file = "taxo_Germany", sep = "\t", row.names = T, quote = FALSE)

#############用py脚本计算spearman
setwd("/share/data1/limin/Airway/03.AMP/network/")
rcorr_result = read.table("Throat_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

China_result <- threshold
China_result$Site <- "China"
China_result <- China_result %>% 
  mutate(group = paste(A, B, sep = "|"))
China_result <- China_result[China_result$r.label == "negative",]

rcorr_result = read.table("USA_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

USA_result <- threshold
USA_result$Site <- "USA"
USA_result <- USA_result %>% 
  mutate(group = paste(A, B, sep = "|"))
USA_result <- USA_result[USA_result$r.label == "negative",]

rcorr_result = read.table("Germany_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

Germany_result <- threshold
Germany_result$Site <- "Germany"
Germany_result <- Germany_result %>% 
  mutate(group = paste(A, B, sep = "|"))
Germany_result <- Germany_result[Germany_result$r.label == "negative",]

rcorr_result = read.table("Sputum_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

Sputum_result <- threshold
Sputum_result$Site <- "Sputum"
Sputum_result <- Sputum_result %>% 
  mutate(group = paste(A, B, sep = "|"))
Sputum_result <- Sputum_result[Sputum_result$r.label == "negative",]

rcorr_result = read.table("Nasopharynx_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

Nasopharynx_result <- threshold
Nasopharynx_result$Site <- "Nasopharynx"
Nasopharynx_result <- Nasopharynx_result %>% 
  mutate(group = paste(A, B, sep = "|"))
Nasopharynx_result <- Nasopharynx_result[Nasopharynx_result$r.label == "negative",]

rcorr_result = read.table("Throat_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

Throat_result <- threshold
Throat_result$Site <- "Throat"
Throat_result <- Throat_result %>% 
  mutate(group = paste(A, B, sep = "|"))
Throat_result <- Throat_result[Throat_result$r.label == "negative",]

library(purrr) 
library(tidyverse) 
setwd("/share/data1/limin/Airway/03.AMP/")
combind_dt1 <- list(China_result[,9:8], USA_result[,9:8],Germany_result[,9:8]) %>%
  reduce(full_join, by = "group")
names(combind_dt1) <- c("group","China","USA","Germany")
#write.table(combind_dt1,file = "Country_negative_amp", sep = "\t", row.names = F, quote = FALSE)

combind_dt2 <- list(Sputum_result[,9:8], Nasopharynx_result[,9:8],Throat_result[,9:8]) %>%
  reduce(full_join, by = "group")
names(combind_dt2) <- c("group","Sputum","Nasopharynx","Throat")
#write.table(combind_dt2,file = "Site_negative_amp", sep = "\t", row.names = F, quote = FALSE)

fun_taxo_no <- taxo2[taxo2$V2 == "Eukaryota"&taxo2$V3 != "k_Fungi",]
fun_taxo_no <- as.data.frame(unique(fun_taxo_no$V9))
names(fun_taxo_no) <- "V9"

combind_dt1$taxo_name <- str_extract(combind_dt1$group, "(?<=\\|).*")
combind_dt1_clean <- combind_dt1[!(combind_dt1$taxo_name %in% fun_taxo_no$V9), ]
combind_dt1_clean <- combind_dt1_clean[,-5]

combind_dt2$taxo_name <- str_extract(combind_dt2$group, "(?<=\\|).*")
combind_dt2_clean <- combind_dt2[!(combind_dt2$taxo_name %in% fun_taxo_no$V9), ]
combind_dt2_clean <- combind_dt2_clean[,-5]

setwd("/share/data1/limin/Airway/03.AMP/")
combind_dt1 = read.table("network/Country_negative_amp", header = T,  sep = "\t")
combind_dt1 <- combind_dt1 %>% 
  filter(!grepl("s_Homo sapiens", group))
combind_dt2 = read.table("network/Site_negative_amp", header = T,  sep = "\t")
combind_dt2 <- combind_dt2 %>% 
  filter(!grepl("s_Homo sapiens", group))
################
#cols <- c( "China", "USA", "Germany")
#cols <- c( "Sputum", "Nasopharynx", "Throat")

result1 <- combind_dt1 %>% 
  separate(col = group, into = c("A", "B"), sep = "\\|", remove = TRUE)

result2 <- combind_dt2 %>% 
  separate(col = group, into = c("A", "B"), sep = "\\|", remove = TRUE)

###############
result_filtered <- result2 %>%
  filter(
    !grepl("[0-9]", B),            
    !grepl("uncultured", B, ignore.case = TRUE),
    !grepl("sp\\.", B, ignore.case = TRUE) 
  )
result_filtered <- result2

df_long <- result_filtered %>%
  pivot_longer(
    cols = starts_with(c( "Sputum","Nasopharynx","Throat")), # "Sputum","Nasopharynx","Throat" "China", "USA", "Germany"
    names_to = "sample_type",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>% # 移除没有样本类型的记录
  select(-value)            # 删除不必要的列

color_mapping <- c(China= "#e31a1c", USA = "#fdbf6f",Germany = "#cab2d6",Sputum = "#a6cee3", Nasopharynx = "#b2df8a", Throat = "#fb9a99")

df_summary <- df_long %>%
  group_by(A, B) %>%
  summarise(
    types = list(sample_type),
    .groups = 'drop'
  ) %>%
  mutate(
    n_types = lengths(types),
    first_type = sapply(types, `[`, 1),
    
    # 使用 case_when 处理多种颜色类型的映射
    fill_color = case_when(
      # 条件1：只有一种类型时，从 color_mapping 中提取对应颜色
      n_types == 1 ~ color_mapping[first_type],
      
      # 条件2：包含 China 和 Germany 的组合，指定第一种颜色
      n_types == 2 & all(c("China", "Germany") %in% unlist(types)) ~ "#1f78b4", 
      
      # 条件3：包含 USA 和 Germany 的组合，指定第二种颜色
      n_types == 2 & all(c("USA", "Germany") %in% unlist(types)) ~ "#33a02c",
      
      # 兜底条件：其他所有未匹配的情况，使用标准灰色
      TRUE ~ "#999999" 
    ),
    
    # 显示文本逻辑保持不变
    display_text = ifelse(n_types == 1, first_type, "Multiple")
  )

ggplot(df_summary, aes(x = A, y = B)) +
  geom_tile(aes(fill = fill_color), color = "white") +
  scale_fill_identity() +
  guides(fill = "legend") +  # 强制显示图例
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),   # 去掉X轴刻度标签
    axis.ticks.x = element_blank(),  # 去掉X轴刻度线
    axis.text.y = element_text(size = 12),
    aspect.ratio = 2,
    
    # 添加绘图区域边框
    panel.border = element_rect(
      fill = NA,          # 必须设置为透明，否则会遮挡图表内容
      color = "black",    # 边框颜色
      size = 0.5,           # 边框粗细
      linetype = "solid"  # 边框线型
    )
  ) +
  labs(x = "AMP", y = "")

ggplot(df_summary, aes(x = A, y = B)) +
  geom_tile(aes(fill = fill_color), color = "white") +
  scale_fill_identity() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    legend.position = "none"
  ) +
  labs(x = "AMP", y = "Species")
