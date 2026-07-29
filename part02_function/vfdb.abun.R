setwd("/share/data1/limin/Airway/01.function/")
load("vfdb.known.RData")
vfdb_dt <- dt_norm
load("rgi.known.RData")
rgi_dt <- dt_norm

sample_map = read.table("../00.data/sample.info.txt", header = T,  sep = "\t")
names(sample_map)[1] <- "Sample"
sample_map <- sample_map%>%
  select(Sample,site,country)
names(sample_map) <- c("Sample","Site","Country")
rownames(sample_map) <- sample_map$Sample
sample_map <- sample_map[sample_map$Site %in% c("Sputum", "Nasopharynx",  "Throat"), ]
sample_map <- sample_map[sample_map$Country %in% c("China", "USA", "Germany"), ]

sample_map2 <- sample_map[sample_map$Site == "BALF", ]
sample_map2 <- sample_map[sample_map$Country == "China", ]

##############
result <- align_dt_sample(vfdb_dt, sample_map2, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

##############
dt2 <- aligned_dt %>%
rownames_to_column(var = "name")

map <- read.table("VFDB_setB.latest.map", header = T, fill = TRUE, sep = "\t",check.names = F)
merge <- merge(dt2,map[,c(1,6)],by.x = "name",by.y = "GeneID")
merge <- merge[,-1]
names(merge)[ncol(merge)] <- "name"

result <- merge %>%
  group_by(name) %>%   # 第一步：按照 name 列进行分组
  summarise(           # 第二步：对组内数据进行汇总计算
    across(where(is.numeric), sum, na.rm = TRUE) # 对所有数值型列求和，并自动忽略 NA 值
  )

result2 <- result %>%
  mutate(
    # 创建新列 row_mean，计算除 name 列外所有数值的行平均
    row_mean = rowMeans(select(., -name), na.rm = TRUE) 
  )

dt_long <- result %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Abundance")

# 合并 sample_map 获取每个样本对应的 Group
dt_grouped <- dt_long %>%
  left_join(aligned_sample, by = "Sample")

mean_abundance <- dt_grouped %>%
  group_by(name, Country) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

# 将结果转为宽格式，便于比较两组
mean_abundance_wide <- mean_abundance %>%
  pivot_wider(names_from = Country, values_from = MeanAbundance) %>%
  replace(is.na(.), 0)  # 缺失值补零（如果某 fungi 在某组中无样本）
#write.table(mean_abundance_wide,file = "BALF_abun_sort.txt", sep = "\t", row.names = T, quote = FALSE)

result <- mean_abundance_wide %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE))) %>%
  mutate(row_idx = row_number()) %>%
  group_by(group = ifelse(row_idx <= 10, name, "Others")) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  rename(name = group) %>%
  select(-row_idx) %>%
  select(name, everything()) %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))                      

####Sputum China
data1 <- result %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))

######Nasopharynx USA
data2 <- result %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))

#########Throat Germany
data3 <- result %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))

data <- data1 %>%
  full_join(data2, by = "name") %>%   # 第一步：将 data1 与 data2 按 name 左连接
  full_join(data3, by = "name")   

names(data) <- c("name","Sputum_China","Sputum_USA","Sputum_Germany","Nasopharynx_China","Nasopharynx_USA","Nasopharynx_Germany","Throat_China","Throat_USA","Throat_Germany")
names(data) <- c("name","China_Nasopharyn","China_Sputum","China_Throat","USA_Nasopharynx","USA_Sputum","USA_Throat","Germany_Nasopharynx","Germany_Sputum","Germany_Throat")
###########
df_long <- data %>%
  pivot_longer(cols = -name, names_to = "Country", values_to = "Abundance") %>%
  arrange(factor(name, levels = unique(result$name))) 

df_long2 <- df_long %>%
  mutate(Abundance100 = Abundance)

original_levels <- unique(data$name)
x_levels <- c(original_levels[original_levels != "Others"], "Others")
y_levels <- c("China","USA","BGermany")

y_levels <- c("Sputum_China","Sputum_USA","Sputum_Germany","Nasopharynx_China","Nasopharynx_USA","Nasopharynx_Germany","Throat_China","Throat_USA","Throat_Germany")
y_levels <- c("name","China_Nasopharyn","China_Sputum","China_Throat","USA_Nasopharynx","USA_Sputum","USA_Throat","Germany_Nasopharynx","Germany_Sputum","Germany_Throat")

library(scales)
df_plot <- df_long2 %>%
  mutate(
    Abundance_raw = Abundance100,  # 保留原始值用于标注
    Abundance_log_neg = Abundance  # 用于填色
  )

p <- ggplot(df_plot, aes(x = factor(name, levels = x_levels), y = factor(Country, levels = rev(y_levels)), fill = Abundance_log_neg)) +
  geom_tile(color = "white") +
  
  # 修改此处：使用 ifelse 判断，如果是 NA (原数据为0) 则不显示文字，否则保留两位小数
  geom_text(aes(label = ifelse(is.na(Abundance_log_neg), "", sprintf("%.2f", Abundance_raw))), 
            size = 3.5, 
            color = "black", 
            na.rm = TRUE) + 
  
  scale_fill_gradientn(
    colours = c("#0077B5","#E6F5FF","#FFE6E6","#DC143C"),
    values = rescale(range(df_plot$Abundance_log_neg, na.rm = TRUE)),
    name = "-Log(Relative Abundance)",
    na.value = "grey90"  # 0值依然填充为浅灰色
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    aspect.ratio = 0.4 
  ) +
  labs(title = "", x = "Function", y = "Country") 

print(p)
