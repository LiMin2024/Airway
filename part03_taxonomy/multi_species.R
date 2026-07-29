bac_data <- data[rownames(data) %in% bac_taxo$V9,]
fun_data <- data[rownames(data) %in% fun_taxo$V9,]
vir_data <- data[rownames(data) %in% vir_taxo$V9,]

all_data <- rbind(bac_data,fun_data,vir_data)
all_taxo  <- rbind(bac_taxo,fun_taxo,vir_taxo)

sample_map = read.table("../00.data/sample.info.txt", sep="\t", header=T, check.names=F)
names(sample_map)[1] <- "Sample"
sample_map <- sample_map%>%
  select(Sample,site,country)
names(sample_map) <- c("Sample","Site","Country")
rownames(sample_map) <- sample_map$Sample

sample_map <- sample_map[sample_map$Site %in% c("Sputum", "Nasopharynx",  "Throat"), ]
sample_map <- sample_map[sample_map$Country %in% c("China", "USA", "Germany"), ]

result <- align_dt_sample(all_data, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

norm_data <- function(dt){
  
  prof <- as.data.frame(apply(dt, 2, function(x) {
    col_sum <- sum(x)
    if (col_sum == 0) {
      return(rep(0, length(x)))  # 列和为0时返回全0
    } else {
      return(x / col_sum * 100)  # 否则正常归一化
    }
  }))
  rownames(prof) = rownames(dt)
  prof
}
dt_norm <- norm_data(aligned_dt)
#save(dt_norm,file="all.taxonomy.profile.norm")

dt <- rownames_to_column(dt_norm,var="name")
tax <- taxo2[,c(6,8)]
tax <- distinct(tax,V7, V9)

dt2 <- merge(dt,tax,by.x="name",by.y="V9")
names(dt2)[ncol(dt2)] <- "name"
dt2 <- dt2[,-1]

dt3 <- dt2 %>%
  group_by(name) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

dt3<- dt3 %>%
  mutate(name = ifelse(name == "Unclassified", "Others", name))

dt_long1 <- dt3 %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "abundance") %>%
  left_join(sample_map %>% select(Sample, Site), by = "Sample")

top10_union1 <- dt_long1 %>%
  group_by(Site, name) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site) %>%
  slice_max(mean_abund, n = 10, with_ties = FALSE) %>%
  pull(name) %>%
  unique()  # 所有组 top 10 的并集

#############
dt_long2 <- dt3 %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "abundance") %>%
  left_join(sample_map %>% select(Sample, Country), by = "Sample")

top10_union2 <- dt_long2 %>%
  group_by(Country, name) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Country) %>%
  slice_max(mean_abund, n = 10, with_ties = FALSE) %>%
  pull(name) %>%
  unique()  # 所有组 top 10 的并集

##################
source("/share/data1/limin/airway_geneset/Analysis/aa/01.composition/composition_color_map.R")
top10_union <- cbind(top10_union1,top10_union2)
total_color1 = c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                 "#6a3d9a", "#ffff99", "#b15928","#8dd3c7",
                 "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
                 "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
                 "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4",
                 "#b2df8a", "#33a02c", "#fb9a99","#ed1299","#09f9f5","#246b93","#cc8e12","#d561dd","#c93f00","#ddd53e","#4aef7b","#e86502","#9ed84e","#39ba30","#6ad157","#8249aa","#99db27","#e07233","#ff523f","#ce2523","#f7aa5d","#cebb10","#03827f","#931635","#373bbf","#a1ce4c","#ef3bb6","#d66551","#1a918f","#ff66fc","#2927c4","#7149af","#57e559","#8e3af4","#f9a270","#22547f","#db5e92","#edd05e","#6f25e8","#0dbc21","#280f7a","#6373ed","#5b910f","#7b34c1","#0cf29a","#d80fc1","#dd27ce","#07a301","#167275","#391c82","#2baeb5","#925bea","#63ff4f")

N_needed <- length(top10_union)
if (N_needed > length(total_color1)) {
  warning("颜色池不足！将循环使用颜色。建议扩展 total_color1。")
  color_pool <- rep(total_color1, ceiling(N_needed / length(total_color1)))
} else {
  color_pool <- total_color1[1:N_needed]
}

# 创建新的 color_mapping：仅针对这些重要物种
color_mapping <- setNames(color_pool, top10_union)
color_mapping["Others"] <- "#D3D3D3"

##################
dt_long <- dt3 %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "abundance") %>%
  left_join(sample_map %>% select(Sample, Country), by = "Sample")

# 计算每组内物种的平均丰度，并选出前N个物种
top_n_species_per_group <- function(df, n = 10) {
  df %>%
    group_by(Country, name) %>%
    summarise(avg_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop') %>%
    group_by(Country) %>%
    top_n(n, wt = avg_abundance)
}
xx <- top_species_result
top_species_result <- top_n_species_per_group(dt_long)
top_species_names <- unique(top_species_result$name)
#save(top_species_result,file="fun_Site.top_species_result.RData")
################
process_data_for_plotting <- function(dt, sample_map, top_species_df) {
  # 1. 转长格式
  dt_long <- dt %>%
    pivot_longer(cols = -name, names_to = "Sample", values_to = "abundance") %>%
    left_join(sample_map %>% select(Sample, Site), by = "Sample")
  
  # 2. 创建带标记的 top 表
  top_key <- top_species_df %>%
    select(Site, name) %>%
    distinct() %>%
    mutate(is_top_ref = TRUE)   # 标记这些是 top 物种
  
  # 3. 左连接并判断
  dt_annotated <- dt_long %>%
    left_join(top_key, by = c("Site", "name")) %>%
    mutate(
      is_top = !is.na(is_top_ref),
      name_clean = ifelse(is_top, name, "Others")
    ) %>%
    select(-name, -is_top_ref, -is_top) %>%
    rename(name = name_clean)
  
  # 4. 汇总
  dt_filtered <- dt_annotated %>%
    group_by(Sample, Site, name) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = 'drop')
  
  return(dt_filtered)
}

dt_processed <- process_data_for_plotting(
  dt = dt3,
  sample_map = sample_map,
  top_species_df = top_species_result   # 传入带 Site 的 top 表
)
#save(dt_processed,file="fun_Site.dt_processed.RData")
####################
dt_wide <- dt_processed %>%
  select(Sample, name, abundance) %>%          # 只保留绘图所需列
  pivot_wider(
    names_from = Sample,
    values_from = abundance,
    values_fill = 0                          # 缺失样本补0
  ) %>%
  column_to_rownames(var = "name")  

sample_map_for_plot <- dt_processed %>%
  select(Sample, Site) %>%
  distinct()  # 去重（每个样本只出现一次）
sample_map_for_plot <- as.data.frame(sample_map_for_plot)
rownames(sample_map_for_plot) <- sample_map_for_plot$Sample

p <- zy_group_compositions(
  dt = dt_wide,
  sample_map = sample_map_for_plot,
  ID = "Sample",
  group = "Site",        # 必须和 sample_map 中的列名一致
  top_N = nrow(dt_wide),       # 已预筛选，直接用全部
  title = "Airway Microbiome Composition by Site",
  taxo.color = color_mapping,
  rescale = FALSE,             # 假设 TPM 不需再归一化
  width = 0.85,
  order_func = "order",        # 或 "cluster"
  order_n = 1
)

p <- p + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)
print(p)

####################ko
tax <- rownames_to_column(taxo2,var="name")
tax <- tax %>%
  filter(.[[2]] != "Eukaryota" | (.[[2]] == "Eukaryota" & .[[3]] == "k_Fungi"))
tax <- tax[,1:2]

map <- read.table("../01.function/kegg.ko.map", sep="\t", header=F, check.names=F)
names(map) <- c("name","ko")

dt3 <- left_join(tax,map,by="name")
dt4 <- dt3[,-1]

names(dt4)[1] <- "kin"
dt4$kin[dt4$kin == "Archaea"] <- "Bacteria"
dt4 <- dt4 %>% filter(!is.na(ko))
setwd("/share/data1/limin/Airway/01.function/")
#save(dt4,file="veen.RData")
data_list <- split( dt4$ko,dt4$kin)

#options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#install.packages("ggvenn", dependencies = TRUE)
library(ggvenn)
p <- ggvenn(data = data_list,              # 传入处理好的列表
            show_elements = FALSE,         # 不显示具体元素名称
            show_percentage = TRUE,        # 显示百分比
            digits = 1,                    # 保留1位小数
            fill_color = c("#4EBBA6", "#F27D52", "#5C5C99"), # 沿用之前的配色
            stroke_size = 0.5,             # 边框粗细
            set_name_size = 5,             # 集合名称字体大小
            text_size = 4) +               # 内部数字字体大小
  ggtitle("Overlap of KOs Across Bacteria, Virus, and Fungi") +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p
