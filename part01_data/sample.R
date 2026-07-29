library(ggplot2)
library(sf)         # 用于处理地图数据
library(dplyr)      # 用于数据处理
library(ggforce)    # 用于绘制饼图 (geom_arc_bar)
library(viridis)    # 用于配色
library(scales)  
library(maps)
library(readxl)
setwd("/share/data1/limin/Airway/00.data/")

data <- read.table("sample.info.txt",header=TRUE,sep="\t")

data_summary <- data %>%
  group_by(country, site) %>%
  summarise(Count = n(), .groups = "drop") %>%
  # 计算每个国家的总数 (用于气泡大小)
  group_by(country) %>%
  mutate(Total_Samples = sum(Count)) %>%
  ungroup()
#write.table(data_summary,"clean_metadata",sep = "\t",quote = F,row.names = F)

data_pie <- data_summary %>%
  select(country, site, Count, Total_Samples) %>%
  pivot_wider(names_from = site, values_from = Count, values_fill = 0)

####
world_map <- map_data("world")
world_map$region <- tolower(world_map$region)
data_summary$country <- tolower(data_summary$country)
# 3. 合并数据 (使用 left_join 更安全)
map_data_merged <- left_join(world_map, data_summary, by = c("region" = "country"))

# 4. 绘图
ggplot() +
  # 绘制地图背景
  geom_polygon(
    data = map_data_merged,
    aes(x = long, y = lat, group = group, fill = is.na(Count)), # 根据是否有数据填充颜色
    color = "grey", # 边框颜色设为灰色，避免黑色线条过于抢眼
    size = 0.1
  ) +
  # 如果你有具体的点数据要画，用 geom_point
  # geom_point(data = data_summary, aes(x = lon, y = lat, size = Total_Samples), color = "red") +
  
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "steelblue"), guide = FALSE) +
  coord_fixed(1) +
  theme_void()

###############
ggplot() +
  geom_polygon(
    data = map_data_merged,
    aes(x = long, y = lat, group = group, fill = Total_Samples),
    color = "grey", size = 0.1
  ) +
  
  # 【关键修改】使用 scale_fill_stepsn 替代 gradientn
  scale_fill_stepsn(
    colors = c("#CBE5D6", "#4292C6", "#08306B","#FDFCF0" ),
    values = scales::rescale(c(100, 500, 1000, 2000)),
    breaks = c(100, 500, 1000, 2000),
    labels = c("100", "500", "1000", "2000"),
    limits = c(0, 10000),
    
    na.value = "white"
  ) +
  
  coord_fixed(1) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text.position = "bottom", # 配合 stepsn 效果更好
    legend.title.position = "top"
  ) +
  guides(fill = guide_colorsteps(
    direction = "horizontal",
    barwidth = 10,
    barheight = 0.5,
    show.limits = TRUE # 显示起止刻度
  ))

####
label_data <- map_data_merged %>%
  filter(!is.na(Count)) %>%
  distinct(region, .keep_all = TRUE)

ggplot() +
  # 1. 绘制地图背景
  geom_polygon(
    data = map_data_merged,
    aes(x = long, y = lat, group = group, fill = is.na(Count)),
    color = "grey",
    size = 0.1
  ) +
  
  # 2. 添加国家名称标注 (核心修改部分)
  geom_text(
    data = label_data, # 使用上面准备好的去重数据
    aes(x = long, y = lat, group = group, label = region), # label = region 指定显示国家名
    size = 3,          # 字体大小，可以根据需要调整
    color = "black",   # 字体颜色
    fontface = "bold", # 字体加粗，让名字更显眼
    alpha = 0.8        # 透明度，防止遮挡地图细节
  ) +
  
  # 3. 颜色设置
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "steelblue"), guide = FALSE) +
  
  # 4. 坐标与主题
  coord_fixed(1.3) +
  theme_void()

##########
library(ggplot2)
library(dplyr)
library(scales) # 用于百分比格式化
library(ggplot2)
library(dplyr)
library(scales)
library(ggplot2)
library(dplyr)
library(scales)
unique(data_summary$country)

data_summary2 <- data_summary %>%
  group_by(country, site) %>%  # 按大洲和采样部位分组
  summarise(
    Count = sum(Count),           # 汇总每个部位的数量
    Total_Samples = first(Total_Samples),  # 取该洲的总样本数（同一洲内相同）
    .groups = "drop"
  ) %>%
  # 重新排序，让结果更清晰
  arrange(country, site)

# 1. 数据预处理
plot_data <- data_summary2 %>%
  group_by(country) %>%
  mutate(
    Country_Total = sum(Count),
    Fraction = Count / Country_Total,
    Y_Position = cumsum(Fraction) - Fraction/2,
    Country_Label = paste0(country, "\n(", Country_Total, ")")
  ) %>%
  ungroup()

# 2. 计算每个部位的总样本数
site_totals <- data_summary %>%
  group_by(site) %>%
  summarise(Site_Total = sum(Count), .groups = "drop")

# 3. 创建自定义图例标签
# 将部位名称和总数合并
legend_labels <- paste(site_totals$site, "(", site_totals$Site_Total, ")")
names(legend_labels) <- site_totals$site

color_pool <- c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                "#6a3d9a", "#b15928","#8dd3c7",
                "#fb8072", "#80b1d3",
                "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
                "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4",
                "#b2df8a", "#33a02c", "#fb9a99")
sites <- site_totals$site
my_colors <- setNames(color_pool[1:length(sites)], sites)
#save(my_colors,file="site.color.RData")
# 5. 绘图
p <- ggplot(plot_data, aes(x = "", y = Fraction, fill = site)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_wrap(~ Country_Label, ncol = 4) +
  geom_text(aes(y = Y_Position, label = percent(Fraction, accuracy = 1)), 
            color = "white", size = 3) +
  theme_void() +
  theme(
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Body Site Distribution by Country", fill = "Site Category") +
  # 应用自定义颜色和标签
  scale_fill_manual(
    values = my_colors, 
    labels = legend_labels, # 使用自定义标签
    name = "site"
  )

# 显示图形
print(p)

###########map.rate
setwd("/share/data1/limin/Airway/00.data/")
data <- read.table("sample.info.txt",header=TRUE,sep="\t")
map <- read.table("map.rate",header=FALSE,sep="\t")#map.rate
names(map) <- c("Sample_ID","rate")
map2 <- map[map$Sample_ID %in% data$Sample_ID,]

not_in_data <- map$Sample_ID[!map$Sample_ID %in% data$Sample_ID]
not_in_map <- data$Sample_ID[!data$Sample_ID %in% map$Sample_ID]

data2 <- left_join(data,map2,by="Sample_ID")
data_clean <- data2[!is.na(data2$rate), ]
data_clean <- data_clean[,c(3,5,18)]
data_clean$rate_num <- as.numeric(gsub("%", "", data_clean$rate))

Air_data_clean <- data_clean
rmgc_data_clean <- data_clean

Air_data_clean$Source <- "iHAMGC"
rmgc_data_clean$Source <- "RMGC"

# 合并数据
combined_data <- rbind(Air_data_clean, rmgc_data_clean)

combined_data2 <- combined_data[combined_data$country != "China"&combined_data$country != "USA"&combined_data$country != "Germany",]
mean_rates <- combined_data2 %>%
  group_by(Source) %>%       # 按 site 和 Source 联合分组
  summarise(
    Mean_Rate = mean(rate_num, na.rm = TRUE),  # 计算平均值，忽略缺失值
    Count = n()                              # 顺便统计每组的样本量（可选）
  ) %>%
  ungroup()  

median_site_table <- combined_data %>%
  group_by(site,Source) %>%
  summarise(Median_Map_Rate = median(rate_num, na.rm = TRUE))

df_wide <- median_country_table %>%
  pivot_wider(names_from = Source, values_from = Median_Map_Rate)
df_ratio <- df_wide %>%
  mutate(Ratio_iHAMGC_vs_RMGC = iHAMGC / RMGC)
df_ratio <- df_ratio[-3,]

# 2. 按 Country 统计 Map Rate 的中位数
median_country_table <- combined_data %>%
  group_by(country,Source) %>%
  summarise(Median_Map_Rate = median(rate_num, na.rm = TRUE))

df_wide <- pivot_wider(
  data = median_country_table, 
  names_from = Source,      # 指定哪一列的值会变成新的列名
  values_from = Median_Map_Rate # 指定用哪一列的数据来填充新列的单元格
)
df_wide <- df_wide %>%
  mutate(
    Fold_Change = iHAMGC / RMGC  # 将 iHAMGC 除以 RMGC，结果存入新列 Fold_Change
  )

p_all <- ggplot(combined_data, aes(x = "All Samples", y = rate_num, fill = Source)) +
  
  # 1. 箱线图：使用 position_dodge 并排错开，隐藏默认异常值
  geom_boxplot(
    alpha = 0.7, 
    color = "black", 
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  # 3. 统一配色（与 p_site_combined 保持一致）
  scale_fill_manual(values = c("iHAMGC" = "lightblue", "RMGC" = "lightcoral")) +
  
  labs(title = "", x = "", y = "Map Rate (%)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 19, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 19, color = "black"),
    axis.title = element_text(size = 19),
    plot.title = element_text(size = 19, hjust = 0.5),
    aspect.ratio = 17,
    legend.position = "none"
  )

p_site_combined <- ggplot(combined_data, aes(x = site, y = rate_num, fill = Source)) +  
  
  # 1. 绘制箱线图：使用 position_dodge() 使不同 Source 的箱子并排错开
  geom_boxplot(
    alpha = 0.7, 
    color = "black", 
    position = position_dodge(width = 0.8),
    outlier.shape = NA  # 【建议】隐藏箱线图自带的异常值，避免和下面的散点重叠
  ) +  
  
  # 3. 自定义颜色
  scale_fill_manual(values = c("iHAMGC" = "lightblue", "RMGC" = "lightcoral")) + 
  
  # 5. 主题美化
  theme_bw() +                         
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(size = 19, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(), 
    aspect.ratio = 1.7,
    legend.position = "none"
  )

p_country_combined <- ggplot(combined_data, aes(x = country, y = rate_num, fill = Source)) +  
  
  # 1. 绘制箱线图：使用 position_dodge() 使不同 Source 的箱子并排错开
  geom_boxplot(
    alpha = 0.7, 
    color = "black", 
    position = position_dodge(width = 0.8),
    outlier.shape = NA  # 【建议】隐藏箱线图自带的异常值，避免和下面的散点重叠
  ) +  
  
  # 3. 自定义颜色
  scale_fill_manual(values = c("iHAMGC" = "lightblue", "RMGC" = "lightcoral")) + 
  
  # 5. 主题美化
  theme_bw() +                         
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(size = 19, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(), 
    aspect.ratio = 1,
    legend.position = "top"
  )

final_plot <- p_all + p_site_combined + p_country_combined
print(final_plot)

###############pip
# 1. 加载绘图包
library(ggplot2)

#airway
df <- data.frame(
  group = c("Bacteria", "Eukaryota","Virus"),
  value = c(15211482, 1533924, 360688)
)

df <- df2

# 3. 数据预处理：计算百分比和标签位置
# 计算百分比（保留两位小数）
df$percentage <- round(df$value / sum(df$value) * 100, 2)
# 计算标签在饼图中的垂直位置（保证标签显示在每个扇形的中间）
df$ypos <- cumsum(df$percentage) - 0.5 * df$percentage

# 4. 绘制饼图
p <- ggplot(df, aes(x = "", y = percentage, fill = group)) +
  # 绘制堆叠柱状图，width=1表示柱子之间没有间隙，color="white"给扇形之间加白色边框
  geom_bar(stat = "identity", width = 1, color = "white") +
  # 核心步骤：将直角坐标系转换为极坐标系，生成饼图
  coord_polar("y", start = 0) +
  # 添加百分比标签，并设定标签颜色为白色，大小为4
  geom_text(aes(y = ypos, label = paste0(percentage, "%")), color = "white", size = 4) +
  # 使用空白主题，去掉背景、网格和坐标轴，让饼图更干净
  theme_void() +
  # 自定义填充颜色（你可以换成自己喜欢的颜色）
  scale_fill_brewer(palette = "Set2")

# 5. 展示图片
print(p)

##############curve
dt = read.table("Air.roc.data.tsv", sep="\t", header=T)

dtf <- dt %>%
  group_by(nspecies, group) %>%
  summarise(value = mean(obs)) %>%
  data.frame()

dtf_filtered <- dtf %>%
  filter(group == "all" | (group == "nosingle"))

p1 <- ggplot(data = dtf_filtered, aes(x = nspecies / 100000, y = value / 100000, color = group)) + # 1. Y轴数值除以1000
  geom_line(size = 1) +
  #geom_point()+ 
  theme_bw() +
  theme(
    aspect.ratio = 0.6,
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.grid = element_blank() # 2. 去掉所有网格线
  ) +
  scale_x_continuous(breaks = seq(0, max(dtf_filtered$nspecies, na.rm = T), by = 1000)) +
  xlab("Number of genes (million)") +
  ylab("Number of non-redundant genes (million)") 
p1

############
data <- read.table("kegg.profiles", header = T,  sep = "\t",check.names = F)
result <- data %>%
  mutate(average = rowMeans(select(., -name), na.rm = TRUE)) %>% 
  select(name, average)
#result <- result[-1,]
#result$average_100 <- (result$average / sum(result$average, na.rm = TRUE)) * 100
sum(result$average, na.rm = TRUE)

map <- read.table("KO_level_A_B_C_D_Description", 
                  header = F, 
                  fill = TRUE, 
                  sep = "\t", 
                  check.names = FALSE, 
                  quote = "") 
map <- map[,c(2,6)]
names(map) <- c("pathway","name")
map_unique <- map[!duplicated(map[, c("pathway", "name")]), ]
merge <- merge(result,map_unique,by="name")
merge <- merge[!grepl("Brite", merge$pathway), ]

df_top50 <- merge %>% 
  arrange(desc(average)) %>%   # 第一步：按 average 降序排列
  slice_head(n = 50) 
df_top50 <- df_top50 %>% 
  mutate(name = gsub("\\s*\\[.*?\\]", "", name))
df_top50$name <- factor(df_top50$name, levels = df_top50$name)

# 4. 定义颜色映射 (参考原图配色)
custom_colors <- c(
  "Metabolism" = "#FFCC99",                 # 浅橙色
  "Environmental Information Processing" = "#CDA4DE", # 紫色
  "Cellular Processes" = "#99CC66",         # 绿色
  "Genetic Information Processing" = "#99CCFF", # 蓝色
  "Human Diseases" = "#FF9999"              # 粉红色
)

# 5. 绘图
p <- ggplot(df_top50, aes(x = name, y = average, fill = pathway)) +
  geom_bar(stat = "identity", width = 0.8) +   # 绘制柱状图
  scale_fill_manual(values = custom_colors) +  # 应用自定义颜色
  theme_bw(base_size = 12) +                   # 使用黑白主题作为基础
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11), # X轴文字倾斜45度
    axis.text.y = element_text(size = 11),  
    legend.position = "right",                     # 图例放在顶部
    legend.title = element_blank(),              # 去掉图例标题
    panel.grid = element_blank(), 
    panel.border = element_blank(),              # 去掉边框
    axis.line = element_line(color = "black"),    # 保留坐标轴线
    aspect.ratio = 0.3
  ) +
  labs(
    x = NULL,                                    # 去掉X轴标题
    y = "Relative abundance (%)",                # Y轴标题
    title = "KEGG Pathway Enrichment Analysis"   # 可选的标题
  )
p

##########summary
sample_map = read_xlsx("part - sample.xlsx",sheet = 1)

agg_df <- sample_map %>%
  group_by(country, site) %>%
  summarise(
    Total_Samples = n_distinct(Sample_ID),       
    Total_BioProjects = n_distinct(BioProject_ID), 
    Cleandata_amount = sum(`Number of bases (Gbp)`, na.rm = TRUE), 
    .groups = "drop"                             
  )
agg_df1 <- agg_df[,1:3]
agg_df2 <- agg_df[,-3]

df_wide <- agg_df1 %>%
  pivot_wider(
    names_from = site,         # 指定哪一列的值作为新的列名
    values_from = Total_Samples, # 指定用哪一列的值来填充新表格
    values_fill = list(Total_Samples = 0), # 将缺失值填充为 0
    names_sort = TRUE          # 按首字母对新生成的列名进行排序
  )

result1 <- agg_df2 %>%
  group_by(country) %>%
  summarise(across(c(Total_BioProjects, Cleandata_amount), sum, na.rm = TRUE))

result2 <- agg_df2 %>%
  group_by(site) %>%
  summarise(across(c(Total_BioProjects, Cleandata_amount), sum, na.rm = TRUE))

