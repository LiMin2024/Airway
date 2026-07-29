setwd("/share/data1/limin/Airway/01.function/")
data <- myread("rgi.profiles")
numeric_cols <- sapply(data, is.numeric)
data$mean <- rowMeans(data[, numeric_cols], na.rm = TRUE)
data <- data%>%
  select(name,mean)

known_sum <- sum(data[data$name != "unknown", "mean"], na.rm = TRUE)
unknown_value <- data[data$name == "unknown", "mean"]
merged_data <- data.frame(
  name = c("unknown", "known"),
  mean = c(unknown_value, known_sum)
)
merged_data$name <- factor(merged_data$name, levels = c("unknown", "known"))

ggplot(merged_data, aes(x = "", y = mean, fill = name)) + # x轴为空字符串，集中显示
  geom_bar(stat = "identity", color="black", position = position_stack(reverse = TRUE)) + # 使用identity统计变换并反转堆叠顺序
  scale_fill_manual(values=c("known"="black", "unknown"="white")) + # 手动设置填充颜色，确保VFDB在上部
  labs(title = "Relative Proportions of Groups",
       x = NULL, # 移除x轴标题
       y = "Proportion (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # y轴标签格式化为百分比
  theme_minimal() +
  theme(axis.text.x = element_blank(), # 移除x轴文本
        axis.ticks.x = element_blank(), # 移除x轴刻度线
        legend.background = element_rect(fill="transparent"), 
        aspect.ratio = 4) # 添加边框

##########
merge <- data[-1,]
map <- read.table("rgi_category", header = T, fill = TRUE, sep = "\t",check.names = F)
map2 <- map[,-1]
map3 <- map2 %>%
  separate_rows('Drug Class', sep = ";\\s*")
map4 <- map3 %>% 
  distinct(`Drug Class`, ARO, .keep_all = TRUE)

merge2 <- merge(merge,map4,by.x = "name",by.y = "ARO")
merge3 <- merge2[,c(3,2)]
names(merge3)[1] <- "name"
merge4 <- merge3 %>%
  group_by(name) %>%           # 按第1列分组
  summarise(mean = sum(mean, na.rm = TRUE))

merge <- merge4[order(merge4$mean, decreasing = TRUE), ]

top10 <- merge[1:10, ]
others_sum <- sum(merge[11:nrow(merge), "mean"])
others_row <- data.frame(
  name = "Others",
  mean = others_sum
)
result <- rbind(top10, others_row)
result$Proportion <- result$mean / sum(result$mean) * 100
result$name <- factor(result$name, levels = result$name)

colors <- colorRampPalette(c("lightblue", "white"))(nrow(result))


result$label_name <- paste0(result$name, " (", sprintf("%.3f", result$mean), ")")

result$label_name <- factor(result$label_name, levels = result$label_name)

# 3. 绘图
p <- ggplot(result, aes(x = "", y = Proportion, fill = label_name)) + # <--- 这里改用 label_name
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = colors, name = NULL) + # name=NULL 去掉图例标题，更清爽
  labs(title = "Relative Proportions of Groups",
       x = NULL,
       y = "Proportion (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",      # 确保图例在右侧
    legend.background = element_rect(fill = "transparent"),
    aspect.ratio = 1.5,             # 根据实际画布调整长宽比
    legend.text = element_text(size = 10) # 可适当调整图例文字大小
  )

print(p)

####################
merge <- data[-1,]
map <- read.table("rgi_category", header = T, fill = TRUE, sep = "\t",check.names = F)
merge1 <- merge(merge,map,by.x = "name",by.y = "ARO")

gene.map <- read.table("rgi.map", header = F, fill = TRUE, sep = "\t",check.names = F)
merge2 <- merge(gene.map,merge1,by.x = "V1",by.y = "ORF_ID")

taxo <- read.table("geneset.taxonomy.known", fill = TRUE, header = FALSE, sep = "\t")
merge3 <- merge(merge2,taxo[,c(1,8)],by = "V1")

merge3 <- merge3[,c(4:6)]
names(merge3)[2:3] <- c("Function","Genus")

merge_dt <- merge3 %>%
  separate_rows(Function, sep = ";\\s*") %>% 
  group_by(Function, Genus) %>%      
  summarise(Count = sum(mean, na.rm = TRUE), .groups = "drop") %>% 
  arrange(desc(Count))
merge_dt <- merge_dt[merge_dt$Genus!= "Unclassified",]

# 2. 计算 Function 和 Genus 的排序（按总流量降序）
func_total <- merge_dt %>%
  group_by(Function) %>%
  summarise(total_count = sum(Count)) %>%
  arrange(desc(total_count))

top_n <- 10
top_functions <- func_total %>%
  slice(1:top_n) %>%
  pull(Function)

# 3. 创建映射：Function -> Group（前10保留，其余为 Others）
df_top10_plus_others <- func_total %>%
  mutate(Group = ifelse(Function %in% top_functions, Function, "Others")) %>%
  select(Function, Group)

top_functions <- df_top10_plus_others %>%
  filter(Group != "Others") %>%
  pull(Group)

func_species_count_filtered <- merge_dt %>%
  mutate(Function = ifelse(Function %in% top_functions, Function, "Others")) %>%
  group_by(Function, Genus) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  ungroup()

###############
species_total <- func_species_count_filtered %>%
  group_by(Genus) %>%
  summarise(total_count = sum(Count), .groups = 'drop') %>%
  arrange(desc(total_count))

top_species <- species_total %>%
  slice(1:20) %>%
  pull(Genus)

func_species_count_final <- func_species_count_filtered %>%
  mutate(Genus = ifelse(Genus %in% top_species, Genus, "Others")) %>%
  group_by(Function, Genus) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  ungroup()

###############
node_totals <- func_species_count_final %>%
  pivot_longer(cols = c(Function, Genus), names_to = "variable", values_to = "node") %>%
  group_by(node) %>%
  summarise(total_flow = sum(Count), .groups = 'drop') %>%
  arrange(desc(total_flow))

# 分离 Others 和非 Others
non_others <- node_totals %>%
  filter(node != "Others") %>%
  pull(node)

# ✅ 构建最终排序：非 Others 降序 + Others 在最后
desired_order <- c(non_others, "Others")

library(ggsankey) 
# 筛选实际存在的节点
sankey_data_temp <- func_species_count_final %>%
  make_long(x = "Function", next_x = "Genus", value = "Count")

existing_nodes <- unique(sankey_data_temp$node)
desired_order <- desired_order[desired_order %in% existing_nodes]
desired_order <- rev(desired_order)
# =================== 5. 构建桑基图数据 ===================
sankey_data <- func_species_count_final %>%
  make_long(x = "Function", next_x = "Genus", value = "Count")

# 设置 node 因子顺序（决定垂直位置）
sankey_data$node <- factor(sankey_data$node, levels = desired_order, ordered = TRUE)

# 设置 x 轴顺序
sankey_data$x <- forcats::fct_relevel(sankey_data$x, "Function", "Genus")

# =================== ✅ 6. 配色：Others 用灰色，其余用彩虹色 ===================
all_nodes <- levels(sankey_data$node)
n_colors <- length(all_nodes)

# 你的彩虹色
node_colors <- c(
  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
  "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#ed1299",
  "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e", "#4aef7b", "#e86502",
  "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f", "#ce2523",
  "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
  "#1a918f", "#ff66fc", "#2927c4", "#7149af", "#57e559", "#8e3af4", "#f9a270", "#22547f",
  "#db5e92", "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f", "#7b34c1",
  "#0cf29a", "#d80fc1", "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5", "#925bea", "#63ff4f"
)

n_rainbow_needed <- if ("Others" %in% all_nodes) {
  n_colors - 1
} else {
  n_colors
}

# 3. 截取或循环颜色 (修正了这里的逻辑)
if (n_rainbow_needed > length(node_colors)) {
  # 如果需要的颜色比色板多，则循环使用色板
  node_colors <- rep(node_colors, length.out = n_rainbow_needed)
} else {
  # 如果色板够用，直接截取
  node_colors <- node_colors[1:n_rainbow_needed]
}

# 4. 构建最终映射
# 找出非-Others的节点
non_others_nodes <- all_nodes[all_nodes != "Others"]

# 确保颜色数量与节点数量一致 (Debug 检查)
if (length(non_others_nodes) != length(node_colors)) {
  stop("颜色数量与节点数量不匹配！请检查 all_nodes 是否有重复值。")
}

# 构建向量
node_color_mapping <- c(
  setNames(node_colors, non_others_nodes),
  "Others" = "gray50"
)


sankey_data$node <- fct_rev(sankey_data$node)
# =================== 7. 绘图 ===================
p <- ggplot(sankey_data,
            aes(x = x,
                next_x = next_x,
                node = node,
                next_node = next_node,
                value = value,
                fill = node)) +
  
  geom_sankey(flow.alpha = 0.8,
              node.color = "black",
              node.linewidth = 0.5,
              show.legend = FALSE) +
  
  # --- 修改点 1: 添加 angle 参数旋转标签 ---
  geom_sankey_label(aes(label = node),
                    hjust = 0.5,
                    size = 5,
                    color = "black",
                    fill = "white",
                    alpha = 0.7,
                    show.legend = FALSE) +
  
  scale_fill_manual(values = node_color_mapping, guide = "none") +
  
  scale_x_discrete(expand = expansion(add = c(0.2, 0.6)),
                   labels = c("x" = "Function", "next_x" = "Species")) +
  
  theme_sankey(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    aspect.ratio = 1.2
  ) +
  
  labs(
    title = "Sankey Diagram: Function → Species",
    subtitle = "Nodes ordered by flow (Others at bottom, gray)",
    x = ""
  )
print(p)

###########
data2 <- data[-1,]
rownames(data2) <- NULL
data3 <- data2 %>%
  column_to_rownames(var = "name") %>% 
  as.data.frame() %>%  
  mutate(across(everything(), as.numeric))

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
dt_norm <- norm_data(data3)
#save(dt_norm,file="rgi.known.RData")