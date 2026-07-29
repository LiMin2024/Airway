library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(readxl)
options(warn=-1)

#根据给定的分组变量 group 和一个阈值条件（可以是绝对数量 cut_num 或者比例 cut_rate），过滤或重新标记样本映射数据框 sample_map 中的组。
#这个函数特别适用于需要对某些组进行最小样本量限制或者将不符合样本量要求的组合并为一个新的类别（如 "less_X"，其中 X 是设定的阈值）的情况。
filter_group <- function(sample_map, group=NA, cut_rate = NA, cut_num=NA){
  cut_off = cut_num
  if(is.finite(cut_rate)){
    if(cut_rate>1 || cut_rate < 0){
      stopifnot("cut_rate should range(0,1)" = 1)
    }else{
      cut_off = sample_map %>%
        dplyr::select(!!sym(group)) %>%
        dplyr::group_by(!!sym(group)) %>%
        dplyr::summarise(count=n()) %>%
        pull(count) %>%
        max() * cut_rate
      cut_off = as.integer(cut_off)
    }
  }
  
  select_grps <- sample_map %>%
    dplyr::select(!!sym(group)) %>%
    dplyr::group_by(!!sym(group)) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::filter(count > cut_off) %>%
    dplyr::pull(!!sym(group))
  message("???????????")
  sample_map %>%
    dplyr::mutate( !!sym(group) := ifelse( !!sym(group) %in% select_grps, !!sym(group), paste("less_",cut_off, sep="")))
}

align_dist_sample <- function(dist, sample_map, ID=NA){
  dist = as.matrix(dist)
  intersect_id = intersect(sample_map[,ID], colnames(dist))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dist = dist[ sample_map[,ID], sample_map[,ID] ]
  list(dist=as.dist(dist), sample_map=sample_map)
}

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

zy_pcoa <- function(dt = NA, sample_map = NA, group = NA, ID = NA, sample.color = NULL,
                    ado_method = "bray", pca_method = "bray",
                    levels = 0.95, star_plot = FALSE, ellipse_plot = TRUE,
                    cut_rate = NA, cut_num = NA,
                    title = "PCoA", x = 1, y = 2, ados = TRUE, mydist = NULL,
                    return_dist = FALSE  # 新增参数：是否返回距离矩阵
){
  
  # 对齐数据
  if (!is.null(mydist)) {
    fmt_profile = align_dist_sample(mydist, sample_map, ID = ID)
    mydist = fmt_profile$dist
  } else {
    fmt_profile = align_dt_sample(dt, sample_map, ID = ID)
    dt = fmt_profile$dt
  }
  sample_map = fmt_profile$sample_map
  
  # 过滤分组（可选）
  if (is.finite(cut_rate) || is.finite(cut_num)) {
    sample_map = filter_group(sample_map, group = group, cut_rate = cut_rate, cut_num = cut_num)
    # 重新对齐 dt（如果 sample_map 被过滤）
    if (!is.null(dt)) {
      common_samples <- intersect(colnames(dt), sample_map[[ID]])
      dt <- dt[, common_samples, drop = FALSE]
      sample_map <- sample_map[sample_map[[ID]] %in% common_samples, , drop = FALSE]
      sample_map <- sample_map[match(colnames(dt), sample_map[[ID]]), , drop = FALSE]
    }
  }
  
  ## 颜色设置
  if (is.null(sample.color)) {
    sample.color = rainbow(length(unique(sample_map[[group]])))
  }
  
  # 统计每组样本数，用于图例
  group_summ <- sample_map %>%
    dplyr::select(all_of(group)) %>%
    dplyr::group_by(across(all_of(group))) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::mutate(new_label = paste(!!sym(group), " (", count, ")", sep = ""))
  
  new_label <- structure(group_summ$new_label, names = as.character(group_summ[[group]]))
  
  message(paste(length(unique(sample_map[[group]])), "of groups to plot"))
  
  # 计算或使用距离矩阵
  if (is.null(mydist)) {
    # 关键：确保 dt 是数值型矩阵
    dt <- as.matrix(dt)
    if (any(is.na(dt))) {
      warning("dt contains NA, replacing with 0.")
      dt[is.na(dt)] <- 0
    }
    if (min(dt) < 0) {
      stop("Abundance values must be non-negative for Bray-Curtis distance.")
    }
    mydist <- vegdist(t(dt), method = pca_method)
  } else {
    mydist <- as.dist(mydist)
  }
  
  # adonis 检验
  ado_r2 <- ado_p <- NA
  if (isTRUE(ados) && length(unique(sample_map[[group]])) > 1) {
    ado <- adonis2(mydist ~ sample_map[[group]])
    ado_r2 <- round(ado$R2[1], digits = 4)
    ado_p <- ado$`Pr(>F)`[1]
  }
  
  # PCoA 分析
  pcoa <- cmdscale(mydist, k = min(10, nrow(as.matrix(mydist)) - 1), eig = TRUE)
  eigs <- signif(pcoa$eig / sum(pcoa$eig), 4) * 100
  point <- pcoa$points
  colnames(point) <- paste0("pcoa.", 1:ncol(point))
  
  xlab <- paste0("PCoA ", x, " (", eigs[x], "%)")
  ylab <- paste0("PCoA ", y, " (", eigs[y], "%)")
  
  substitle <- parse(text = paste0("'R'^2 == '", ado_r2, "'~~italic(p) == '", ado_p, "'"))
  
  dm <- merge(as.data.frame(point), sample_map, by.x = 'row.names', by.y = ID)
  
  # 绘图
  p1 <- ggscatter(data = dm, 
                  x = paste0("pcoa.", x), 
                  y = paste0("pcoa.", y),
                  color = group,
                  star.plot = star_plot,
                  ellipse.level = levels, 
                  ellipse = ellipse_plot) +
    theme_bw() +
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black", linewidth = 0.25),
          panel.border = element_rect(colour = "black", linewidth = 0.25)) +
    scale_fill_manual(values = sample.color, guide = "none") +
    scale_color_manual(values = sample.color, labels = new_label) +
    labs(x = xlab, y = ylab, title = title, subtitle = substitle)
  
  # 构建返回结果
  result <- list(
    plot = p1,
    new_label = new_label,
    pcoa_points = point,
    sample_map = sample_map
  )
  
  # 如果要求返回距离矩阵
  if (return_dist) {
    result$dist <- mydist  # 添加距离矩阵
  }
  
  return(result)
}

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

##############
result <- align_dt_sample(dt, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map
aligned_dt2 <- t(aligned_dt)
#write.table(aligned_dt2,file = "all.kegg.aligned_dt", sep = "\t", row.names = T, quote = FALSE)

##################adonis
dist_matrix <- read.table("kegg.dist.txt", header = TRUE, row.names = 1,fill = TRUE)
kegg_mydist <- as.dist(dist_matrix)  # 转为 dist 对象

dist_matrix <- read.table("rgi.dist.txt", header = TRUE, row.names = 1,fill = TRUE)
rgi_mydist <- as.dist(dist_matrix)  # 转为 dist 对象

mydist <- dist_matrix

###########
table(sample_map$Site,sample_map$Country)

data_samples_raw <- sample_map[sample_map$Site == "Throat", ]#"Sputum", "BALF", "Nasopharynx", "Oropharynx", 
#                                                                 "Throat", "Nose"
data_samples_raw <- sample_map[sample_map$Country == "Germany", ]#"China", "USA", "Germany"
data_sample_ids <- trimws(as.character(data_samples_raw$Sample))

# 2. 匹配索引并提取子矩阵
matrix_ids <- trimws(rownames(as.matrix(dist_matrix)))
valid_samples <- data_sample_ids[data_sample_ids %in% matrix_ids]
row_indices <- match(valid_samples, matrix_ids)
row_indices <- row_indices[!is.na(row_indices)]

# 提取子矩阵
sub_mat <- as.matrix(dist_matrix)[row_indices, row_indices]

# 3. 【核心新增】使用 valid_samples 重新对齐元数据表！
# 这一步确保 meta 数据的行数、顺序与提取出的距离矩阵完全一致
data_samples <- data_samples_raw[match(valid_samples, data_sample_ids), ] 

# 4. 清洗矩阵 (处理 NA 和负数)
if(sum(is.na(sub_mat)) > 0) {
  cat("警告：检测到距离矩阵中有缺失值(NA)，已自动替换为 0。\n")
  sub_mat[is.na(sub_mat)] <- 0
}
sub_mat[sub_mat < 0] <- 0
data_dist_matrix_clean <- as.dist(sub_mat)

# 5. 验证行数是否绝对相等（建议保留此行作为安全检查）
cat("距离矩阵行数:", nrow(sub_mat), "\n")
cat("元数据表行数:", nrow(data_samples), "\n")

#运行 adonis.R

########################
load("Country_Germany_ado.Rdata")
ado_r2 <- round(ado_result$R2[1], 4)
ado_p <- ado_result$`Pr(>F)`[1]

mydist <- data_dist_matrix_clean
# 格式化 p 值（比如 p < 0.001）
p_label <- ifelse(ado_p < 0.001, "< 0.001", round(ado_p, 3))
pcoa_plot_result <- zy_pcoa(
  dt = NULL,                  
  mydist = mydist,            
  sample_map = data_samples,
  group = "Site",          
  ID = "Sample",           
  
  sample.color = NULL,        
  pca_method = "bray",        
  title = "Oropharynx",    # 建议加上标题
  
  x = 1, y = 2,               
  star_plot = FALSE,
  ellipse_plot = TRUE,
  levels = 0.95,              
  
  ados = FALSE,               # 关键：关闭内部 adonis 计算
  
  return_dist = FALSE
)

#save(pcoa_plot_result,file="kegg/Germany_PCoA.Rdata")

##############
#load("kegg/Germany_PCoA.Rdata")
p <- pcoa_plot_result$plot+
  theme(axis.text.x = element_text(size = 17), # 设置x轴文本大小
        axis.text.y = element_text(size = 17), # 设置y轴文本大小
        axis.title.x = element_text(size = 19), # 设置x轴标题大小
        axis.title.y = element_text(size = 19),
        aspect.ratio = 1 )

target_sites <- c("Sputum", "Nasopharynx", "Throat")
named_colors <- c(
  "Sputum" = "#a6cee3", 
  "Nasopharynx" = "#b2df8a", 
  "Throat" = "#fb9a99"
)

filtered_data <- p$data %>% 
  dplyr::filter(Site %in% target_sites)

# 3. 整合绘图与美化操作
p_final <- p %+% filtered_data +   
  geom_point(aes(color = Site), size = 1) +             
  scale_color_manual(values = named_colors) +  
  labs(
    title = "China",         # 设置主标题文字
    subtitle = bquote(R^2 == .(ado_r2) ~ "," ~ italic(p) == .(p_label))) +
  theme(
    axis.text.x = element_text(size = 17), 
    axis.text.y = element_text(size = 17), 
    axis.title.x = element_text(size = 19), 
    axis.title.y = element_text(size = 19),
    aspect.ratio = 1
  )

print(p_final)
p1 <- p_final
p2 <- p_final
p3 <- p_final
library(patchwork)
combined_plot <- p1 + p2 + p3

ggsave("kegg.site.β.pdf", plot = combined_plot, width = 20, height = 20, dpi = 300)

#############
#load("kegg/Throat_PCoA.Rdata")

my_colors2 <- c("#e31a1c", "#fdbf6f", "#cab2d6")

p <- pcoa_plot_result$plot+
  theme(axis.text.x = element_text(size = 17), # 设置x轴文本大小
        axis.text.y = element_text(size = 17), # 设置y轴文本大小
        axis.title.x = element_text(size = 19), # 设置x轴标题大小
        axis.title.y = element_text(size = 19),
        aspect.ratio = 1 )

target_sites <- c("China", "USA", "Germany")
named_colors <- c(
  "China" = "#e31a1c", 
  "USA" = "#fdbf6f", 
  "Germany" = "#cab2d6"
)

filtered_data <- p$data %>% 
  dplyr::filter(Country %in% target_sites)

# 3. 整合绘图与美化操作
p_final <- p %+% filtered_data +   
  geom_point(aes(color = Country), size = 1) +             
  scale_color_manual(values = named_colors) +  
  labs(
    title = "Nasopharynx",         # 设置主标题文字
    subtitle = bquote(R^2 == .(ado_r2) ~ "," ~ italic(p) == .(p_label))) +
  theme(
    axis.text.x = element_text(size = 17), 
    axis.text.y = element_text(size = 17), 
    axis.title.x = element_text(size = 19), 
    axis.title.y = element_text(size = 19),
    aspect.ratio = 1
  )

print(p_final)

p1 <- p_final
p2 <- p_final
p3 <- p_final
library(patchwork)
combined_plot <- p1 + p2 + p3

ggsave("kegg.country.β.pdf", plot = combined_plot, width = 20, height = 20, dpi = 300)

#############
load("China_PCoA.Rdata")
dm <- filtered_data 

# 更新因子水平顺序（根据实际需要调整这6个的顺序）
dm$Site <- factor(dm$Site, 
                  levels = c("Sputum", "BALF", "Nasopharynx", "Oropharynx", "Throat", "Nose"))

# 使用对应的 6 色配色方案
my_colors_site <- c(
  "Sputum" = "#a6cee3", 
  "BALF" = "#1f78b4", 
  "Nasopharynx" = "#b2df8a", 
  "Oropharynx" = "#33a02c", 
  "Throat" = "#fb9a99", 
  "Nose" = "#6a3d9a"
)

# 生成两两比较的组合
groups <- levels(dm$Site)
comparisons <- combn(groups, 2, simplify = FALSE)

max_y2 <- max(dm$pcoa.2, na.rm = TRUE)
y_positions_side <- seq(max_y2 * 1.2, max_y2 * 2.5, length.out = length(comparisons))

# 动态计算下方水平箱线图 (PCoA1) 的位置
max_x1 <- abs(min(dm$pcoa.1, na.rm = TRUE))
y_positions_horizontal <- seq(max_x1 * 1.2, max_x1 * 2.5, length.out = length(comparisons))

# ====== 绘制右侧竖直箱线图 (PCoA2) ======
p_side <- ggplot(dm, aes(x = Site, y = pcoa.2, fill = Site)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) + 
  geom_signif(
    comparisons = comparisons,
    y_position = y_positions_side,
    map_signif_level = TRUE,
    textsize = 3, tip_length = 0.015, linewidth = 0.5, vjust = -0.2
  ) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    aspect.ratio = 2
  ) +
  labs(y = "PCoA2")

# ====== 绘制下方水平箱线图 (PCoA1) ======
p_bottom <- ggplot(dm, aes(x = Site, y = pcoa.1, fill = Site)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) +
  geom_signif(
    comparisons = comparisons,
    y_position = y_positions_horizontal,
    map_signif_level = TRUE,
    textsize = 3, tip_length = 0.015, linewidth = 0.5, vjust = -0.2
  ) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    aspect.ratio = 0.5
  ) +
  labs(x = "Site") +
  coord_flip()

# ====== 绘制中间的主散点图 ======
p_final <- ggplot(dm, aes(x = pcoa.1, y = pcoa.2, color = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.5) +
  scale_color_manual(values = my_colors_site) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x.bottom = element_blank(),  # 移除底部x轴标题（因x轴移到顶部）
    aspect.ratio = 1
  ) +
  scale_x_continuous(position = "top") +  # x轴刻度/标签移到顶部
  labs(x = "PCoA 1", y = "PCoA 2")

final_plot <- (
  (p_final + p_side) /        
    (p_bottom + plot_spacer())     
) +
  plot_layout(
    widths = c(9, 1),   
    heights = c(9, 3)   
  ) +
  theme(plot.margin = margin(10, 10, 10, 10))

print(final_plot)

####################无显著性差异版
load("China_PCoA.Rdata")
dm <- filtered_data 

# 更新因子水平顺序（根据实际需要调整这6个的顺序）
dm$Site <- factor(dm$Site, 
                  levels = c("Sputum", "BALF", "Nasopharynx", "Oropharynx", "Throat", "Nose"))

my_colors_site <- c(
  "Sputum" = "#a6cee3", "BALF" = "#1f78b4", "Nasopharynx" = "#b2df8a",
  "Oropharynx" = "#33a02c", "Throat" = "#fb9a99", "Nose" = "#6a3d9a"
)

# ====== 1. 绘制右侧竖直箱线图 (PCoA2) ======
p_side <- ggplot(dm, aes(x = Site, y = pcoa.2, fill = Site)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 30, l = 5.5),  # 侧图边距
    aspect.ratio = 2
  ) +
  labs(y = "PCoA2")

# ====== 2. 绘制下方水平箱线图 (PCoA1) ======
p_bottom <- ggplot(dm, aes(x = pcoa.1, y = Site, fill = Site)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5.5, r = 30, b = 5.5, l = 55),  # 底图边距（重点调整！）
    aspect.ratio = 0.5  # 可根据显示效果微调（如0.45/0.55）
  ) +
  labs(x = "PCoA1")

# ====== 3. 绘制主图 (PCoA Scatter) ======
p_final <- ggplot(dm, aes(x = pcoa.1, y = pcoa.2, color = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 1) +
  scale_color_manual(values = my_colors_site) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x.bottom = element_blank(),  # 移除底部x轴标题（因x轴移到顶部）
    aspect.ratio = 1
  ) +
  scale_x_continuous(position = "top") +  # x轴刻度/标签移到顶部
  labs(x = "PCoA 1", y = "PCoA 2")

# ====== 4. 拼图与对齐 ======
final_plot <- (
  (p_final + p_side) /
    (p_bottom + plot_spacer())
) +
  plot_layout(
    widths = c(1, 0.3), heights = c(1, 0.3),
    guides = "collect"
  ) +
  theme(plot.margin = margin(10, 10, 10, 10))

# 显示图形
print(final_plot)