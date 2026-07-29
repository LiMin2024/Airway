myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}
setwd("/share/data1/limin/Airway/02.taxo/")
data <- myread("taxonomy.profiles.known.norm")

sample_map = read.table("../00.data/sample.info.txt", header = T,  sep = "\t")
names(sample_map)[1] <- "Sample"
sample_map <- sample_map%>%
  select(Sample,site,country)
names(sample_map) <- c("Sample","Site","Country")
rownames(sample_map) <- sample_map$Sample

sample_map <- sample_map[sample_map$Site %in% c("Sputum", "Nasopharynx",  "Throat"), ]
sample_map <- sample_map[sample_map$Country %in% c("China", "USA", "Germany"), ]

##############
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

result <- align_dt_sample(data, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=F, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

taxo <- myread("geneset.taxonomy.known.gz")
taxo2 <- taxo[taxo$V9!="Unclassified",]

bac_taxo <- taxo2[taxo2$V2 == "Bacteria",]
fun_taxo <- taxo2[taxo2$V3 == "k_Fungi",]
vir_taxo <- taxo2[taxo2$V2 == "Viruses",]

bac_data <- data[rownames(data) %in% bac_taxo$V9,]
fun_data <- data[rownames(data) %in% fun_taxo$V9,]
vir_data <- data[rownames(data) %in% vir_taxo$V9,]

alpha_result <- calculate_alpha_by_sample(vir_data, sample_map)

#save(alpha_result,file="vir.alpha_result.RData")
load("04.multi/vir.alpha_result.RData")
data <- alpha_result[,c(1:3)]

#############SampleType
data2 <- merge(sample_map,data,by="Sample")

data3 <- merge(sample_map,data,by="Sample")

#################
# ================== 0. 准备环境与数据 ==================
library(ggplot2)
library(rstatix)      # 用于 wilcox_test 和显著性标注
library(tidyr)        # 用于 pivot_longer 数据重塑
library(dplyr)        # 数据处理
library(patchwork)    # 可选：用于多图拼接

# 假设您的原始数据框为 df，包含列: Country, Site, obs, shannon
# 定义颜色变量 (请根据您的实际色号替换)
my_colors1 <- c("#e31a1c", "#fdbf6f", "#cab2d6")
my_colors2 <- c("#a6cee3", "#b2df8a", "#fb9a99")  

data2$Country <- factor(data2$Country, 
                        levels = c("China", "USA", "Germany"))
data3$Site <- factor(data3$Site, 
                     levels = c("Sputum", "Nasopharynx", "Throat"))
# 定义需要分析的类别
countries <- c("China", "USA", "Germany")
sites <- c("Sputum", "Nasopharynx", "Throat")

theme_pub_base <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 19),
    axis.text.y = element_text(size = 17),
    axis.title = element_text(size = 19),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.margin = margin(0, 0, 0, 0) # 先清零，由 patchwork 统一控制外边距
  )

# ================== 1. 封装通用绘图函数 (修复版) ==================
plot_alpha_diversity_split <- function(data, group_var, index_vars, color_values) {
  
  # --- 1. 数据重塑 ---
  df_long <- data %>%
    select(all_of(c(group_var, index_vars))) %>%
    pivot_longer(cols = all_of(index_vars), names_to = "Index", values_to = "Value")
  
  df_long[[group_var]] <- factor(df_long[[group_var]], levels = unique(df_long[[group_var]]))
  
  # --- 2. 计算统计结果 ---
  stat_results <- df_long %>%
    group_by(Index) %>%
    wilcox_test(as.formula(paste("Value ~", group_var))) %>%
    add_significance("p") %>%
    add_xy_position(x = group_var, dodge = 0, step.increase = 0.1)
  
  # --- 3. 绘制图 A：显著性差异图 ---
  p_sig <- ggplot() +
    stat_pvalue_manual(
      stat_results,
      label = "p.signif",
      tip_length = 0.01,
      hide.ns = TRUE,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ Index, scales = "free_y") +
    labs(y = NULL, x = NULL) +
    scale_fill_manual(values = color_values) +
    theme_pub_base +
    theme(
      # 【关键】彻底隐藏 Y 轴相关元素，防止与下方图错位
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    coord_cartesian(expand = TRUE) # 移除 fixed ratio，避免报错
  
  # --- 4. 绘制图 B：小提琴分布图 ---
  p_dist <- ggplot(df_long, aes(x = .data[[group_var]], y = Value, fill = .data[[group_var]])) +
    geom_violin(trim = FALSE, alpha = 0.8, size = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA, fill = "white", size = 0.5) +
    facet_wrap(~ Index, scales = "free_y") +
    scale_fill_manual(values = color_values) +
    labs(x = group_var, y = "Diversity Index") +
    guides(fill = "none") +
    theme_pub_base +
    theme(
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      # 顶部负边距，让主图向上吸附到显著性图底部
      plot.margin = margin(t = -10, b = 5, l = 5, r = 5)
    )
  
  # --- 5. 拼接并对齐 ---
  # heights = c(1, 4) 意味着下方图的高度约为上方的 4 倍。
  # 由于通常小提琴图的 Y 轴数值跨度远大于显著性线的跨度，
  # 设置 1:4 的比例在视觉上最容易接近“高是宽的 2 倍”的效果。
  combined_plot <- (p_sig / p_dist) + 
    plot_layout(
      heights = c(1, 4),
      ncol = 1,
      guides = 'collect'
    ) + 
    theme(plot.margin = margin(10, 10, 10, 10))
  
  return(combined_plot)
}

# ================== 执行绘图循环 (按 Country 绘制 Site) ==================
plots_by_country <- list()

for (cnt in countries) {
  sub_data <- data2 %>% filter(Country == cnt & Site %in% sites)
  
  if(nrow(sub_data) > 0) {
    p <- plot_alpha_diversity_split(
      data = sub_data,
      group_var = "Site",
      index_vars = c("obs", "shannon"),
      color_values = my_colors2
    ) +
      plot_annotation(
        title = paste("Alpha Diversity in", cnt),
        theme = theme(plot.title = element_text(size = 20, face = "bold"))
      )
    
    plots_by_country[[cnt]] <- p
    print(p)
    
    ggsave(
      filename = paste0("Alpha_Diversity_Site_", cnt, ".pdf"),
      plot = p,
      width = 10,   # 宽度固定为 10
      height = 12,  # 高度设为 12，配合内部 1:4 布局，整体更接近 2:1 的视觉感受
      dpi = 300
    )
  }
}

# ================== 执行绘图循环 (按 Site 绘制 Country) ==================
plots_by_site <- list()

for (st in sites) {
  sub_data <- data3 %>% filter(Site == st & Country %in% countries)
  
  if(nrow(sub_data) > 0) {
    p <- plot_alpha_diversity_split(
      data = sub_data,
      group_var = "Country",
      index_vars = c("obs", "shannon"),
      color_values = my_colors1
    ) +
      plot_annotation(
        title = paste("Alpha Diversity in", st),
        theme = theme(plot.title = element_text(size = 20, face = "bold"))
      )
    
    plots_by_site[[st]] <- p
    print(p)
    
    ggsave(
      filename = paste0("Alpha_Diversity_Country_", st, ".pdf"),
      plot = p,
      width = 10,
      height = 12,
      dpi = 300
    )
  }
}

