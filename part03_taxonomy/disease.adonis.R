myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}
setwd("/share/data1/limin/Airway/02.taxo/")
dist_matrix <- myread("all.taxo.dist.txt")
mydist <- as.dist(dist_matrix)

sample_map = read.table("../00.data/sample.info.txt", header = T,  sep = "\t")
names(sample_map)[1] <- "Sample"
sample_map <- sample_map %>% 
  filter(
    !is.na(Disease),          # 1. 排除真正的 NA 缺失值
    str_squish(Disease) != "" # 2. 去除所有空格后，再判断是否不等于空字符串
  )

sample_map$Disease <- sample_map$Group
rownames(sample_map) <- sample_map$Sample
names(sample_map)[c(3,5)] <- c("Site","Country")
#sample_map <- sample_map[sample_map$Site %in% c("Sputum", "Nasopharynx",  "Throat"), ]
#sample_map <- sample_map[sample_map$Country %in% c("China", "USA", "Germany"), ]
#
xx <- unique(sample_map$Country)
# 提前将原始距离矩阵转为普通矩阵，避免在循环中反复转换
raw_dist_mat <- as.matrix(dist_matrix)
matrix_ids <- trimws(rownames(raw_dist_mat))

# 创建空列表用于存储所有结果
results_list <- list()

# ========== 2. 按国家 (Country) 循环分析疾病状态 (Disease) ==========
target_countries <- unique(sample_map$Country)

for (country in target_countries) {
  
  # --- 步骤 A: 筛选当前国家的样本 ---
  data_samples_raw <- sample_map[sample_map$Country == country, ]
  
  if(nrow(data_samples_raw) == 0) {
    cat(paste0("⚠️ 警告：未找到 Country = '", country, "' 的样本，已跳过。\n"))
    next 
  }
  
  # --- 步骤 B: 匹配索引并提取子矩阵 ---
  valid_samples <- intersect(trimws(as.character(data_samples_raw$Sample)), matrix_ids)
  
  if(length(valid_samples) < 5) { # 增加安全阈值：样本太少无法做 Adonis
    cat(paste0("⚠️ 警告：", country, " 有效样本数过少 (", length(valid_samples), ")，已跳过。\n"))
    next
  }
  
  row_indices <- match(valid_samples, matrix_ids)
  # 【关键】确保元数据和距离矩阵严格对齐
  data_samples <- data_samples_raw[match(valid_samples, data_samples_raw$Sample), ] 
  
  sub_mat <- raw_dist_mat[row_indices, row_indices]
  
  # --- 步骤 C: 清洗矩阵 ---
  sub_mat[is.na(sub_mat)] <- 0
  sub_mat[sub_mat < 0] <- 0
  data_dist_matrix_clean <- as.dist(sub_mat)
  
  # --- 步骤 D: 运行 adonis2 分析 (按国家比较 Disease) ---
  # 【修改点】：去掉了 parallel=4，防止嵌套并行导致内存爆炸；公式改为 ~ Disease
  ado_result <- tryCatch({
    adonis2(data_dist_matrix_clean ~ Disease, 
            data = data_samples, 
            permutations = 999)
  }, error = function(e) {
    cat(paste0("❌ ", country, " Adonis 运行失败: ", e$message, "\n"))
    return(NULL)
  })
  
  if(is.null(ado_result)) next
  
  # --- 步骤 E: 计算校正 R² ---
  r2 <- ado_result$R2[1]
  df_model <- ado_result$Df[1]
  n_subset <- length(valid_samples)
  
  adj_r2_val <- tryCatch({
    RsquareAdj(x = r2, n = n_subset, m = df_model)$adj.r.squared
  }, error = function(e) NA)
  
  ado_result$adj_R2 <- adj_r2_val 
  
  # 将结果存入列表
  results_list[[paste0("Country_", country)]] <- ado_result
  
  cat(paste0("✔️ 成功完成 ", country, " 的分析 (有效样本数: ", n_subset, ", Adj.R2: ", round(adj_r2_val, 4), ")\n"))
}

# ========== 3. 统一保存结果 ==========
save(results_list, file = "taxo_Countries_Disease_Ado_Results.Rdata")

############
results_list <- list() 

target_sites <- unique(sample_map$Site)

for (Site in target_sites) {
  
  # --- 步骤 A: 筛选当前国家的样本 ---
  data_samples_raw <- sample_map[sample_map$Site == Site, ]
  
  if(nrow(data_samples_raw) == 0) {
    cat(paste0("⚠️ 警告：未找到 Site = '", Site, "' 的样本，已跳过。\n"))
    next 
  }
  
  # --- 步骤 B: 匹配索引并提取子矩阵 ---
  valid_samples <- intersect(trimws(as.character(data_samples_raw$Sample)), matrix_ids)
  
  if(length(valid_samples) < 5) { # 增加安全阈值：样本太少无法做 Adonis
    cat(paste0("⚠️ 警告：", Site, " 有效样本数过少 (", length(valid_samples), ")，已跳过。\n"))
    next
  }
  
  row_indices <- match(valid_samples, matrix_ids)
  # 【关键】确保元数据和距离矩阵严格对齐
  data_samples <- data_samples_raw[match(valid_samples, data_samples_raw$Sample), ] 
  
  sub_mat <- raw_dist_mat[row_indices, row_indices]
  
  # --- 步骤 C: 清洗矩阵 ---
  sub_mat[is.na(sub_mat)] <- 0
  sub_mat[sub_mat < 0] <- 0
  data_dist_matrix_clean <- as.dist(sub_mat)
  
  # --- 步骤 D: 运行 adonis2 分析 (按国家比较 Disease) ---
  # 【修改点】：去掉了 parallel=4，防止嵌套并行导致内存爆炸；公式改为 ~ Disease
  ado_result <- tryCatch({
    adonis2(data_dist_matrix_clean ~ Disease, 
            data = data_samples, 
            permutations = 999)
  }, error = function(e) {
    cat(paste0("❌ ", Site, " Adonis 运行失败: ", e$message, "\n"))
    return(NULL)
  })
  
  if(is.null(ado_result)) next
  
  # --- 步骤 E: 计算校正 R² ---
  r2 <- ado_result$R2[1]
  df_model <- ado_result$Df[1]
  n_subset <- length(valid_samples)
  
  adj_r2_val <- tryCatch({
    RsquareAdj(x = r2, n = n_subset, m = df_model)$adj.r.squared
  }, error = function(e) NA)
  
  ado_result$adj_R2 <- adj_r2_val 
  
  # 将结果存入列表
  results_list[[paste0("Site_", Site)]] <- ado_result
  
  cat(paste0("✔️ 成功完成 ", Site, " 的分析 (有效样本数: ", n_subset, ", Adj.R2: ", round(adj_r2_val, 4), ")\n"))
}

# ========== 3. 统一保存结果 ==========
save(results_list, file = "taxo_Site_Disease_Ado_Results.Rdata")

#############
extract_adonis_stats <- function(result_obj, group_name) {
  # result_obj: adonis2 的结果对象
  # group_name: 当前循环的组名（如 "Sputum" 或 "China"）
  
  # 获取自变量名称（即公式左边的变量，通常是 Disease 或 Country）
  factor_name <- rownames(result_obj)[1]
  
  # 提取关键数值
  r2_val <- result_obj$R2[1]
  p_val <- result_obj$`Pr(>F)`[1]
  
  # 返回一个单行数据框
  return(data.frame(
    Group = group_name,       # 分组名称（如 Site 或 Country）
    Factor = factor_name,     # 解释变量（如 Disease）
    R2 = r2_val,              # 解释度
    P_value = p_val           # P值
  ))
}

# --- 2. 批量提取并合并 ---
library(purrr)
load("taxo_Countries_Disease_Ado_Results.Rdata")
plot_data1 <- imap_dfr(results_list, extract_adonis_stats)
plot_data1$Group <- sub(".*_", "", plot_data1$Group)
load("taxo_Site_Disease_Ado_Results.Rdata")
plot_data2 <- imap_dfr(results_list, extract_adonis_stats)
plot_data2$Group <- sub(".*_", "", plot_data2$Group)

plot_data <- plot_data2
#target_groups <- c("China", "USA", "Germany", "Sputum", "Nasopharynx", "Throat")
#plot_data <- plot_data[plot_data$Group %in% target_groups, ]

p <- ggplot(plot_data, aes(x = Group, y = R2, fill = P_value)) +
  
  geom_bar(stat = "identity", color = NA) +
  
  coord_flip() +
  
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_blank(),
    
    legend.title = element_blank() # 去除图例标题，让图例直接显示 p 值范围
  ) +
  
  # 【修改点2】使用连续渐变配色方案代替离散的 brewer 色板
  scale_fill_gradient(
    low = "#a6cee3",   # p 值较小（显著）时的颜色
    high = "#fb9a99",  # p 值较大（不显著）时的颜色
    trans = "log10"    # 💡 强烈建议：p值跨度大时，使用 log10 转换使颜色过渡更自然
  ) +
  
  labs(x = "Group / Site", 
       y = expression(R^2),
       fill = "P-value") # 可选：为图例添加一个清晰的标题

p <- p + 
  geom_text(
    aes(label = paste0("italic(p)==", round(P_value, 4))), # 💡 核心修改：拼接字符串并设置斜体
    hjust = 1,                   # 文字靠右对齐（在柱子内部）
    size = 4,                    # 字体大小
    color = "red",
    parse = TRUE                 # 💡 关键参数：必须设为TRUE，R才会解析并渲染斜体语法
  )

print(p)

