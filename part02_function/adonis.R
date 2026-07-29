# 1. 定义需要分析的 Site 列表
target_sites <- c("Sputum", "Nasopharynx",  "Throat")

# 2. 提前将原始距离矩阵转为普通矩阵，避免在循环中反复转换（提升性能）
raw_dist_mat <- as.matrix(dist_matrix)
matrix_ids <- trimws(rownames(raw_dist_mat))

# 3. 开始循环处理每个 Site
for (site in target_sites) {
  
  # --- 步骤 A: 筛选当前 Site 的样本 ---
  data_samples_raw <- sample_map[sample_map$Site == site, ]
  
  # 安全检查：如果该 Site 没有样本，跳过本次循环
  if(nrow(data_samples_raw) == 0) {
    cat(paste0("警告：在 sample_map 中未找到 Site = '", site, "' 的样本，已跳过。\n"))
    next 
  }
  
  # --- 步骤 B: 匹配索引并提取子矩阵 ---
  data_sample_ids <- trimws(as.character(data_samples_raw$Sample))
  valid_samples <- data_sample_ids[data_sample_ids %in% matrix_ids]
  row_indices <- match(valid_samples, matrix_ids)
  row_indices <- row_indices[!is.na(row_indices)]
  
  # 【关键】使用有效样本重新对齐元数据表！确保行数严格一致
  data_samples <- data_samples_raw[match(valid_samples, data_sample_ids), ] 
  
  # 提取子矩阵
  sub_mat <- raw_dist_mat[row_indices, row_indices]
  
  # --- 步骤 C: 清洗矩阵 (处理 NA 和负数) ---
  if(sum(is.na(sub_mat)) > 0) {
    sub_mat[is.na(sub_mat)] <- 0
  }
  sub_mat[sub_mat < 0] <- 0
  data_dist_matrix_clean <- as.dist(sub_mat)
  
  # --- 步骤 D: 运行 adonis2 分析 ---
  # 注意：此时自变量改为了 Country（根据你的需求）
  ado_result <- adonis2(data_dist_matrix_clean ~ Country, 
                        data = data_samples, 
                        permutations = 999, parallel = 4) # 建议先用较小的并行数防崩溃
  
  # --- 步骤 E: 计算校正 R² (使用当前子集的实际样本数 n_subset) ---
  r2 <- ado_result$R2
  df_model <- ado_result$Df 
  n_subset <- length(valid_samples)  # 【关键修正】必须使用当前子集的样本数
  
  adj_r2 <- RsquareAdj(r2, m = df_model, n = n_subset)
  
  # 【修复点】直接将 adj_r2 存入结果中，不要加 $r.squared
  ado_result$adj_R2 <- adj_r2 
  
  # --- 步骤 F: 保存为 RData ---
  save_file_name <- paste0("Site_", site, "_ado.Rdata")
  save(ado_result, file = save_file_name)
  
  cat(paste0("✔️ 成功完成 ", site, " 的分析 (有效样本数: ", n_subset, ")，已保存至 ", save_file_name, "\n"))
}

############
# 1. 定义需要分析的国家列表
target_countries <- c("China", "USA", "Germany")

# 2. 提前将原始距离矩阵转为普通矩阵，避免在循环中反复转换（提升性能）
raw_dist_mat <- as.matrix(dist_matrix)
matrix_ids <- trimws(rownames(raw_dist_mat))

# 3. 开始循环处理每个国家
for (country in target_countries) {
  
  # --- 步骤 A: 筛选当前国家的样本 ---
  data_samples_raw <- sample_map[sample_map$Country == country, ]
  
  # 安全检查：如果该国家没有样本，跳过本次循环
  if(nrow(data_samples_raw) == 0) {
    cat(paste0("警告：在 sample_map 中未找到 Country = '", country, "' 的样本，已跳过。\n"))
    next 
  }
  
  # --- 步骤 B: 匹配索引并提取子矩阵 ---
  data_sample_ids <- trimws(as.character(data_samples_raw$Sample))
  valid_samples <- data_sample_ids[data_sample_ids %in% matrix_ids]
  row_indices <- match(valid_samples, matrix_ids)
  row_indices <- row_indices[!is.na(row_indices)]
  
  # 【关键】使用有效样本重新对齐元数据表！确保行数严格一致
  data_samples <- data_samples_raw[match(valid_samples, data_sample_ids), ] 
  
  # 提取子矩阵
  sub_mat <- raw_dist_mat[row_indices, row_indices]
  
  # --- 步骤 C: 清洗矩阵 (处理 NA 和负数) ---
  if(sum(is.na(sub_mat)) > 0) {
    sub_mat[is.na(sub_mat)] <- 0
  }
  sub_mat[sub_mat < 0] <- 0
  data_dist_matrix_clean <- as.dist(sub_mat)
  
  # --- 步骤 D: 运行 adonis2 分析 ---
  # 比较同一国家内不同部位(Site)的差异
  ado_result <- adonis2(data_dist_matrix_clean ~ Site, 
                        data = data_samples, 
                        permutations = 999, parallel = 4) # 建议先用较小的并行数防崩溃
  
  # --- 步骤 E: 计算校正 R² (使用当前子集的实际样本数 n_subset) ---
  r2 <- ado_result$R2
  df_model <- ado_result$Df 
  n_subset <- length(valid_samples)  # 【关键修正】必须使用当前子集的样本数
  
  adj_r2 <- RsquareAdj(r2, m = df_model, n = n_subset)
  
  # 直接将 adj_r2 存入结果中
  ado_result$adj_R2 <- adj_r2 
  
  # --- 步骤 F: 保存为 RData (【修复点】全部替换为 country 变量) ---
  save_file_name <- paste0("Country_", country, "_ado.Rdata")
  save(ado_result, file = save_file_name)
  
  cat(paste0("✔️ 成功完成 ", country, " 的分析 (有效样本数: ", n_subset, ")，已保存至 ", save_file_name, "\n"))
}
