setwd("/share/data1/limin/Airway/02.taxo/")
getwd()
rm(list=ls())
load("04.multi/vir_Country.top_species_result.RData")

library(Maaslin2)
library(permute)
library(lattice)
library(stringr)
library(tidyr)
library(dplyr)
library(tibble)

myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

result <- read.table("03.data/all.Country.data", sep="\t", header=T, check.names=F)

result <- result[,-c(2:10)]
label <- result[,1:2]
result <- result[,-2]
names(result) <- c("species","superkingdom","phylum","class","order","family","genus")

####################tree

taxonomy_df <- result

char_cols <- sapply(taxonomy_df, is.character)

#options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#install.packages("data.tree")
library(data.tree)
library(ape)

# 构建路径列表
tree_nodes <- list()
for (i in 1:nrow(taxonomy_df)) {
  row <- taxonomy_df[i, ]
  
  path <- paste0(
    row$superkingdom, "/", 
    row$phylum, "/",
    row$class, "/",
    row$order, "/",
    row$family, "/",
    row$genus, "/",
    row$species
  )
  
  node <- unlist(strsplit(path, "/"))
  tree_nodes[[row$species]] <- as.list(node)
}

# 合并为路径字符串
paths <- lapply(tree_nodes, function(x) paste(x, collapse = "/"))
path_df <- data.frame(path = unlist(paths), stringsAsFactors = FALSE)

clean_path_df <- data.frame(
  path = path_df$path,
  stringsAsFactors = FALSE
)

# 提取叶子名（最后一个 / 后面的内容）
clean_path_df$id <- sub(".*/", "", clean_path_df$path)

# 定义函数：将路径添加到树中
add_path_to_tree <- function(tree, parts, leaf_name) {
  if (length(parts) == 1) {
    tree[[leaf_name]] <- list()
  } else {
    current <- parts[1]
    rest <- parts[-1]
    
    if (is.null(tree[[current]])) {
      tree[[current]] <- list()
    }
    
    # 递归处理子路径，并将结果正确赋值回当前节点
    tree[[current]] <- add_path_to_tree(tree[[current]], rest, leaf_name)
  }
  return(tree)
}

# 定义函数：将树结构转为 Newick 格式
tree_to_newick <- function(tree) {
  if (is.null(tree) || length(tree) == 0) return("")
  
  children <- names(tree)
  child_strings <- sapply(children, function(child) {
    child_tree <- tree[[child]]
    if (length(child_tree) == 0) {
      return(child)
    } else {
      return(paste0("(", tree_to_newick(child_tree), ")", child))
    }
  })
  
  if (length(child_strings) == 0) return("")
  
  return(paste0("(", paste(child_strings, collapse = ","), ")"))
}

# 初始化空树
tree <- list()

# 遍历所有行，添加路径到树中
for (i in seq_len(nrow(clean_path_df))) {
  parts <- unlist(strsplit(clean_path_df$path[i], "/"))
  leaf_name <- clean_path_df$id[i]
  tree <- add_path_to_tree(tree, parts, leaf_name)
}

# 转换为 Newick 格式
newick_str <- paste0(tree_to_newick(tree), ";")

#writeLines(newick_str, "Country.newick")

###################分支颜色
species_to_phylum <- taxonomy_df[, c("species", "phylum")]
species_to_phylum <- as.data.frame(species_to_phylum)


# === 自定义颜色列表 ===
total_color1 = c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                 "#6a3d9a", "#ffff99", "#b15928","#8dd3c7",
                 "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
                 "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
                 "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4",
                 "#b2df8a", "#33a02c", "#fb9a99","#ed1299",
                 "#09f9f5","#246b93","#cc8e12","#d561dd",
                 "#c93f00","#ddd53e","#4aef7b","#e86502",
                 "#9ed84e","#39ba30","#6ad157","#8249aa",
                 "#99db27","#e07233","#ff523f","#ce2523",
                 "#f7aa5d","#cebb10","#03827f","#931635",
                 "#373bbf","#a1ce4c","#ef3bb6","#d66551",
                 "#1a918f","#ff66fc","#2927c4","#7149af",
                 "#57e559","#8e3af4","#f9a270","#22547f",
                 "#db5e92","#edd05e","#6f25e8","#0dbc21",
                 "#280f7a","#6373ed","#5b910f","#7b34c1",
                 "#0cf29a","#d80fc1","#dd27ce","#07a301",
                 "#167275","#391c82","#2baeb5","#925bea","#63ff4f")

# === 提取唯一 phylum ===
unique_phyla <- unique(species_to_phylum$phylum)
n_phyla <- length(unique_phyla)
n_colors <- length(total_color1)

# 检查颜色是否足够
if (n_colors < n_phyla) {
  stop("错误：颜色数量不足！有 ", n_phyla, " 个门，但只有 ", n_colors, " 个颜色。")
}

# 使用你的自定义颜色（按顺序或随机）
assigned_colors <- total_color1[1:n_phyla]  # 按顺序取前 n_phyla 个颜色
# 或者随机打乱：assigned_colors <- sample(total_color1, n_phyla)

# 正确创建 phylum_to_color 数据框（不要覆盖！）
phylum_to_color <- data.frame(
  phylum = unique_phyla,
  color = assigned_colors,
  stringsAsFactors = FALSE
)

# === 现在可以安全合并 ===
data <- left_join(species_to_phylum, phylum_to_color, by = "phylum")

# 再与 clean_path_df 合并
merge_data <- merge(clean_path_df, data, by.x = "id", by.y = "species")
merge_data <- merge_data[3:4]

df <- cbind(merge_data[1], clade = "clade", merge_data[2:ncol(merge_data)])
df <- cbind(df[1:3], normal = "normal")
df <- cbind(df[1:4], size = "3")
length(unique(df$phylum))
table(df$phylum)
# 打开连接（写入文件）
file_conn <- file("branch_color.txt", "w")

# 写入前三行
writeLines("TREE_COLORS", file_conn)
writeLines("SEPARATOR TAB", file_conn)  
writeLines("DATA", file_conn)

# 写入数据（用制表符分隔）
write.table(df, file_conn, 
            sep = "\t",             # 使用制表符分隔
            row.names = FALSE,      # 不写行名
            col.names = FALSE,      # 不重复写列名（因为前面已有 DATA）
            quote = FALSE)          # 不给字符串加引号

# 关闭连接
close(file_conn)

############
color_df <- read.table("branch_color.legend.txt", header = FALSE, sep = "\t", comment.char = "")
color_df <- color_df[,c(1,3)]
color_df <- color_df[!duplicated(color_df), ]

# 创建命名向量
my_colors <- setNames(color_df$V3, color_df$V1)

# 绘制纯图例
ggplot(color_df, aes(x = 1, y = V1, color = V1)) + 
  geom_point() + # 这里用固定的 x=1 和 V1 作为 y 轴来触发图例生成
  scale_color_manual(values = my_colors) +
  labs(color = "Phylum") +
  theme(
    axis.title = element_blank(),   # 隐藏坐标轴标题
    axis.text.y = element_blank(),  # 隐藏 Y 轴文字
    axis.ticks = element_blank(),   # 隐藏刻度线
    panel.background = element_blank() # 隐藏背景网格
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) # 调大色块大小

################
tax <- label

file_name <- "labels.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入注释和头部信息
writeLines(
  c(
    "LABELS",
    "#use this template to change the leaf labels, or define/change the internal node names (displayed in mouseover popups)",
    "#lines starting with a hash are comments and ignored during parsing",
    "",
    "#=================================================================#",
    "#          MANDATORY SETTINGS              #",
    "#=================================================================#",
    "#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file (except in the SEPARATOR line, which uses space).",
    "SEPARATOR TAB",
    "#SEPARATOR SPACE",
    "#SEPARATOR COMMA",
    "",
    "#=================================================================#",
    "#    Actual data follows after the \"DATA\" keyword       #",
    "#=================================================================#",
    "DATA"
  ),
  con = file_conn
)

# 5. 写入数据部分（用制表符分隔）
write.table(tax, file_conn,
            sep = "\t",           # 使用 TAB 分隔
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

close(file_conn)

##############abun
mean_abundance_wide = read.table("Sputum_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(taxonomy_df,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,-c(2:7)]

merge_data2 <- column_to_rownames(merge_data,var="species")
row_scaled <- t(apply(merge_data2, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

row_scaled_df <- rownames_to_column(as.data.frame(row_scaled), var = "name")

file_name <- "Sputum_abun_scale.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入头部信息
writeLines(
  c(
    "DATASET_HEATMAP",
    "SEPARATOR SPACE",
    "DATASET_LABEL Effect Size Heatmap",
    "# 渐变颜色设置",
    "COLOR_MIN lightblue",
    "COLOR_MAX #1f78b4",
    "USE_MID_COLOR 1",
    "COLOR_MID #FFFFFF",
    "STRIP_WIDTH 25",
    "MARGIN 0",
    "SHOW_TREE 1",
    "FIELD_LABELS China Germany USA ",
    "DATA"
  ),
  con = file_conn
)
#Sputum Nasopharynx Throat
#China Germany USA
write.table(row_scaled_df, file_conn,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. 关闭连接
close(file_conn)

#############
mean_abundance_wide = read.table("Nasopharynx_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(taxonomy_df,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,-c(2:7)]

merge_data2 <- column_to_rownames(merge_data,var="species")
row_scaled <- t(apply(merge_data2, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

row_scaled_df <- rownames_to_column(as.data.frame(row_scaled), var = "name")

file_name <- "Nasopharynx_abun_scale.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入头部信息
writeLines(
  c(
    "DATASET_HEATMAP",
    "SEPARATOR SPACE",
    "DATASET_LABEL Effect Size Heatmap",
    "# 渐变颜色设置",
    "COLOR_MIN #b2df8a",
    "COLOR_MAX #33a02c",
    "USE_MID_COLOR 1",
    "COLOR_MID #FFFFFF",
    "STRIP_WIDTH 25",
    "MARGIN 0",
    "SHOW_TREE 1",
    "FIELD_LABELS China Germany USA ",
    "DATA"
  ),
  con = file_conn
)
#Sputum Nasopharynx Throat
#China Germany USA
write.table(row_scaled_df, file_conn,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. 关闭连接
close(file_conn)

############
mean_abundance_wide = read.table("Throat_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(taxonomy_df,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,-c(2:7)]

merge_data2 <- column_to_rownames(merge_data,var="species")
row_scaled <- t(apply(merge_data2, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

row_scaled_df <- rownames_to_column(as.data.frame(row_scaled), var = "name")

file_name <- "Throat_abun_scale.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入头部信息
writeLines(
  c(
    "DATASET_HEATMAP",
    "SEPARATOR SPACE",
    "DATASET_LABEL Effect Size Heatmap",
    "# 渐变颜色设置",
    "COLOR_MIN #fcd5d4",
    "COLOR_MAX #fb9a99",
    "USE_MID_COLOR 1",
    "COLOR_MID #FFFFFF",
    "STRIP_WIDTH 25",
    "MARGIN 0",
    "SHOW_TREE 1",
    "FIELD_LABELS China Germany USA ",
    "DATA"
  ),
  con = file_conn
)
#Sputum Nasopharynx Throat
#China Germany USA
write.table(row_scaled_df, file_conn,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. 关闭连接
close(file_conn)

my_colors1 <- c("lightblue", "#1f78b4", "#b2df8a", "#33a02c","#fcd5d4", "#fb9a99")  

#############
mean_abundance_wide = read.table("China_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(taxonomy_df,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,-c(2:7)]

merge_data2 <- column_to_rownames(merge_data,var="species")
row_scaled <- t(apply(merge_data2, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

row_scaled_df <- rownames_to_column(as.data.frame(row_scaled), var = "name")

file_name <- "China_abun_scale.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入头部信息
writeLines(
  c(
    "DATASET_HEATMAP",
    "SEPARATOR SPACE",
    "DATASET_LABEL Effect Size Heatmap",
    "# 渐变颜色设置",
    "COLOR_MIN pink",
    "COLOR_MAX #e31a1c",
    "USE_MID_COLOR 1",
    "COLOR_MID #FFFFFF",
    "STRIP_WIDTH 25",
    "MARGIN 0",
    "SHOW_TREE 1",
    "FIELD_LABELS Nasopharynx Sputum Throat",
    "DATA"
  ),
  con = file_conn
)
#Sputum Nasopharynx Throat
#China Germany USA
write.table(row_scaled_df, file_conn,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. 关闭连接
close(file_conn)

#############
mean_abundance_wide = read.table("USA_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(taxonomy_df,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,-c(2:7)]

merge_data2 <- column_to_rownames(merge_data,var="species")
row_scaled <- t(apply(merge_data2, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

row_scaled_df <- rownames_to_column(as.data.frame(row_scaled), var = "name")

file_name <- "USA_abun_scale.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入头部信息
writeLines(
  c(
    "DATASET_HEATMAP",
    "SEPARATOR SPACE",
    "DATASET_LABEL Effect Size Heatmap",
    "# 渐变颜色设置",
    "COLOR_MIN ##fee0b7",
    "COLOR_MAX #fdbf6f",
    "USE_MID_COLOR 1",
    "COLOR_MID #FFFFFF",
    "STRIP_WIDTH 25",
    "MARGIN 0",
    "SHOW_TREE 1",
    "FIELD_LABELS Nasopharynx Sputum Throat",
    "DATA"
  ),
  con = file_conn
)
#Sputum Nasopharynx Throat
#China Germany USA
write.table(row_scaled_df, file_conn,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. 关闭连接
close(file_conn)
my_colors2 <- c("pink","#e31a1c", "#fee0b7","#fdbf6f", "#cab2d6","#6a3d9a")

############
mean_abundance_wide = read.table("Germany_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(taxonomy_df,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,-c(2:7)]

merge_data2 <- column_to_rownames(merge_data,var="species")
row_scaled <- t(apply(merge_data2, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

row_scaled_df <- rownames_to_column(as.data.frame(row_scaled), var = "name")

file_name <- "Germany_abun_scale.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入头部信息
writeLines(
  c(
    "DATASET_HEATMAP",
    "SEPARATOR SPACE",
    "DATASET_LABEL Effect Size Heatmap",
    "# 渐变颜色设置",
    "COLOR_MIN #cab2d6",
    "COLOR_MAX #6a3d9a",
    "USE_MID_COLOR 1",
    "COLOR_MID #FFFFFF",
    "STRIP_WIDTH 25",
    "MARGIN 0",
    "SHOW_TREE 1",
    "FIELD_LABELS Nasopharynx Sputum Throat",
    "DATA"
  ),
  con = file_conn
)

write.table(row_scaled_df, file_conn,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. 关闭连接
close(file_conn)
