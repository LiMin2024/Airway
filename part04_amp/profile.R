setwd("/share/data1/limin/Airway/03.AMP/")
myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

###############
data <- myread("profs.amp")
amp_profile <- data
taxo <- myread("../02.taxo/geneset.taxonomy.known.gz")
taxo2 <- rownames_to_column(taxo,var="gene")
taxo2 <- taxo2[,-3]
names(taxo2) <- c("gene","domain","phylum","class","order","family","genus","species")

amp <- rownames_to_column(data,var="gene")
amp <- as.data.frame(amp$gene)
names(amp) <- "gene"

merge <- left_join(amp,taxo2,by="gene")
all_na <- apply(merge[, -1], 1, function(x) all(is.na(x)))
merge[all_na, -1] <- "Unclassified"

#write.table(merge,file = "amp_taxon.txt", sep = "\t", row.names = F, quote = FALSE)

##################
data <- merge[,c(1,7)]
data[, 1] <- "AMP"
data <- data %>%
  count(genus, sort = TRUE) 
names(data) <- c("genus","count")

data_grouped <- data %>%
  mutate(status = ifelse(genus %in% c("Unknown", "Unclassified"), "Unknown", "Known")) %>%
  group_by(status) %>%
  summarise(total = sum(count), .groups = 'drop')

pie(data_grouped$total, labels = data_grouped$status, main = "Known vs Unknown", col = c("#e31a1c", "lightgray"))

##################
colors <- c(
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
  "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#ed1299",
  "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e", "#4aef7b", "#e86502",
  "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f", "#ce2523",
  "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
  "#1a918f", "#ff66fc", "#2927c4", "#7149af", "#57e559", "#8e3af4", "#f9a270", "#22547f",
  "#db5e92", "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f", "#7b34c1",
  "#0cf29a", "#d80fc1", "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5", "#925bea", "#63ff4f"
)
known_data <- data %>%
  filter(!genus %in% c("Unknown", "Unclassified")) %>%
  arrange(desc(count))

top_known <- known_data %>% slice(1:15)
others_sum <- sum(known_data$count) - sum(top_known$count)

top_known <- top_known %>% add_row(genus = "Others", count = others_sum)

top_known <- top_known %>%
  arrange(desc(count)) %>%   # 按丰度降序
  mutate(genus = factor(genus, levels = genus))  

ggplot(top_known, aes(x="", y=count, fill=genus)) +
  geom_bar(stat="identity", width=1, color="white") + # 使用白色边框使分割更清晰
  coord_polar("y", start=0) +
  theme_void() +
  labs(title="Top 15 Known Species vs Others", fill="Species") +
  scale_fill_manual(values = colors[1:nrow(top_known)]) # 从您提供的颜色列表中选择颜色
