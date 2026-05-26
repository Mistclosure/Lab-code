# ==============================================================================
# Script 5: Monocytes 亚群比例分析 (Cold_4C vs RT_25C) - X轴为注释版本
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(qs)
library(ggplot2)

print("🚀 开始执行 Script 5: 计算并可视化 Monocyte 亚群比例...")

# 1. 加载 Script 4 保存的单核细胞数据
if(file.exists('pbmc_monocytes_sub-clustered.qs')){
  mono_cells <- qread('pbmc_monocytes_sub-clustered.qs')
} else {
  stop("❌ 未找到 pbmc_monocytes_sub-clustered.qs，请确保已成功运行 Script 4。")
}

# ------------------------------------------------------------------------------
# 2. 提取 Metadata 并计算比例
# ------------------------------------------------------------------------------
print("📊 正在计算 Cold_4C 和 RT_25C 中各 Cluster 占单核细胞的比例...")

# 获取细胞的 metadata
meta_data <- mono_cells@meta.data

# 过滤目标组别，统计各组中各亚群的细胞数，并计算百分比
# 【修改】在 group_by 中加入 mono_annotation，以便将其保留到最终的统计表中
prop_data <- meta_data %>%
  filter(Group %in% c("Cold_4C", "RT_25C")) %>%
  group_by(Group, mono_cluster_id, mono_annotation) %>% 
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Group) %>%
  mutate(proportion = cell_count / sum(cell_count) * 100) # 计算百分比

# 【新增】锁定 mono_annotation 的因子水平，使其强制按照 1~6 的生物学顺序排列，而不是按字母排列
annotation_levels <- c(
  "Activated/Mature Monocytes",
  "Pro-inflammatory classic Monocytes",
  "FKBP5+ classic Monocytes",
  "FKBP5+ non-classic Monocytes",
  "Antigen-Presenting Monocytes",
  "ISG high Monocytes"
)
prop_data$mono_annotation <- factor(prop_data$mono_annotation, levels = annotation_levels)

# 确保输出目录存在
if(!dir.exists("pictures")) dir.create("pictures")

# 定义统一的配色方案 (与 Seurat 默认离散色板对齐，共6个亚群)
# 我们将颜色映射保留给 mono_cluster_id，确保图表颜色和总图严格一致
my_colors <- scales::hue_pal()(6)
names(my_colors) <- paste0("Mono_", 1:6)

# ------------------------------------------------------------------------------
# 3. 绘制并导出图 1: Cold_4C 亚群比例图
# ------------------------------------------------------------------------------
print("🖼️ 正在生成 Cold_4C 亚群比例图...")

prop_cold <- prop_data %>% filter(Group == "Cold_4C")

# 【修改】x = mono_annotation
p_cold <- ggplot(prop_cold, aes(x = mono_annotation, y = proportion, fill = mono_cluster_id)) +
  geom_col(color = "black", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", proportion)), vjust = -0.8, size = 4) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # 顶部留白，防止文字溢出
  theme_classic() +
  labs(title = "Proportion of Monocyte Clusters (Cold_4C)",
       x = "Monocyte Annotation",
       y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"), # 优化长文本的对齐
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 30), # 增加左侧边距防止长标签被截断
        legend.position = "none") # x轴已有标签，隐藏图例

file_name_cold <- file.path("pictures", "Proportion_Monocytes_Cold_4C_Annotated.png")
# 【修改】由于文字较长，略微增加了图片宽度并调整了高度比例，防止文字挤压
ggsave(filename = file_name_cold, plot = p_cold, width = 9, height = 6.5, dpi = 300)

# ------------------------------------------------------------------------------
# 4. 绘制并导出图 2: RT_25C 亚群比例图
# ------------------------------------------------------------------------------
print("🖼️ 正在生成 RT_25C 亚群比例图...")

prop_rt <- prop_data %>% filter(Group == "RT_25C")

# 【修改】x = mono_annotation
p_rt <- ggplot(prop_rt, aes(x = mono_annotation, y = proportion, fill = mono_cluster_id)) +
  geom_col(color = "black", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", proportion)), vjust = -0.8, size = 4) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_classic() +
  labs(title = "Proportion of Monocyte Clusters (RT_25C)",
       x = "Monocyte Annotation",
       y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 30),
        legend.position = "none")

file_name_rt <- file.path("pictures", "Proportion_Monocytes_RT_25C_Annotated.png")
ggsave(filename = file_name_rt, plot = p_rt, width = 9, height = 6.5, dpi = 300)

print(paste("✅ 比例分析完成！带有注释标签的两张组别占比图已保存至", file.path(getwd(), "pictures"), "目录下。"))
# ==============================================================================
# ------------------------------------------------------------------------------
# 5. 绘制并导出图 3: 合并对比柱状图 (Grouped Bar Chart, 整体比例)
# ------------------------------------------------------------------------------
print("🖼️ 正在生成 Cold_4C 与 RT_25C 的并排对比柱状图...")

# 提取用于对比的两组数据 (直接使用脚本前半部分已经算好的 prop_data)
# 确保因子的顺序正确
prop_data$mono_annotation <- factor(prop_data$mono_annotation, levels = annotation_levels)

# 为两个温度组定义颜色 (冷色调代表 Cold，暖色调代表 RT)
group_colors <- c("Cold_4C" = "#3C5488FF", "RT_25C" = "#DC0000FF")

# 绘制并排对比图
p_combined <- ggplot(prop_data, aes(x = mono_annotation, y = proportion, fill = Group)) +
  # position_dodge(width = 0.8) 让两组的柱子并排靠在一起
  geom_col(position = position_dodge(width = 0.8), color = "black", alpha = 0.85, width = 0.7) +
  # 在柱子上添加具体的百分比数值
  geom_text(aes(label = sprintf("%.1f%%", proportion)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # 顶部留白防止文字被裁
  theme_classic() +
  labs(title = "Comparison of Monocyte Subgroups Proportion",
       x = "Monocyte Annotation",
       y = "Percentage of Monocytes (%)",
       fill = "Temperature Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 30),
        legend.position = "top", # 图例放顶部更好看
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 11))

# 保存图片
file_name_combined <- file.path("pictures", "Proportion_Monocytes_Combined_Annotated.png")
ggsave(filename = file_name_combined, plot = p_combined, width = 9, height = 6.5, dpi = 300)

print(paste("✅ 对比柱状图生成完毕！已保存至:", file_name_combined))
# ------------------------------------------------------------------------------
# 6. 绘制并导出图 4: 基于 Fkbp5 表达与否的亚群比例堆叠图 (复刻附图风格)
# ------------------------------------------------------------------------------
print("🖼️ 正在生成基于 Fkbp5 表达的分组堆叠柱状图...")

# 提取 Fkbp5 基因的表达数据
# 注意：小鼠基因通常是首字母大写 (Fkbp5)。如果你的矩阵中全大写，请改成 "FKBP5"
gene_target <- "Fkbp5" 

# 检查基因是否存在于矩阵中
if (gene_target %in% rownames(mono_cells)) {
  
  # 提取该基因的表达量矩阵 (默认提取归一化后的 data 层)
  expr_data <- FetchData(mono_cells, vars = gene_target)
  
  # 在 meta.data 中新增一列，判断表达量是否大于 0
  mono_cells$Fkbp5_status <- ifelse(expr_data[[gene_target]] > 0, "Fkbp5+", "Fkbp5-")
  
  # 将状态转化为因子，控制 Y 轴上下显示的顺序 (这里让 Fkbp5+ 显示在上面，Fkbp5- 在下面)
  mono_cells$Fkbp5_status <- factor(mono_cells$Fkbp5_status, levels = c("Fkbp5-", "Fkbp5+"))
  
  # 重新提取带有新标签的 meta.data 并统计比例
  prop_fkbp5 <- mono_cells@meta.data %>%
    group_by(Fkbp5_status, mono_annotation) %>%
    summarise(cell_count = n(), .groups = 'drop') %>%
    group_by(Fkbp5_status) %>%
    mutate(proportion = cell_count / sum(cell_count) * 100) # 计算百分比 (0-100)
  
  # 固定图例中亚群的顺序 (复用之前的 annotation_levels)
  prop_fkbp5$mono_annotation <- factor(prop_fkbp5$mono_annotation, levels = annotation_levels)
  
  # 为了复刻附图：定义一套区分度高的颜色 (如果你想用之前的 my_colors 也可以直接替换)
  # 这里提供一套类似你附图中的多色调色板
  cluster_colors <- c("#D51F26", "#616599", "#3C8E99", "#44B84A", "#4B8263", "#A8538E")
  names(cluster_colors) <- annotation_levels
  
  # 绘制水平堆叠图
  p_stacked <- ggplot(prop_fkbp5, aes(x = Fkbp5_status, y = proportion, fill = mono_annotation)) +
    # 堆叠柱状图，添加黑色边框线，线条改细一点
    geom_col(color = "black", linewidth = 0.3, width = 0.6) + 
    # 翻转坐标轴，变成水平图 (x变y，y变x)
    coord_flip() + 
    # 填充颜色
    scale_fill_manual(values = cluster_colors) +
    # X轴刻度设为 0-100，不留白
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, by = 20)) + 
    theme_classic() +
    labs(title = "Cluster proportions",
         x = NULL, # 去掉 Y 轴（翻转后的X轴）标题
         y = "Percentage (%)",
         fill = "Cluster") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"), # 标签如 Fkbp5+
      axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
      axis.title.x = element_text(size = 13, face = "bold"),
      legend.position = "right", # 图例放右边
      legend.title = element_blank(), # 去除图例标题
      axis.line.y = element_blank(),  # 去掉左侧多余的竖线
      axis.ticks.y = element_blank(), # 去掉左侧刻度小短线
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  # 保存图片
  file_name_stacked <- file.path("pictures", "Proportion_Monocytes_Fkbp5_Stacked.png")
  ggsave(filename = file_name_stacked, plot = p_stacked, width = 9, height = 4.5, dpi = 300)
  
  print(paste("✅ Fkbp5 堆叠图生成完毕！已保存至:", file_name_stacked))
  
} else {
  print(paste("❌ 错误：在数据中未找到名为", gene_target, "的基因，请检查大小写是否拼写正确（如 'FKBP5' 或 'Fkbp5'）。"))
}