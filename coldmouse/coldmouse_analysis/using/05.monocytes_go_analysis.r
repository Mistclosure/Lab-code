# ==============================================================================
# Script 5: Monocytes 亚群 GO 通路富集分析
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# 加载必要的 R 包
library(Seurat)
library(tidyverse)
library(qs)
library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠基因注释数据库

# ------------------------------------------------------------------------------
# 1. 读取单核细胞精细聚类对象
# ------------------------------------------------------------------------------
print("🚀 步骤1: 正在加载 Monocytes 子集对象...")
mono_cells <- qread("pbmc_monocytes_sub-clustered.qs")

# 确保以我们自定义的亚群编号 (Mono_1 到 Mono_6) 作为当前分组标准
Idents(mono_cells) <- "mono_cluster_id"

# ------------------------------------------------------------------------------
# 2. 提取各个亚群的 Marker 基因
# ------------------------------------------------------------------------------
print("🚀 步骤2: 正在计算各亚群的特异性 Marker 基因...")
# 注意：GO富集通常只需要上调的基因 (only.pos = TRUE)
markers <- FindAllMarkers(mono_cells, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25,
                          verbose = FALSE)

# 过滤出显著的 Marker 基因 (校正后 p 值 < 0.05)
sig_markers <- markers %>% filter(p_val_adj < 0.05)

# ------------------------------------------------------------------------------
# 3. 逐亚群运行 GO (Biological Process) 富集分析
# ------------------------------------------------------------------------------
print("🚀 步骤3: 正在执行 GO 生物学过程 (BP) 富集分析...")

# 创建一个空列表用于存储每个亚群的富集结果
go_results_list <- list()

# 获取所有唯一的亚群名称 (Mono_1, Mono_2...)
clusters <- sort(unique(sig_markers$cluster))

for (cluster_name in clusters) {
  print(paste("   正在分析亚群:", cluster_name))
  
  # 提取当前亚群的显著基因名 (SYMBOL)
  genes <- sig_markers %>% 
    filter(cluster == cluster_name) %>% 
    pull(gene)
  
  # 如果该亚群没有足够的显著基因，则跳过
  if (length(genes) < 5) {
    print(paste("   ⚠️ 亚群", cluster_name, "的显著基因少于5个，跳过富集。"))
    next
  }
  
  # 运行 GO 富集分析
  go_enrich <- enrichGO(gene          = genes,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = 'SYMBOL', # 直接使用基因名称匹配
                        ont           = "BP",     # BP: Biological Process (生物学过程)
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2,
                        readable      = FALSE)    # 因为 keyType 已经是 SYMBOL，无需转换
  
  # 如果结果不为空，提取数据框并加上 Cluster 标签
  if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
    # 过滤掉不显著的通路
    go_res_df <- go_enrich@result %>% filter(p.adjust < 0.05)
    
    if (nrow(go_res_df) > 0) {
      # 添加亚群名称列，放在第一列方便后续查看
      go_res_df <- go_res_df %>% 
        mutate(Cluster = cluster_name) %>%
        dplyr::select(Cluster, everything())
      
      go_results_list[[cluster_name]] <- go_res_df
    }
  }
}

# ------------------------------------------------------------------------------
# 4. 合并并导出结果 (已新增注释列映射)
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在整合结果并导出为 CSV...")

# 将所有亚群的结果列表按行合并成一个大的数据框
if (length(go_results_list) > 0) {
  all_go_results <- bind_rows(go_results_list)
  
  # 【重点修改】提取 Seurat 对象中的 Cluster 与 Annotation 的映射关系，并去重
  anno_mapping <- unique(mono_cells@meta.data[, c("mono_cluster_id", "mono_annotation")])
  
  # 【重点修改】将注释信息合并进去，并调整列顺序 (Cluster 第一列，mono_annotation 第二列，其余靠后)
  all_go_results <- all_go_results %>%
    left_join(anno_mapping, by = c("Cluster" = "mono_cluster_id")) %>%
    dplyr::select(Cluster, mono_annotation, everything())
  
  # 导出路径
  output_file <- file.path("files", "GO_BP_Enrichment_Monocytes_Subgroups.csv")
  
  # 确保 files 文件夹存在
  if (!dir.exists("files")) dir.create("files")
  
  write.csv(all_go_results, file = output_file, row.names = FALSE)
  print(paste("✅ 所有亚群 GO 富集分析完成！结果已保存至:", output_file))
  
} else {
  print("⚠️ 警告：没有找到任何显著富集的 GO 通路。")
}