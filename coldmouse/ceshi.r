# ==============================================================================
# 单独检查：针对 sc_by_tissue[[1]] 的 SingleR 标注验证
# ==============================================================================

# 1. 锁定第一个组织对象
obj <- sc_by_tissue[[1]]
tissue_name <- names(sc_by_tissue)[1]
print(paste("🔍 正在检查组织:", tissue_name))

# 2. 【关键】强制汇合所有样本的 Layer (Seurat v5 必需)
# 这一步如果不做，GetAssayData 提取的矩阵可能是不完整的
obj <- JoinLayers(obj)

# 3. 提取归一化后的矩阵
# 确保 layer = "data" 提取的是 LogNormalized 后的数据
expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

# 4. 运行 SingleR
# 提示：请确保环境中心已经运行过 mouse_ref <- celldex::MouseRNAseqData()
print(">>> 正在运行 SingleR 标注 (针对单组织测试)...")
pred_res <- SingleR(test = expr_mat, 
                    ref = mouse_ref, 
                    labels = mouse_ref$label.main)

# 5. 将结果写入 Metadata
# 直接赋值，顺序与 obj 的列名 (Barcode) 一一对应
obj$SingleR.labels <- pred_res$labels

# ==============================================================================
# 结果检查区域
# ==============================================================================

# A. 检查 Metadata 是否成功写入
print("📊 1. Metadata 前 6 行预览 (检查 SingleR.labels 列)：")
print(head(obj@meta.data))

# B. 查看标注出的细胞类型统计 (最直观判断是否有结果)
print("🧬 2. 所有识别出的细胞类型分布统计：")
tag_table <- table(obj$SingleR.labels)
print(tag_table)

# C. 验证髓系细胞过滤逻辑
myeloid_keywords <- c("Macrophages", "Monocytes", "Granulocytes", "Dendritic cells", "Microglia", "Neutrophils")
is_myeloid <- obj$SingleR.labels %in% myeloid_keywords
is_myeloid[is.na(is_myeloid)] <- FALSE

print(paste("✅ 3. 在该组织中命中髓系关键词的细胞数:", sum(is_myeloid)))

# 如果细胞数 > 0，可以尝试看一眼具体的 Barcode 匹配情况
if(sum(is_myeloid) > 0){
  print("📌 命中髓系关键词的前 5 个标签示例：")
  print(head(obj$SingleR.labels[is_myeloid]))
} else {
  warning("❌ 警告：未在该组织中识别到任何髓系细胞，请核对关键字名称与 SingleR 标注结果是否匹配。")
}