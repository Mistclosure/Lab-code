# ==========================================
# TElocal 结果合并脚本 (优化版)
# ==========================================

# 1. 获取工作路径
wd <- '/mnt/disk1/qiuzerui/downloads/Phf8_GSE212779'
filesname <- 'Phf8_GSE212779_TElocal_locus_counts.csv' # 建议文件名体现出是 locus 水平

setwd(wd)
if (!require("data.table")) install.packages("data.table")
library(data.table)

# 获取所有 .cntTable 文件
# 这里的 pattern 匹配我们之前脚本输出的 "_TElocal.cntTable"
files <- list.files(path = paste0(wd, '/counts'), 
                    pattern = "_TElocal\\.cntTable$", 
                    recursive = TRUE, 
                    full.names = TRUE)

if (length(files) == 0) {
  stop("❌ 未找到任何 TElocal.cntTable 文件，请检查路径和文件名！")
}

# 2. 读取并整合
te_list <- lapply(files, function(f){
  # 使用 fread 读取，TElocal 输出通常没有表头，或者第一列是 ID
  dt <- fread(f, header = FALSE) 
  
  # 【优化点】提取样本名：去掉路径，并去掉 "_TElocal.cntTable" 后缀
  # 这样你的列名就是纯粹的样本 ID，如 "L1MKL2609676-Scr_1_Mixt"
  sample <- gsub(".*配合/(.*)_TElocal\\.cntTable$", "\\1", f)
  # 如果上面的 gsub 没起作用，可以用下面更通用的这行：
  sample <- gsub("_TElocal\\.cntTable$", "", basename(f))
  
  colnames(dt) <- c("RepeatID", sample)
  setkey(dt, RepeatID) # 设置 key 加速后续合并
  return(dt)
})

# 3. 合并为一个大表
# 使用 merge.data.table 效率更高
te_merged <- Reduce(function(x, y) merge(x, y, by="RepeatID", all=TRUE), te_list)

# 将 NA 填为 0（表示该位点在某些样本中没有 read 覆盖）
te_merged[is.na(te_merged)] <- 0 

# 4. 输出文件
write.csv(te_merged, filesname, quote=FALSE, row.names = FALSE)

cat(paste0("✅ 已成功合并 ", length(files), " 个样本\n"))
cat(paste0("💾 结果保存至: ", wd, "/", filesname, "\n"))