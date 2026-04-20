# 导入所需的包，如果没有安装 stringr 请先运行 install.packages("stringr")
library(stringr)
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE231559/GSE231559_RAW')
# 1. 根据附图6定义样本 GSM 编号到组别名称的映射关系
sample_mapping <- c(
  "GSM7290760" = "L1N",
  "GSM7290761" = "L1T",
  "GSM7290762" = "C1N",
  "GSM7290763" = "C1T",
  "GSM7290764" = "L2N",
  "GSM7290765" = "L3N",
  "GSM7290766" = "L4N",
  "GSM7290767" = "L4T",
  "GSM7290768" = "C2N",
  "GSM7290769" = "C2T",
  "GSM7290770" = "L5N",
  "GSM7290771" = "C3N",
  "GSM7290772" = "C3T",
  "GSM7290773" = "C4T",
  "GSM7290774" = "C5T",
  "GSM7290775" = "L6T",
  "GSM7290776" = "L7N",
  "GSM7290777" = "C6T",
  "GSM7290778" = "L8T1",
  "GSM7290779" = "L8T2",
  "GSM7290780" = "L9N",
  "GSM7290781" = "L9T",
  "GSM7290782" = "L10T",
  "GSM7290783" = "L11T",
  "GSM7290784" = "L11N",
  "GSM7290785" = "L12T"
)

# 2. 获取当前目录下所有的原始单细胞文件
# 假设您的文件后缀都包含 .gz
all_files <- list.files(pattern = "^GSM.*\\.gz$")

if(length(all_files) == 0) {
  stop("当前目录下未找到匹配的 GSM 文件，请检查工作目录 (setwd) 是否正确！")
}

# 3. 遍历文件，自动创建文件夹、重命名并移动
for (file in all_files) {
  # 利用正则表达式提取文件开头部分的 GSM 号 (例如 "GSM7290760")
  gsm_id <- str_extract(file, "^GSM[0-9]+")
  
  if (!is.na(gsm_id) && gsm_id %in% names(sample_mapping)) {
    
    # 获取映射的目标文件夹名称 (例如 "L1N")
    folder_name <- sample_mapping[[gsm_id]]
    
    # 如果文件夹不存在，则新建该文件夹
    if (!dir.exists(folder_name)) {
      dir.create(folder_name)
    }
    
    # 判断原始文件的类别，并赋予 10x 读取所需的标准名称
    if (grepl("barcodes\\.tsv\\.gz$", file)) {
      new_name <- "barcodes.tsv.gz"
    } else if (grepl("features\\.tsv\\.gz$", file) || grepl("genes\\.tsv\\.gz$", file)) {
      new_name <- "features.tsv.gz"
    } else if (grepl("matrix\\.mtx\\.gz$", file)) {
      new_name <- "matrix.mtx.gz"
    } else {
      # 如果有其他无关文件，则跳过
      next 
    }
    
    # 构建新文件的完整相对路径
    new_path <- file.path(folder_name, new_name)
    
    # 将文件复制并重命名到新文件夹中
    # 注意：这里使用的是 file.copy，原始文件仍会保留，以防出错。
    # 如果您确信无误且想节省硬盘空间，可以将其改为 file.rename(file, new_path)
    success <- file.copy(from = file, to = new_path, overwrite = TRUE)
    
    if(success){
      cat(sprintf("成功: [%s] 已转换为 -> %s\n", file, new_path))
    } else {
      cat(sprintf("失败: 无法处理文件 [%s]\n", file))
    }
  }
}

cat("所有文件分类和重命名处理完毕！\n")