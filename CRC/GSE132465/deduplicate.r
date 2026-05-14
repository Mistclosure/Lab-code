# 1. 读取数据
df <- read.csv("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/files/CRC_Proliferation_Invasion_Metastasis_Genes.csv", check.names = FALSE)

# 2. 对第一列去重，!duplicated 表示保留非重复项（默认保留第一条）
df_unique <- df[!duplicated(df[, 1]), ]

# 3. 保存结果
write.csv(df_unique, "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/files/CRC_Genes_Deduplicated.csv", row.names = FALSE)