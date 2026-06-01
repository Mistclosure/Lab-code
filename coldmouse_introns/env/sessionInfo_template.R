# ==============================================================================
# sessionInfo_template.R: 记录 R 会话信息
# ==============================================================================
# 运行此脚本记录当前 R 会话的包版本信息。
# 输出保存到 env/cmintro_detected_packages.txt
# ==============================================================================

# 记录 sessionInfo
session_info <- sessionInfo()

# 保存到文件
output_file <- file.path(getwd(), "cmintro_detected_packages.txt")
sink(output_file)
cat("=== Session Info ===\n")
print(session_info)
cat("\n=== Loaded Packages ===\n")
loaded <- session_info$loadedOnly
for (pkg_name in names(loaded)) {
  cat(paste0(pkg_name, ": ", loaded[[pkg_name]]$Version, "\n"))
}
sink()

message("Session info 已保存到: ", output_file)
