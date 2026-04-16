# 1. 安装通讯包
if (!require("reticulate", quietly = TRUE)) install.packages("reticulate")

# 2. 自动安装 Python 依赖
# 这会默认安装到 reticulate 的虚拟环境中，或者你当前激活的 conda 环境
reticulate::py_install(c("leidenalg", "igraph", "pandas"))

# 3. 安装 R 端的轻量级封装（可选，但推荐）
if (!require("leiden", quietly = TRUE)) install.packages("leiden")