# ==============================================================================
# utils_paths.R: 路径和配置管理工具
# ==============================================================================
# 提供项目路径获取、配置加载、样本映射读取等基础功能。
# 所有脚本通过此模块获取路径，避免硬编码 setwd()。
# ==============================================================================

#' 加载项目配置文件
#'
#' @param config_path 配置文件路径，默认为项目 config/config.yaml
#' @return 配置列表
load_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    # 尝试从 rstudioapi 获取路径
    tryCatch({
      script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
      config_path <- file.path(dirname(script_dir), "config", "config.yaml")
    }, error = function(e) {
      # 命令行模式：使用 COLDMOUSE_INTRONS_HOME 环境变量或默认路径
      project_home <- Sys.getenv("COLDMOUSE_INTRONS_HOME",
                                 unset = "/home/zerui/code/coldmouse_introns")
      config_path <<- file.path(project_home, "config", "config.yaml")
    })
  }
  if (!file.exists(config_path)) {
    stop("配置文件未找到: ", config_path)
  }
  yaml::read_yaml(config_path)
}

#' 读取样本映射表
#'
#' @param sample_map_path 样本映射文件路径
#' @return data.frame
read_sample_map <- function(sample_map_path = NULL) {
  if (is.null(sample_map_path)) {
    sample_map_path <- file.path(getwd(), "config", "sample_map.csv")
  }
  if (!file.exists(sample_map_path)) {
    stop("样本映射文件未找到: ", sample_map_path)
  }
  read.csv(sample_map_path, stringsAsFactors = FALSE)
}

#' 获取项目输出目录
#'
#' @param config 配置列表
#' @param subdir 子目录名 (objects/files/plots/logs)
#' @return 完整路径
get_output_dir <- function(config, subdir = NULL) {
  base_dir <- config$paths$output_dir
  if (!is.null(subdir)) {
    dir_path <- file.path(base_dir, subdir)
  } else {
    dir_path <- base_dir
  }
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  dir_path
}

#' 获取结果目录
#'
#' @param config 配置列表
#' @param subdir 子目录名 (objects/files/plots)
#' @return 完整路径
get_results_dir <- function(config, subdir = NULL) {
  base_dir <- config$paths$results_dir
  if (!is.null(subdir)) {
    dir_path <- file.path(base_dir, subdir)
  } else {
    dir_path <- base_dir
  }
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  dir_path
}

#' 设置线程和内存策略
#'
#' @param config 配置列表
setup_performance <- function(config) {
  max_size <- config$performance$future_globals_maxSize
  if (!is.null(max_size)) {
    # 解析为字节数
    if (is.character(max_size)) {
      multiplier <- 1
      size_str <- max_size
      if (grepl("G$", size_str, ignore.case = TRUE)) {
        multiplier <- 1024^3
        size_str <- gsub("G$", "", size_str, ignore.case = TRUE)
      } else if (grepl("M$", size_str, ignore.case = TRUE)) {
        multiplier <- 1024^2
        size_str <- gsub("M$", "", size_str, ignore.case = TRUE)
      }
      max_size <- as.numeric(size_str) * multiplier
    }
    options(future.globals.maxSize = max_size)
    message("设置 future.globals.maxSize = ", max_size)
  }
  n_cores <- config$performance$n_cores
  if (!is.null(n_cores)) {
    n_cores <- min(n_cores, parallel::detectCores())
    message("使用 ", n_cores, " 个核心")
  }
}
