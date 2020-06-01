
rm(list = ls())

source("./project_support.r")

tic.clearlog()

dir_init("./temp")

tic("fit all models")
datasets_available <- list.files("./data", full.names = TRUE)
analyses_to_run <- datasets_available
analyses_to_run <- gsub("\\./data", "./site_scripts", analyses_to_run)
analyses_to_run <- gsub("\\.csv", ".r", analyses_to_run)
analyses_to_run <- as.list(analyses_to_run)
x <- mclapply(1:length(analyses_to_run), function(z) source(analyses_to_run[[z]]), mc.cores = n_cores)
toc(log = TRUE)

tic("prep output table")
available_reports <- list.files("./temp", pattern = "\\.json$", full.names = TRUE)
available_reports %>% purrr::map(read_json, simplifyVector = TRUE) %>%
  dplyr::bind_rows() %>% as.data.frame() -> res
toc(log = TRUE)

write.csv(res, "./temp/table2.csv", row.names = FALSE)

tic.log(format = TRUE)
msg_log <- unlist(tic.log())
msg_log <- paste0("- ", msg_log)
header <- c(
  "project: wealth transmission model-fitting", 
  "events:")
msg_log <- c(header, msg_log)
writeLines(msg_log, "./log.txt")
