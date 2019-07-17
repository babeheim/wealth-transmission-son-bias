
rm(list = ls())

source("./project_support.r")

tic.clearlog()

dir_init("./temp")

tic("fit all models")
datasets_available <- list.files("./data", full.names = TRUE)
analyses_to_run <- datasets_available
analyses_to_run <- gsub("\\./data", "./site_scripts", analyses_to_run)
analyses_to_run <- gsub("\\.csv", ".r", analyses_to_run)
x <- sapply(analyses_to_run, source)
toc(log = TRUE)

tic("prep output table")
available_reports <- list.files("./temp", pattern = "\\.json$", full.names = TRUE)
res <- lapply(available_reports, function(z) read_json(z, simplifyVector=TRUE))
res <- fromJSON(as.character(toJSON(res)), simplifyVector=TRUE)
res <- yamltools::vectorize(res)
toc(log = TRUE)

write.csv(res, "./temp/table1.csv", row.names = FALSE)

tic.log(format = TRUE)
msg_log <- unlist(tic.log())
msg_log <- paste0("- ", msg_log)
header <- c(
  "project: wealth transmission model-fitting", 
  "events:")
msg_log <- c(header, msg_log)
writeLines(msg_log, "./log.txt")
