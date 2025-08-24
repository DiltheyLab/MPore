library(dplyr)
library(openxlsx)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

files <- list.files(output_dir, pattern = "_blast_results.txt$", full.names = TRUE)
data_list <- list()
Work_List <- list()
Work_list_check <- list()
thresholds <- c(0.001, 0.0001, 0.00001, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18, 1e-19, 1e-20, 1e-21, 1e-22, 1e-23, 1e-24, 1e-25)

for (file in files) {
  identifier <- sub("_blast_results.txt$", "", basename(file))
  data <- read.csv(file, sep = "\t", header = TRUE)

  data_filtered <- data %>%
    filter(evalue < 1e-25, pident >= 40)

  data_list <- append(data_list, identifier)
  Work_List[[identifier]] <- data_filtered
  Work_list_check[[identifier]] <- data_filtered
} 

identifier_df <- data.frame(Isolate = unlist(data_list), stringsAsFactors = FALSE)
for (i in seq_along(thresholds)) {
  col_name <- paste0("e_", sprintf("%02d", i + 2))
  evalue_List <- list()
  for (identifier in names(Work_List)) {
    df <- Work_List[[identifier]]
    filtered_df <- df %>% filter(evalue <= thresholds[i])
    unique_genes_count <- length(unique(filtered_df$Gene))
    evalue_List <- append(evalue_List, unique_genes_count)
  }
  identifier_df[[col_name]] <- unlist(evalue_List)
}

filtered_DataFrames <- list()

for (i in seq_along(thresholds)) {
  threshold <- thresholds[i]
  evalue_List <- list()
  for (identifier in names(Work_List)) {
    df <- Work_List[[identifier]]
    filtered_df <- df %>% 
      filter(evalue <= threshold) %>%
      group_by(Gene) %>%
      slice_min(evalue) %>%
      distinct() %>%
      ungroup()
    
    filtered_DataFrames[[paste0(identifier, "_e_", sprintf("%02d", i + 2))]] <- filtered_df
    unique_genes_count <- nrow(filtered_df)
    evalue_List <- append(evalue_List, unique_genes_count)
  }
  identifier_df[[paste0("e_", sprintf("%02d", i + 2))]] <- unlist(evalue_List)
}

#Get e_25 candidates
filtered_e_25 <- filtered_DataFrames[grep("_e_25$", names(filtered_DataFrames))]
unique_enzymes <- unique(unlist(lapply(filtered_e_25, function(df) df$Enzyme)))
Isolates <- sub("_e_25$", "", names(filtered_e_25))
enzyme_presence_df <- data.frame(Isolates = Isolates)
for (enzyme in unique_enzymes) {
  enzyme_presence_df[[enzyme]] <- 0
}
for (i in seq_along(Isolates)) {
  isolate <- names(filtered_e_25)[i]
  df <- filtered_e_25[[isolate]]
  present_enzymes <- unique(df$Enzyme)
  enzyme_presence_df[i, present_enzymes] <- 1
}

enzyme_presence_df_2 <- data.frame(Isolates = Isolates)
for (enzyme in unique_enzymes) {
  enzyme_presence_df_2[[enzyme]] <- NA
}
for (isolate_name in names(filtered_e_25)) {
  print(paste("Isolate name in filtered_e_25:", isolate_name))
  print(head(filtered_e_25[[isolate_name]]))
}

suppressWarnings({
  for (i in seq_along(filtered_e_25)) {
    isolate <- names(filtered_e_25)[i]
    df <- filtered_e_25[[isolate]]
    isolate_name <- sub("_e_25$", "", isolate)
    for (enzyme in unique(df$Enzyme)) {
      e_value <- df[df$Enzyme == enzyme, "evalue"]
      if (length(e_value) > 0) {
        enzyme_presence_df_2[enzyme_presence_df_2$Isolates == isolate_name, enzyme] <- e_value[1]
      }
    }
  }
})

filtered_df <- enzyme_presence_df_2
all_results <- data.frame()

for (iso in Isolates) {
  isolate_row <- filtered_df %>% filter(Isolates == iso)
  
  present_enzymes <- isolate_row %>%
    select(-Isolates) %>%
    select(where(~ !is.na(.))) %>%
    colnames()
  
  if (length(present_enzymes) == 0) next
  df <- Work_list_check[[iso]]
  
  result_df <- df %>%
    filter(Enzyme %in% present_enzymes)
  if (nrow(result_df) > 0) {
    result_df$Isolate <- iso
    all_results <- bind_rows(all_results, result_df)
    all_results <- all_results %>%
      arrange(Enzyme, Isolate, evalue) %>%
      distinct(Enzyme, Isolate, .keep_all = TRUE)
  }
}
if (nrow(all_results) > 0) {
  write.xlsx(all_results, file.path(output_dir, "All_Isolates_gene_loci.xlsx"))
  write.csv(all_results, file.path(output_dir, "All_Isolates_gene_loci.csv"), row.names = FALSE)
}

output_values_file <- file.path(output_dir, "Mtase_presence_e_25_values.xlsx")
output_file <- file.path(output_dir, "Mtase_presence_e25.xlsx")
csv_file <- file.path(output_dir, "Mtase_presence_e_25_values.csv")
write.xlsx(enzyme_presence_df_2, file = output_values_file, rowNames = FALSE)
write.xlsx(enzyme_presence_df, file = output_file, rowNames = FALSE)
write.csv(enzyme_presence_df_2, file= csv_file, row.names = FALSE)

