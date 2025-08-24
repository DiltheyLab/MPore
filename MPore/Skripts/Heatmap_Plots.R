library(pheatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Biostrings)
library(tidyr)
library(tibble)
library(ggnewscale)
library(patchwork)
library(scales)
library(argparse)
library(gridExtra)
library(stringr)


add_marker_to_motif <- function(motif, position, strand) {
  if (strand  == "-"){
    motif <- reverseComplement(DNAString(motif))
  }
  if (position >= 1 && position <= nchar(motif)) {
    part_before <- substr(motif, 1, position - 1)
    part_at_position <- substr(motif, position, position)
    part_after <- substr(motif, position + 1, nchar(motif))
    marked_motif <- paste0(part_before, "\"", part_at_position, "\"", part_after)
    return(marked_motif)
  }
  return(motif)
}



args <- commandArgs(trailingOnly = TRUE)
Output_2 <- args[1]
Context <- args[2]
Gene_df <- args[3]






file_pattern <- "enzyme_presence_df.csv"
file_paths <- list.files(path = Output_2, 
                         pattern = file_pattern, 
                         full.names = TRUE, 
                         recursive = TRUE)


filtered_paths <- file_paths[grepl("plots_6mA/enzyme_presence_df\\.csv$", file_paths)]
gene_df <- read.csv(Gene_df)

all_dataframes <- list()
Enzyme_combined_data <- data.frame()
combined_df <- data.frame()
combined_VRE <- data.frame()

for (path in filtered_paths) {
  Dataframe <- read.csv(path)
  colnames(Dataframe) <- sub("^X", "", colnames(Dataframe))  
  Enzyme_subset <- Dataframe[, c("Enzyme",colnames(Dataframe)[2:6])]
  Enzyme_combined_data <- rbind(Enzyme_combined_data, Enzyme_subset)
  
  all_dataframes[[length(all_dataframes) + 1]] <- Dataframe
}
combined_VRE <- Reduce(function(x, y) merge(x, y, by = "Enzyme", all = TRUE), all_dataframes)

subset_list <- list()
for (df in all_dataframes) {
  
  subset <- df[, c("Enzyme", colnames(df)[7:ncol(df)])]
  subset_list[[length(subset_list) +1]] <- subset
}
combined_df <- subset_list[[1]]
if (length(subset_list) >= 2) {
  for (i in 2:length(subset_list)) {
    combined_df <- merge(
      combined_df,
      subset_list[[i]],
      by = "Enzyme",
      all = TRUE,
      suffixes = c("", paste0("_VRE", i))
    )
  }
}
combined_df <- combined_df[!duplicated(combined_df$Enzyme), ]
Methyltypes <- na.omit(combined_VRE$MethylationType)

isolates <- colnames(combined_df)[-1] 

if (all(grepl("^barcode\\d+$", isolates))) {
  isolates <- isolates[order(as.numeric(sub("barcode", "", isolates)))]
} else {
  isolates <- isolates[order(isolates)]
}


heatmap_matrix1 <- matrix(0, nrow = nrow(combined_df), ncol = length(isolates),
                          dimnames = list(combined_df$Enzyme, isolates))
colnames(heatmap_matrix1) <- gsub ("barcode", "Isolate", colnames(heatmap_matrix1))

for (isolate in isolates) {
  heatmap_matrix1[, isolate] <- as.numeric(!is.na(combined_df[[isolate]]))
}
combined_long <- melt(combined_df, id.vars = "Enzyme")
colnames(combined_long) <- c("Enzyme", "Isolate", "E_value")
heatmap_data_matrix <- melt(heatmap_matrix1)
colnames(heatmap_data_matrix) <- c("Enzyme", "Isolate", "Presence")
heatmap_data_matrix_A <- merge(heatmap_data_matrix, combined_long, by = c("Enzyme", "Isolate"), all.x = TRUE)



isolates <- colnames(combined_df)[-1] 

Path_to_Score_Data = Output_2
load_isolate <- gsub("Isolate", "barcode", isolates)
if ("5mC" %in% Methyltypes){
  file_names <- list()
  for (isolate in load_isolate) {
    pattern <- paste0("Sample_DF_", isolate, "_5mC.csv")
    files <- list.files(path = Output_2, pattern = pattern, full.names = TRUE)
    file_names[[isolate]] <- files
    
    
    
      
    
    
  }
  
  result_df <- data.frame(Key = character(0))
  
  for (isolate in names(file_names)) {
    for (file_path in file_names[[isolate]]) {
      df <- read.csv(file_path)
      key_column <- df$Key
      avg_methylation_column <- df$Avg_Methylation
      if (nrow(result_df) == 0) {
        result_df <- data.frame(Key = key_column)
      }
      result_df[[isolate]] <- sapply(result_df$Key, function(key_value) {
        matched_value <- avg_methylation_column[match(key_value, key_column)]
        if (length(matched_value) > 0) {
          return(matched_value)
        } else {
          return(NA) 
        }
      })
    }
  }
  
  result_df <- result_df[grep("C", result_df$Key), ]
  result_df$Key <- gsub("\\s*\\(\\+\\)$", "", result_df$Key)
  result_df$Key <- gsub("\\s*\\(\\-\\)$", "", result_df$Key)
  
  
  enzyme_5mC <- Enzyme_combined_data[Enzyme_combined_data$MethylationType == "5mC", ]
  enzyme_5mC$MarkedMotif <- mapply(add_marker_to_motif, enzyme_5mC$Motif, enzyme_5mC$Position, enzyme_5mC$Strand)
  
  Motifs_5mC <- c(enzyme_5mC$MarkedMotif)
  result_df <- result_df[result_df$Key %in% Motifs_5mC,]
  
  
  result_df_combined <- result_df %>%
    group_by(Key) %>%
    summarise(across(where(is.numeric), mean))
  
  filtered_dataframes <- list()
  for (df in all_dataframes) {
    filtered_dataframes[[length(filtered_dataframes) + 1]] <- 
      df[, c("Enzyme", "EnzymeUniqueMethylationActivity", "Motif", "MethylationType", "Position", "Strand")]
  }
  combined_enzymes <- do.call(rbind, filtered_dataframes)
  
  
  combined_enzymes <- unique(combined_enzymes)
  combined_enzyme_5mC <- combined_enzymes[combined_enzymes$MethylationType == "5mC", ]
  
  
  combined_enzyme_5mC$Motif <- mapply(function(motif, position) {
    if (position > 0 & position <= nchar(motif)) {
      marked_char <- substr(motif, position, position)
      new_motif <- paste0(substr(motif, 1, position - 1), '"', marked_char, '"', substr(motif, position + 1, nchar(motif)))
      return(new_motif)
    } else {
      return(motif)
    }
  }, combined_enzyme_5mC$Motif, combined_enzyme_5mC$Position)
  
  heatmap_data_5mC <- result_df_combined[, -which(colnames(result_df_combined) == "Key")] 
  heatmap_data_5mC <- as.data.frame(heatmap_data_5mC)
  rownames(heatmap_data_5mC) <- result_df_combined$Key
  
  
  
  
  heatmap_matrix_5mC <- as.matrix(heatmap_data_5mC)
  colnames(heatmap_matrix_5mC) <- gsub("barcode", "Isolate", colnames(heatmap_matrix_5mC))
  
  
  marking_matrix <- matrix(0, nrow = nrow(heatmap_matrix_5mC), ncol = ncol(heatmap_matrix_5mC))
  rownames(marking_matrix) <- rownames(heatmap_matrix_5mC)
  colnames(marking_matrix) <- colnames(heatmap_matrix_5mC)
  
  
  for (i in 1:nrow(combined_enzyme_5mC)) {
    enzyme <- combined_enzyme_5mC$Enzyme[i]
    motif  <- combined_enzyme_5mC$Motif[i]
    
    if (enzyme %in% rownames(heatmap_matrix1)) {
      active_isolates <- colnames(heatmap_matrix1)[heatmap_matrix1[enzyme,] > 0]
    }
    if (motif %in% rownames(marking_matrix)){
      marking_matrix[motif, active_isolates] <- 1
    }
  }
  
  heatmap_df <- melt(heatmap_matrix_5mC)
  colnames(heatmap_df) <- c("Motif", "Isolate", "Methylation_Score")
  
  heatmap_df$Border <- as.factor(as.vector(marking_matrix))
  heatmap_df$Isolate <- gsub("barcode", "Isolate", heatmap_df$Isolate)
  heatmap_df$Methyltype <- "5mC"
}

file_names_6mA <- list()
if ("6mA" %in% Methyltypes) {
  file_names_6mA <- list()
  for (isolate in load_isolate) {
    pattern <- paste0("Sample_DF_", isolate, "_6mA.csv")
    files <- list.files(path = Output_2, pattern = pattern, full.names = TRUE)
    file_names_6mA[[isolate]] <- files
    
    
    
      
    
    
  }
  
  result_df_6mA <- data.frame(Key = character(0))
  
  for (isolate in names(file_names)) {
    for (file_path in file_names_6mA[[isolate]]) {
      df <- read.csv(file_path)
      key_column <- df$Key
      avg_methylation_column <- df$Avg_Methylation
      if (nrow(result_df_6mA) == 0) {
        result_df_6mA <- data.frame(Key = key_column)
      }
      result_df_6mA[[isolate]] <- sapply(result_df_6mA$Key, function(key_value) {
        matched_value <- avg_methylation_column[match(key_value, key_column)]
        if (length(matched_value) > 0) {
          return(matched_value)
        } else {
          return(NA) 
        }
      })
    }
  }
  
  result_df_6mA <- result_df_6mA[grep("A", result_df_6mA$Key), ]
  result_df_6mA$Key <- gsub("\\s*\\(\\+\\)$", "", result_df_6mA$Key)
  result_df_6mA$Key <- gsub("\\s*\\(\\-\\)$", "", result_df_6mA$Key)
  
  
  enzyme_6mA <- Enzyme_combined_data[Enzyme_combined_data$MethylationType == "6mA", ]
  enzyme_6mA$MarkedMotif <- mapply(add_marker_to_motif, enzyme_6mA$Motif, enzyme_6mA$Position, enzyme_6mA$Strand)
  
  Motifs_6mA <- c(enzyme_6mA$MarkedMotif)
  result_df_6mA <- result_df_6mA[result_df_6mA$Key %in% Motifs_6mA,]
  
  
  result_df_combined_6mA <- result_df_6mA %>%
    group_by(Key) %>%
    summarise(across(where(is.numeric), mean))
  
  filtered_dataframes <- list()
  for (df in all_dataframes) {
    filtered_dataframes[[length(filtered_dataframes) + 1]] <- 
      df[, c("Enzyme", "EnzymeUniqueMethylationActivity", "Motif", "MethylationType", "Position", "Strand")]
  }
  combined_enzymes <- do.call(rbind, filtered_dataframes)
  
  
  combined_enzymes <- unique(combined_enzymes)
  combined_enzyme_6mA <- combined_enzymes[combined_enzymes$MethylationType == "6mA", ]
  
  
  combined_enzyme_6mA$Motif <- mapply(function(motif, position, strand) {
    if (strand == "-"){
      motif <- reverseComplement(DNAString(motif))
    }
    if (position > 0 & position <= nchar(motif)) {
      marked_char <- substr(motif, position, position)
      new_motif <- paste0(substr(motif, 1, position - 1), '"', marked_char, '"', substr(motif, position + 1, nchar(motif)))
      return(new_motif)
    } else {
      return(motif)
    }
  }, combined_enzyme_6mA$Motif, combined_enzyme_6mA$Position, combined_enzyme_6mA$Strand)
  
  heatmap_data_6mA <- result_df_combined_6mA[, -which(colnames(result_df_combined_6mA) == "Key")] 
  heatmap_data_6mA <- as.data.frame(heatmap_data_6mA)
  rownames(heatmap_data_6mA) <- result_df_combined_6mA$Key
  
  
  
  
  heatmap_matrix_6mA <- as.matrix(heatmap_data_6mA)
  colnames(heatmap_matrix_6mA) <- gsub ("barcode", "Isolate", colnames(heatmap_matrix_6mA))
  
  
  marking_matrix_6mA <- matrix(0, nrow = nrow(heatmap_matrix_6mA), ncol = ncol(heatmap_matrix_6mA))
  rownames(marking_matrix_6mA) <- rownames(heatmap_matrix_6mA)
  colnames(marking_matrix_6mA) <- colnames(heatmap_matrix_6mA)
  
  for (i in 1:nrow(combined_enzyme_6mA)) {
    enzyme <- combined_enzyme_6mA$Enzyme[i]
    motif  <- combined_enzyme_6mA$Motif[i]
    
    if (enzyme %in% rownames(heatmap_matrix1)) {
      active_isolates <- colnames(heatmap_matrix1)[heatmap_matrix1[enzyme,] > 0]
    }
    if (motif %in% rownames(marking_matrix_6mA)){
      marking_matrix_6mA[motif, active_isolates] <- 1
    }
  }
  heatmap_df_6mA <- melt(heatmap_matrix_6mA)
  colnames(heatmap_df_6mA) <- c("Motif", "Isolate", "Methylation_Score")
  
  heatmap_df_6mA$Border <- as.factor(as.vector(marking_matrix_6mA))
  heatmap_df_6mA$Isolate <- gsub("barcode", "Isolate", heatmap_df_6mA$Isolate)
  heatmap_df_6mA$Methyltype <- "6mA"
}

file_names_4mC <- list()
if ("4mC" %in% Methyltypes) {
  file_names_4mC <- list()
  for (isolate in load_isolate) {
    pattern <- paste0("Sample_DF_", isolate, "_4mC.csv")
    files <- list.files(path = Output_2, pattern = pattern, full.names = TRUE)
    file_names_4mC[[isolate]] <- files
    
    
    
      
    
    
  }
  
  result_df_4mC <- data.frame(Key = character(0))
  
  for (isolate in names(file_names_4mC)) {
    for (file_path in file_names_4mC[[isolate]]) {
      df <- read.csv(file_path)
      key_column <- df$Key
      avg_methylation_column <- df$Avg_Methylation
      if (nrow(result_df_4mC) == 0) {
        result_df_4mC <- data.frame(Key = key_column)
      }
      result_df_4mC[[isolate]] <- sapply(result_df_4mC$Key, function(key_value) {
        matched_value <- avg_methylation_column[match(key_value, key_column)]
        if (length(matched_value) > 0) {
          return(matched_value)
        } else {
          return(NA) 
        }
      })
    }
  }
  
  result_df_4mC <- result_df_4mC[grep("C", result_df_4mC$Key), ]
  result_df_4mC$Key <- gsub("\\s*\\(\\+\\)$", "", result_df_4mC$Key)
  result_df_4mC$Key <- gsub("\\s*\\(\\-\\)$", "", result_df_4mC$Key)
  
  
  enzyme_4mC <- Enzyme_combined_data[Enzyme_combined_data$MethylationType == "4mC", ]
  enzyme_4mC$MarkedMotif <- mapply(add_marker_to_motif, enzyme_4mC$Motif, enzyme_4mC$Position, enzyme_4mC$Strand)
  
  Motifs_4mC <- c(enzyme_4mC$MarkedMotif)
  result_df_4mC <- result_df_4mC[result_df_4mC$Key %in% Motifs_4mC,]
  
  
  result_df_combined_4mC <- result_df_4mC %>%
    group_by(Key) %>%
    summarise(across(where(is.numeric), mean))
  
  filtered_dataframes <- list()
  for (df in all_dataframes) {
    filtered_dataframes[[length(filtered_dataframes) + 1]] <- 
      df[, c("Enzyme", "EnzymeUniqueMethylationActivity", "Motif", "MethylationType", "Position", "Strand")]
  }
  combined_enzymes <- do.call(rbind, filtered_dataframes)
  
  combined_enzymes <- unique(combined_enzymes)
  combined_enzyme_4mC <- combined_enzymes[combined_enzymes$MethylationType == "4mC", ]
  
  
  combined_enzyme_4mC$Motif <- mapply(function(motif, position, strand) {
    if (strand == "-"){
      motif <- reverseComplement(DNAString(motif))
    }
    if (position > 0 & position <= nchar(motif)) {
      marked_char <- substr(motif, position, position)
      new_motif <- paste0(substr(motif, 1, position - 1), '"', marked_char, '"', substr(motif, position + 1, nchar(motif)))
      return(new_motif)
    } else {
      return(motif)
    }
  }, combined_enzyme_4mC$Motif, combined_enzyme_4mC$Position, combined_enzyme_4mC$Strand)
  
  heatmap_data_4mC <- result_df_combined_4mC[, -which(colnames(result_df_combined_4mC) == "Key")] 
  heatmap_data_4mC <- as.data.frame(heatmap_data_4mC)
  rownames(heatmap_data_4mC) <- result_df_combined_4mC$Key
  
  
  
  
  heatmap_matrix_4mC <- as.matrix(heatmap_data_4mC)
  colnames(heatmap_matrix_4mC) <- gsub ("barcode", "Isolate", colnames(heatmap_matrix_4mC))
  
  
  marking_matrix_4mC <- matrix(0, nrow = nrow(heatmap_matrix_4mC), ncol = ncol(heatmap_matrix_4mC))
  rownames(marking_matrix_4mC) <- rownames(heatmap_matrix_4mC)
  colnames(marking_matrix_4mC) <- colnames(heatmap_matrix_4mC)
  
  for (i in 1:nrow(combined_enzyme_4mC)) {
    enzyme <- combined_enzyme_4mC$Enzyme[i]
    motif  <- combined_enzyme_4mC$Motif[i]
    
    if (enzyme %in% rownames(heatmap_matrix1)) {
      active_isolates <- colnames(heatmap_matrix1)[heatmap_matrix1[enzyme,] > 0]
    }
    if (motif %in% rownames(marking_matrix_4mC)){
      marking_matrix_4mC[motif, active_isolates] <- 1
    }
  }
  heatmap_df_4mC <- melt(heatmap_matrix_4mC)
  colnames(heatmap_df_4mC) <- c("Motif", "Isolate", "Methylation_Score")
  
  heatmap_df_4mC$Border <- as.factor(as.vector(marking_matrix_4mC))
  heatmap_df_4mC$Isolate <- gsub("barcode", "Isolate", heatmap_df_4mC$Isolate)
  heatmap_df_4mC$Methyltype <- "4mC"
}

combined_heatmap_data <- data.frame()
if (exists("heatmap_df")) {
  combined_heatmap_data <- rbind(combined_heatmap_data, heatmap_df)
  print("heatmap_df found. Data processed.")
} else {
  print("heatmap_df not found.")
}
if (exists("heatmap_df_6mA")) {
  combined_heatmap_data <- rbind(combined_heatmap_data, heatmap_df_6mA)
  print("heatmap_df_6mA found and appended.")
} else {
  print("heatmap_df_6mA not found.")
}
if (exists("heatmap_df_4mC")) {
  combined_heatmap_data <- rbind(combined_heatmap_data, heatmap_df_4mC)
  print("heatmap_df_4mC found and appended.")
} else {
  print("heatmap_df_4mC not found.")
}

for (i in combined_heatmap_data$Methyltype) {
  assign(paste0("data_", i), combined_heatmap_data %>% filter(Methyltype == i))
}


p <- ggplot()

if (exists("data_5mC")) {
  p <- p + 
    geom_tile(data = data_5mC, aes(x = Isolate, y = Motif, fill = Methylation_Score), 
              color = "grey", linewidth = 0.5)+
    geom_tile(data = subset(data_5mC, Border == "1"),
              aes(x = Isolate, y = Motif), color = "black", linewidth = 1, fill = NA)+
    scale_fill_gradient(low = "white", high = "blue",
                        name = "5mC Score", limits = c(0, 100), na.value = "grey")
}
if (exists("data_6mA")) {
  p <- p + 
    new_scale_fill() +
    geom_tile(data = data_6mA, aes(x = Isolate, y = Motif, fill = Methylation_Score), 
              color = "grey", linewidth = 0.5)+
    geom_tile(data = subset(data_6mA, Border == "1"),
              aes(x = Isolate, y = Motif), color = "black", linewidth = 1, fill = NA)+
    scale_fill_gradient(low = "white", high = "red",
                        name = "6mA Score", limits = c(0,100), na.value = "grey")
}
if (exists("data_4mC")) {
  p <- p + 
    new_scale_fill() + 
    geom_tile(data = data_4mC, aes(x = Isolate, y = Motif, fill = Methylation_Score), 
              color = "grey", linewidth = 0.5)+
    geom_tile(data = subset(data_4mC, Border == "1"),
              aes(x = Isolate, y = Motif), color = "black", linewidth = 1, fill = NA)+
    scale_fill_gradient(low = "white", high = "green",
                        name = "4mC Score", limits = c(0, 100), na.value = "grey")
}


p <- suppressMessages(p + 
  labs(title = "Isolate / methylation (Cohort)", x = "Isolate", y = "Motif") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25, color = "black", family = "Arial"),
    axis.text.y = element_text(color = "black", size = 25, family = "Arial"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  ) +
  coord_fixed(ratio = 1))


ggsave(
  file.path(Output_2,paste0("heatmap_Methylation_Score.png")),
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)

path_beta_data = Output_2


beta_file_names <- list()
Number_mod = c("4mC_1", "4mC_2", "5mC_1", "5mC_2", "6mA_1", "6mA_2")
for (number in Number_mod) {
  pattern <- paste0("Beta_coef_p_values_filt_", number, ".csv")
  files <- list.files(path = Output_2, pattern = pattern, full.names = TRUE)
  beta_file_names[[number]] <- files
}
beta_data <- list()
for (key in names(beta_file_names)) {
  combined_data <- list()
  methyltype <- strsplit(key, "_")[[1]][1] 
  for (file in beta_file_names[[key]]) {
    if (file.exists(file)) {
      file_data <- read.csv(file)
      if (nrow(file_data) > 0) {  
        file_data$methyltype_source <- methyltype
        combined_data <- append(combined_data, list(file_data))
      } else {
        message(paste("Skipping empty file:", file))
      }
    } else {
      message(paste("File not found:", file))
    }
  }
  if (length(combined_data) > 0){
    beta_data[[key]] <- do.call(rbind, combined_data)
  }
}
beta_data <- beta_data[!sapply(beta_data, function(df) all(is.na(df$Isolate)))]
Big_beta_df <- do.call(rbind, beta_data)
rownames(Big_beta_df) <- NULL

suppressWarnings(filtered_df1 <- Big_beta_df %>%
  inner_join(combined_enzymes, 
             by = c("Enzyme" = "Enzyme", "methyltype_source" = "MethylationType")))

filtered_df1 <- filtered_df1 %>%
  group_by(Enzyme, Motif) %>%
  mutate(
    strand_count = n_distinct(Strand),
    methyl_count = n_distinct(methyltype_source)
  ) %>%
  ungroup() %>%
  mutate(
    Enzyme = if_else(
      strand_count > 1 & methyl_count > 1,  
      paste(Enzyme, "_", methyltype_source, sep = ""),  
      Enzyme 
    )
  ) %>%
  select(-strand_count, -methyl_count)

Motifs_for_column <- unique(combined_heatmap_data$Motif)
motif_order <- Motifs_for_column
Enzyme_for_rows <- unique(filtered_df1$Enzyme)
Dataframe_binary_enzyme_vs_motifs <- data.frame(
  Enzyme = Enzyme_for_rows,
  matrix(0, nrow = length(Enzyme_for_rows),
         ncol = length(Motifs_for_column),
         dimnames = list(NULL, Motifs_for_column)
         
  ),
  check.names = FALSE
)
expanded_df <- merge(Dataframe_binary_enzyme_vs_motifs, 
                     data.frame(Isolate = isolates), 
                     by = NULL)
expanded_df$Isolate <- gsub ("barcode", "Isolate", expanded_df$Isolate)
expanded_df <- expanded_df[, c("Enzyme", "Isolate", colnames(Dataframe_binary_enzyme_vs_motifs)[-1])]
mutated_enzyme_names_all_data <- list()
for (df in all_dataframes) {
  df <- df %>%
    group_by(Enzyme, Motif) %>%
    mutate(
      strand_count = n_distinct(Strand),
      methyl_count = n_distinct(MethylationType)
    ) %>%
    ungroup() %>%
    mutate(
      Enzyme = if_else(strand_count > 1 & methyl_count >1 , 
                       paste(Enzyme, "_", MethylationType , sep = ""), 
                       Enzyme)
    ) %>%
    select(-strand_count, -methyl_count)
  mutated_enzyme_names_all_data[[length(mutated_enzyme_names_all_data) +1]] <- df
}
subset_list_mutated <- list()
for (df in mutated_enzyme_names_all_data) {
  
  subset <- df[, c("Enzyme", colnames(df)[7:ncol(df)])]
  subset_list_mutated[[length(subset_list_mutated) +1]] <- subset
}
combined_df <- subset_list_mutated[[1]]
if (length(subset_list_mutated) >= 2) {
  for (i in 2:length(subset_list_mutated)) {
    combined_df <- merge(
      combined_df,
      subset_list[[i]],
      by = "Enzyme",
      all = TRUE,
      suffixes = c("", paste0("_VRE", i))
    )
  }
  
}
combined_df <- combined_df[!duplicated(combined_df$Enzyme), ]


heatmap_matrix <- matrix(0, nrow = nrow(combined_df), ncol = length(isolates),
                         dimnames = list(combined_df$Enzyme, isolates))
colnames(heatmap_matrix) <- gsub ("Isolate", "barcode", colnames(heatmap_matrix))

for (isolate in load_isolate) {
  heatmap_matrix[, isolate] <- as.numeric(!is.na(combined_df[[isolate]]))
}
colnames(heatmap_matrix) <- gsub ("barcode", "Isolate", colnames(heatmap_matrix))


heatmap_long_check <- as.data.frame(as.table(as.matrix(heatmap_matrix)))
colnames(heatmap_long_check) <- c("Enzyme", "Isolate", "Value")
valid_combinations <- subset(heatmap_long_check, Value == 1)
filtered_expanded_df <- merge(expanded_df, valid_combinations[, c("Enzyme", "Isolate")], 
                              by = c("Enzyme", "Isolate"))

filtered_df1$MarkedMotif <- mapply(add_marker_to_motif,filtered_df1$Motif, filtered_df1$Position, filtered_df1$Strand)

for (i in 1:nrow(filtered_df1)) {
  current_enzyme <- filtered_df1$Enzyme[i]
  current_motif  <- filtered_df1$MarkedMotif[i]
  
  if (current_motif %in% colnames(filtered_expanded_df)) {
    filtered_expanded_df[filtered_expanded_df$Enzyme == current_enzyme, current_motif] <- 1
  }
}
filtered_df1$Isolate <- gsub("barcode", "Isolate", filtered_df1$Isolate)
filtered_expanded_df <- merge(filtered_expanded_df, filtered_df1[,c("Isolate", "Enzyme", "beta_coefficient", "methyltype_source")],
                              by = c("Isolate", "Enzyme"),
                              all.x = TRUE)
filtered_expanded_df <- filtered_expanded_df[, c("Enzyme","Isolate", "beta_coefficient","methyltype_source", colnames(Dataframe_binary_enzyme_vs_motifs))]

long_df <- filtered_expanded_df %>%
  pivot_longer(cols = -c(Enzyme, Isolate, beta_coefficient, Enzyme.1, methyltype_source),    
               names_to = "Motif",                                         
               values_to = "BinaryValue")  


long_df <-  long_df %>% 
  mutate(AdjustedBinary = BinaryValue)

for (i in seq_len(nrow(long_df))) {
  motif_i <- long_df$Motif[i]
  isolate_i <- long_df$Isolate[i]
  enzyme_i <- long_df$Enzyme[i]
  
  if (long_df$BinaryValue[i] == 1) {
    for (j in seq_len(nrow(long_df))) {
      if (i != j &&
          long_df$Isolate[j] == isolate_i &&
          long_df$Enzyme[j] == enzyme_i &&
          grepl(motif_i, long_df$Motif[j], fixed = TRUE)) {
        long_df$AdjustedBinary[j] <- 1
      }
    }
  }
}
long_df <- long_df %>%
  mutate(Enzyme_clean = str_trim(Enzyme) %>% sub("_.*$", "", .))	
long_df <- long_df %>%
  left_join(gene_df %>% select(Enzyme, Isolate, Gene), by = c("Enzyme_clean"= "Enzyme", "Isolate"))

beta_df <- filtered_expanded_df %>%
  distinct(Enzyme, Isolate, beta_coefficient, methyltype_source)


beta_df$beta_coefficient[beta_df$beta_coefficient <= 0.1] <- 0
significant_pairs <- beta_df[beta_df$beta_coefficient > 2.021, c("Enzyme", "Isolate")]

long_df_filtered <- semi_join(long_df, significant_pairs, by = c("Enzyme", "Isolate"))
beta_df_filtered <- semi_join(beta_df, significant_pairs, by=c("Enzyme", "Isolate"))

p <- suppressMessages(p + 
  labs(title = "Target site recognition motifs methylation\n(all isolates)", x = "Isolate", y = "Motif") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25, color = "black", family = "Arial"),
    axis.text.y = element_text(color = "black", size = 25, family = "Arial"),
    axis.title.x = element_text(size = 25, face = "bold", color = "black", family = "Arial"),
    axis.title.y = element_text(size = 25, face = "bold", color = "black", family = "Arial"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 25, color = "black", family = "Arial"),   
    legend.title = element_text(size = 25, color = "black", family = "Arial"),
    plot.title = element_text(size = 25, face = "bold", color = "black", family = "Arial"),
    panel.grid = element_blank()
  ) +
  scale_y_discrete(limits = rev(motif_order))+
  coord_fixed(ratio = 1))




heatmap_data_matrix_new <- melt(heatmap_matrix1)
colnames(heatmap_data_matrix_new) <- c("Enzyme", "Isolate", "Presence")

suppressWarnings(heatmap_data_matrix_new <- heatmap_data_matrix_new %>%
  left_join(combined_enzymes %>% select(Enzyme, Motif, Strand, MethylationType, Position), by = "Enzyme"))


heatmap_data_matrix_new$marked_motif <- mapply(add_marker_to_motif,heatmap_data_matrix_new$Motif, heatmap_data_matrix_new$Position, heatmap_data_matrix_new$Strand)

heatmap_data_matrix_new <- heatmap_data_matrix_new %>%
  left_join(
    combined_heatmap_data %>% select(Isolate, Motif, Methylation_Score),
    by = c("Isolate", "marked_motif" = "Motif") 
  )
heatmap_data_matrix_new <- heatmap_data_matrix_new %>%
  filter(!is.na(Methylation_Score))

beta_df <- beta_df_filtered %>%
  mutate(Enzyme = sub("_.*", "", Enzyme))

heatmap_data_matrix_new <- heatmap_data_matrix_new %>%
  left_join(
    beta_df %>% select(Enzyme, Isolate, methyltype_source, beta_coefficient),
    by = c("Enzyme", "Isolate", "MethylationType" = "methyltype_source") 
  )
heatmap_data_matrix_new <- heatmap_data_matrix_new %>%
  mutate(beta_coefficient = replace_na(beta_coefficient, 0))
heatmap_data_matrix_new <- heatmap_data_matrix_new %>%
  mutate(Border = if_else(beta_coefficient >= 2, 1, 0))

joined_df <- heatmap_data_matrix_new %>%
  left_join(gene_df %>% select(Gene, Enzyme, Isolate), by = c("Enzyme", "Isolate"))
enzyme_groups <- joined_df %>%
  filter(!is.na(Gene)) %>%
  group_by(Gene, Strand, MethylationType, Motif, Position) %>%
  mutate(n_enz = n_distinct(Enzyme)) %>%
  filter(n_enz > 1) %>%
  mutate(Best_Enzyme = Enzyme[which.max(beta_coefficient)]) %>%
  ungroup() %>%
  select(Gene, Strand, Motif, MethylationType, Enzyme, Best_Enzyme) %>%
  distinct()
enzyme_mapping <- enzyme_groups %>%
  mutate(Grouped_name = paste0(Best_Enzyme, "*")) %>%
  select(Enzyme, Grouped_name) %>%
  distinct()
final_df <- joined_df %>%
  left_join(enzyme_mapping, by = "Enzyme") %>%
  mutate(Grouped_name = if_else(is.na(Grouped_name), Enzyme, Grouped_name))

combined_df$Enzyme <- sub("_.*$", "",combined_df$Enzyme)
combined_long <- combined_df %>%
  pivot_longer(
    cols = -Enzyme,
    names_to = "Isolate",
    values_to = "E_value"
  )
final_df <- final_df %>%
  left_join(combined_long, by= c("Enzyme", "Isolate"))

final_df <- final_df %>%
  mutate(E_value_clean = E_value,
         E_value_label = ifelse(is.na(E_value), "not found", NA))


panel_A_ggplot <- ggplot(final_df, aes(x = Isolate, y = Grouped_name, fill = E_value)) +
  geom_tile(data = final_df, aes(x = Isolate, y = Grouped_name, fill = E_value), 
            color = "black", linewidth = 0.5) +
 
           
 
  scale_fill_gradient(
    low = "lightblue",  
    high = "blue", 
    na.value = NA,
    labels = scales::label_scientific(digits = 3),
    name = "E-value",
    guide = guide_colorbar(order = 1)
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = filter(final_df, is.na(E_value_clean)),
  aes(fill = E_value_label),
  color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("not found" = "grey"),
    name = "",
    guide = guide_legend(order = 2, override.aes = list(color = "black")))+
  theme(legend.text = element_text(size = 14))+
  labs(title = "Candidate methylase presence\n(all isolates)", size = 20, x = "Isolate", y = "Methylase") +
  theme_minimal(base_size = 15) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25, color = "black", family = "Arial"),  
    axis.text.y = element_text(size = 25, color = "black", family = "Arial"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.title = element_text(size = 30, face = "bold", color = "black", family = "Arial"),
    legend.position = "right",
    legend.text = element_text(size = 25, color = "black", family = "Arial"),   
    legend.title = element_text(size = 25, color = "black", family = "Arial"),
    axis.title = element_text(size = 25, family = "Arial", face = "bold", color = "black")
  ) +
  coord_fixed(ratio = 1)


isolates <- unique(long_df_filtered$Isolate)
for (isolate in isolates) {
  
  
  filtered_data <- long_df_filtered %>% filter(Isolate == isolate)
  filtered_beta_df <- beta_df_filtered %>% filter(Isolate == isolate)
  filtered_score_Data <- combined_heatmap_data %>% filter(Isolate == isolate)
  
  
  filtered_data$Motif <- factor(filtered_data$Motif, levels = motif_order)
  filtered_score_Data$Motif <- factor(filtered_score_Data$Motif, levels = motif_order)
  
  real_motifs <- filtered_data %>%
    filter(AdjustedBinary == 1) %>%
    select(Enzyme, Motif) %>%
    distinct()
  filtered_data <- suppressWarnings(filtered_data %>%
    left_join(real_motifs, by = "Enzyme"))
  
  
  df_gene_grouped <- filtered_data %>%
    group_by(Isolate, Gene, methyltype_source, beta_coefficient, AdjustedBinary, Motif.x) %>%
    mutate(
      Enzyme_count = n_distinct(Enzyme),
    
      Best_Enzyme = Enzyme[which.max(beta_coefficient)][1], 
    
      Enzyme_star = if_else(
        Enzyme_count > 1,
        paste0(Best_Enzyme, "*"),
        dplyr::first(Enzyme)  
      )
    ) %>%
    ungroup()
  
  df_look <- df_gene_grouped %>%
    group_by(Isolate, beta_coefficient, methyltype_source, Motif.x, AdjustedBinary) %>%
    mutate(
      Grouped_Enzyme = str_c(unique(Enzyme_star), collapse = "_")
    ) %>%
    ungroup()
  
  df_look_beta <- df_look %>%
    select(Enzyme, Isolate, methyltype_source, beta_coefficient, Grouped_Enzyme) %>%
    distinct() %>%
    mutate(beta_coefficient = ifelse(beta_coefficient < 0.1, 0, beta_coefficient)) %>%
    distinct(Grouped_Enzyme, .keep_all = TRUE)
  
  
  heatmap_plot <- ggplot() +
    geom_tile(data = df_look, aes(x = Motif.x, y = Grouped_Enzyme, fill = factor(AdjustedBinary)), color = "grey") + 
    geom_text(data = df_look, aes(x = Motif.x, y = Grouped_Enzyme, label = AdjustedBinary), color = "black", size = 15) +
    scale_fill_manual(values = c("0" = "white", "1" = "lightblue"), labels = c("0" = "No Presence", "1" = "Presence")) +
    theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size =15, color = "black", family = "Arial"), 
          axis.title = element_blank(),
          panel.spacing = unit(2, "lines")) +
    labs(title = paste("Enzyme activity in isolate", isolate), fill = "Presence_Indicator") +
    scale_y_discrete(labels = function(x) gsub('"', '', x)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none", plot.title = element_text(size = 25, face = "bold", color = "black", family = "Arial"))
  
  
  barplot_plot <- ggplot(data = df_look_beta, aes(x = Grouped_Enzyme, y = beta_coefficient)) +
    geom_bar(stat = "identity", fill = "lightgreen", color = "black", width = 0.7) +
    scale_y_reverse() +
    coord_flip() +
    scale_x_discrete() + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25, color = "black", family = "Arial"), axis.text.y = element_text(size=25, color = "black", family = "Arial"), axis.title.y = element_blank(),
          plot.title= element_blank(), axis.title.x = element_text(size=25 , face = "bold", color = "black", family = "Arial", hjust = 0))+
    labs(y = "Enzyme activity(beta)")
  
  
  barplot_plot_score <- ggplot(data = filtered_score_Data, aes(x = Motif, y = Methylation_Score)) +
    geom_bar(stat = "identity", fill = "pink", color = "black", width = 0.7) +
    theme_minimal() +
    scale_y_reverse(limits = c(100, 0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25, color = "black", family = "Arial"),
          axis.title.y = element_text(size=25, face= "bold", color = "black", family = "Arial", hjust = 1),
          axis.text.y = element_text(size=25),
          axis.title.x = element_blank(),
          plot.title = element_blank()) +
    labs(y = "Average methylation Score")
  
  
  top_row <- (barplot_plot + heatmap_plot) + plot_layout(widths = c(1, 3))
  bottom_row <- (plot_spacer() + barplot_plot_score) + plot_layout(widths = c(1, 3.5))
  combined_plot <- top_row / bottom_row + plot_layout(heights = c(2, 1))
  
  panel_B <- p
  
  final_plot <- (panel_A_ggplot + panel_B) / combined_plot 
  
  
  
  final_plot
  

  ggsave(file.path(Output_2,paste0("multipanel_plot", "_", isolate ,".png")), plot = final_plot, height = 20, width = 30, units = "in")
}
