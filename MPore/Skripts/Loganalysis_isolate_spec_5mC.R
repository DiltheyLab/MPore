library(jsonlite)
library(dplyr)
library(tidyr)
library(stringr)
library(Biostrings)
library(tidyverse)
library(glmnet)
library(openxlsx)
library(Matrix)
load_fasta <- function(fasta_path) {
  fasta_seqs <- readDNAStringSet(fasta_path)
  seqs_list  <- as.list(fasta_seqs)
  return(seqs_list)
}

reverse_complement <- function(seq) {
  return(reverseComplement(DNAString(seq)))
}


args <- commandArgs(trailingOnly = TRUE)
#input_dir <- "/Users/azlannisar/Desktop/mount/"
input_dir <- args[1]
output_dir <- args[2]
#enzyme_df_path <- "/home/azlan/VRE_Data/TSV_Enzyme.csv"
enzyme_df_path <- args[3]
dataframe_path <- args[4]
#json_files <- list.files(input_dir, pattern = "updated_25_6mA_split_1.json", full.names = TRUE)
#updated_dict <- fromJSON(json_files, simplifyDataFrame = FALSE)
MTase_presence <- args[5]



# Define the main function
main <- function(input_dir, output_dir, enzyme_df_path, dataframe_path, MTase_presence) {
  
  # Load the enzyme data once
  tsv_enzyme_info <- read.csv(enzyme_df_path, sep = "\t")
  Data_ref <- read.csv(dataframe_path, stringsAsFactors = FALSE)
  MTase_df <- read.csv(MTase_presence)
  MTase_df <- MTase_df[ ,c(TRUE, grepl("^M[0-9]*\\..+", colnames(MTase_df)[-1]))]  
  fasta_dict <- list()
  for (i in seq_len(nrow(Data_ref))) {
    file_name <- Data_ref$File_name[i]
    fasta_path <- Data_ref$Reference_path[i]
    
    if (file.exists(fasta_path)) {
      fasta_dict[[file_name]] <- load_fasta(fasta_path)
    } else {
      warning(paste0("Datei nicht gefunden:", fasta_path))
    }
  }
  #trimm whitespaces 
  fasta_dict <- lapply(fasta_dict, function(contigs) {
    names(contigs) <- trimws(names(contigs))
    contigs <- DNAStringSet(contigs)
    names(contigs) <- sub(" .*", "", names(contigs))
    rev_contigs <- DNAStringSet(lapply(contigs, reverse_complement))
    names(rev_contigs) <- paste0(names(contigs), "_rev")
    return(c(contigs, rev_contigs))
  })
  print(fasta_dict)
  
  # Find JSON files with the pattern "updated_25_6mA_split_{number}.json"
  json_files <- list.files(input_dir, pattern = "updated_25_5mC_split_[0-9]+\\.json", full.names = TRUE)
  
  
  if (length(json_files) == 0) {
    json_files <- file.path(input_dir, "updated_25_5mC.json")
  }
  
  for (json_file in json_files) {
    # Extract the file number for output naming
    if (basename(json_file) == "updated_25_5mC.json") {
      file_number <- "all"
    } else {
      # Extract the file number from the split file name
      file_number <- str_extract(basename(json_file), "_split_\\d+") %>%
        str_extract("\\d+")
    }
    updated_dict <- fromJSON(json_file, simplifyDataFrame = TRUE)
    
    updated_dict <- lapply(updated_dict, function(df) {
      keep_first <- 1:9
      match_cols <- which(grepl("^M[0-9]*\\..+", colnames(df)))
      print(match_cols)
      keep_cols <- unique(c(keep_first, match_cols))
      df[ , keep_cols, drop = FALSE]
    })    


    # Process 6mA Data by filtering motifs containing "A"
    for (key in names(updated_dict)) {
      updated_dict[[key]] <- updated_dict[[key]] %>%
        filter(grepl("C", motif)) %>%  
        mutate(motif = gsub('"', '', motif))  
    }
    
    for (key in names(updated_dict)) {
      df <- updated_dict[[key]] %>%
        dplyr::filter(coverage >=10) %>%
        dplyr::distinct(Base, pos, strand, contig, .keep_all = TRUE)
      updated_dict[[key]] <- df
    }
    
    lookup_list <- list() 
    for (key in names(updated_dict)) {
      df_check <- updated_dict[[key]]
      fasta_seqs <- fasta_dict[[key]]
      
      fasta_flat <- lapply(fasta_seqs, as.character)
      fasta_lengths <- sapply(fasta_flat, nchar)
      
      
      unique_loci <- df_check %>%
        dplyr::select(pos, strand, contig) %>%
        dplyr::distinct()
      
      get_context_vec <- function(pos, strand, contig) {
        mapply(function(p, s, c) {
          if (s == "+") {
            seq <- fasta_flat[[c]]
            paste0(
              substr(seq, p - 1, p - 1),
              substr(seq, p + 1, p + 1)
            )
          } else {
            rev_c <- paste0(c, "_rev")
            rev_seq <- fasta_flat[[rev_c]]
            rc_pos <- fasta_lengths[[rev_c]] - p +1
            paste0(
              substr(rev_seq, rc_pos - 1, rc_pos - 1),
              substr(rev_seq, rc_pos + 1, rc_pos + 1)
            )
          }
        }, pos, strand, contig, USE.NAMES = FALSE)
      }
      unique_loci$context <- get_context_vec(unique_loci$pos + 1, unique_loci$strand, as.character(unique_loci$contig))
      lookup_list[[key]] <- unique_loci
      print(paste("Done for", key))
    }
    
    
    base_pairs <- c("AA", "AC", "AG", "AT", 
                    "CA", "CC", "CG", "CT", 
                    "GA", "GC", "GG", "GT", 
                    "TA", "TC", "TG", "TT")
    for (key in names(updated_dict)) {
      df <- updated_dict[[key]]
      
      for (pair in base_pairs) {
        df[paste0("Context_", pair)] <- 0
      }
      updated_dict[[key]] <- df
    }
    
    
    #Filling time 
    for (key in names(updated_dict)) {
      df_check <- updated_dict[[key]]
      lookup_check <- lookup_list[[key]]
      
      context_vector <- lookup_check$context
      
      
      index <- 1:nrow(lookup_check)
      context_columns <- paste0("Context_", context_vector)
      context_columns_in_df <- grep("^Context_", colnames(df_check), value = TRUE)
      
      unique_context <- paste0("Context_", unique(lookup_check$context)) 
      columns_to_remove <- base::setdiff(unique_context, context_columns_in_df)
      
      df_check_matrix <- as.matrix(df_check)
      df_check_matrix[cbind(index, match(context_columns, colnames(df_check)))] <- 1
      df_check <- as.data.frame(df_check_matrix)
      print("DONe with Filling")
      
      
      #df_check <- df_check[, !(colnames(df_check) %in% columns_to_remove)]
      
      updated_dict[[key]] <- df_check
      print(paste("Done for", key))
    }
    
    
    tsv_enzyme_info$Signature <- paste(
      tsv_enzyme_info$Motif,
      tsv_enzyme_info$MethylationType,
      tsv_enzyme_info$Strand,
      tsv_enzyme_info$Position,
      sep = "_"
    )
  
    enzyme_signature_map <- setNames(tsv_enzyme_info$Signature, tsv_enzyme_info$Enzyme) 
    #Sum up colinear entries based on signature entries from tsv rebase 
    for (i in seq_along(updated_dict)) {
      df <- updated_dict[[i]]
    
      enzyme_columns <- names(df)[
        !(names(df) %in% c("Base", "pos", "strand", "score", "motif", "coverage", "modified", "contig", "sample")) &
          !grepl("^Context:", names(df))
      ]
      enzyme_columns <- enzyme_columns[enzyme_columns %in% names(enzyme_signature_map)]
      enzyme_signatures <- enzyme_signature_map[enzyme_columns]
      signature_groups <- split(names(enzyme_signatures), enzyme_signatures)
    
      #Sum up colinear groups based on signature 
      for (group in signature_groups) {
        if (length(group) > 1) {
          combined_col_name <- paste(group, collapse = "_")
          df[[combined_col_name]] <- do.call(pmax, df[group])
          df[group] <- NULL
        }
      }
    
      updated_dict[[i]] <- df
      print(paste0("Done with ", i))
    }

    for (key in names(updated_dict)) {
      df <- updated_dict[[key]] 
      df$mod_status <- ifelse(as.numeric(df$modified) / as.numeric(df$coverage) >= 0.8, 1, 0) 
      df$check_status <- ifelse(round((as.numeric(df$modified) / as.numeric(df$coverage)) * 100, 2) == as.numeric(df$score), 1, 0) 
      df <- df[as.numeric(df$check_status) == 1,]
      updated_dict[[key]] <- df  
    }
    
    # Define logistic regression function
    calculate_pi <- function(df, beta) {
      enzyme_columns <- intersect(names(df), names(beta))
      logit_pi <- apply(df[enzyme_columns], 1, function(x) {
        sum(x * beta[enzyme_columns], na.rm = TRUE) + beta["(Intercept)"]
      })
      df$logit_pi <- logit_pi
      df$pi <- exp(logit_pi)/(1 + exp(logit_pi))
      
      return(df)
    }

    print("Done with Collinear check")

    updated_dict <- lapply(updated_dict, function(df) {
      df[ , colSums(df == 0) != nrow(df), drop = FALSE]
    })

    get_enzyme_cols <- function(df) {
      non_enzyme_cols <- c("Base", "pos", "strand", "score", "motif", "contig", "coverage", "modified",
                           "sample", "mod_status", "check_status")
      non_enzyme_cols <- c(non_enzyme_cols, grep("^Context_", colnames(df), value = TRUE))
      
      # Enzyme columns = all others
      base::setdiff(colnames(df), non_enzyme_cols)
    }
    combine_updated_dicts <- function(updated_dict) {
      
      # First: get all enzyme column names per isolate
      enzyme_lists <- lapply(updated_dict, get_enzyme_cols)
      all_enzymes <- unique(unlist(enzyme_lists))
      
      # Prepare renamed dataframes
      renamed_dfs <- list()
      
      for (isolate in names(updated_dict)) {
        df <- updated_dict[[isolate]]
        enzymes <- enzyme_lists[[isolate]]
        
        # Rename all enzyme columns (shared or unique) with isolate suffix
        new_names <- colnames(df)
        for (col in enzymes) {
          new_names[new_names == col] <- paste0(col, ".", isolate)
        }
        colnames(df) <- new_names
        
        # Add Origin column
        df$Origin <- isolate
        renamed_dfs[[isolate]] <- df
      }
      
      # Create full column union
      all_cols <- unique(unlist(lapply(renamed_dfs, colnames)))
      
      # Align all dataframes to same structure
      for (isolate in names(renamed_dfs)) {
        df <- renamed_dfs[[isolate]]
        for (col in base::setdiff(all_cols, colnames(df))) {
          df[[col]] <- "0"
        }
        df <- df[, all_cols]
        renamed_dfs[[isolate]] <- df
      }
      
      # Combine all into one
      combined_df <- dplyr::bind_rows(renamed_dfs)
      return(combined_df)
    }
    combined_df <- combine_updated_dicts(updated_dict)
    print("Done with creating singular file")
    
    #One Hot Encoding 
    origin_matrix <- model.matrix(~ Origin -1, data = combined_df)
    colnames(origin_matrix) <- sub("^Origin", "", colnames(origin_matrix))
    #combined_df <- cbind(combined_df[ , !colnames(combined_df) %in% "Origin", drop = FALSE], origin_matrix)
    combined_df <- cbind(combined_df, origin_matrix)

    print("starting the glmnet process")
    beta_coefficients2 <- list()
    p_values2 <- list()
    Std_Error2 <- list()
  
    for (key in names(updated_dict)) {
      df <- updated_dict[[key]]
      filter_Mtase <- MTase_df %>%
        dplyr::filter(Isolates == key) %>%
        dplyr::select(-Isolates) %>%
        dplyr::select(where(~ any(!is.na(.)))) %>%
        names()
    
      enzyme_columns <- names(df)[!(names(df) %in% c("Base", "pos", "strand", "score", "motif", "coverage", "modified", "pi", "contig", "sample", "mod_status", "check_status"))]
      cat("enzyme_columns for", key, ":\n")
      print(enzyme_columns)
      context_columns <- enzyme_columns[grepl("^Context_", enzyme_columns)]
      filtered_enzyme_columns <- enzyme_columns[sapply(enzyme_columns, function(col) any(sapply(filter_Mtase, function(e) grepl(e, col))))]
      enzyme_columns <- c(filtered_enzyme_columns, context_columns)

      cat("enzyme_columns for", key, ":\n")
      print(enzyme_columns)
    
      X <- df[, enzyme_columns, drop = FALSE]
      X <- data.frame(lapply(X, function(col) as.numeric(as.character(col))))
      colnames(X) <- gsub("^X(?=\\d)", "", colnames(X), perl = TRUE)
      ncol(X) 
      any(is.na(as.matrix(X)))
      y_matrix <- cbind(
        unmodified = (as.numeric(df$coverage) - as.numeric(df$modified)),
        modified = as.numeric(df$modified)
      )
      X_matrix <- Matrix(as.matrix(X), sparse = TRUE)
    
      cv_fit <- cv.glmnet(X_matrix, y_matrix, alpha = 1, nfolds = 10, family= "binomial")
      coefs <- coef(cv_fit, s = "lambda.min")
      beta_vec <- as.vector(coefs)
      names(beta_vec) <- rownames(coefs)
      Intercept <- beta_vec["(Intercept)"]
      #beta_vec <- beta_vec[beta_vec > 0]
      beta_vec <- c("(Intercept)" = Intercept, beta_vec)
    
      if (length(beta_vec) == 0) {
        beta_coefficients2[[key]] <- c("(Intercept)" = NA)
        print(paste0("Skipping ", key, " due to no valid predictors"))
        next
      }
      if (length(beta_vec) == 0) {
        beta_coefficients2[[key]] <- c("(Intercept)" = NA)
      } else {
        beta_coefficients2[[key]] <- beta_vec
      }
    
      print(paste0("Done for ", key))
    }
   #beta_vec <- beta_vec[beta_vec > 0.5]
    rows2 <- list()
    for (isolate in names(beta_coefficients2)) {
      coefficients <- beta_coefficients2[[isolate]]
      for (enzyme in names(coefficients)) {
        if (enzyme != "(Intercept)") {
          row <- data.frame(
            Isolate = isolate,
            Enzyme = enzyme,
            beta_coefficient = coefficients[enzyme],
            stringsAsFactors = FALSE
          )
          rows2[[length(rows2) + 1]] <- row
        }
      }
    }
    result_df <- do.call(rbind, rows2)
    expanded_rows <- list()
  
    for (i in 1:nrow(result_df)) {
    
      current_row <- result_df[i, ]
    
      if (grepl("_", current_row$Enzyme) && !grepl("^Context",current_row$Enzyme)) {
        enzymes <- strsplit(current_row$Enzyme, "_")[[1]]
        for (enzyme in enzymes) {
        
          new_row <- current_row
          new_row$Enzyme <- enzyme
          expanded_rows[[length(expanded_rows) + 1]] <- new_row
        }
      } else {
        expanded_rows[[length(expanded_rows) + 1]] <- current_row
      }
    }
    expanded_result_df_filt <- do.call(rbind, expanded_rows)
    row.names(expanded_result_df_filt) <- NULL    
    #write.xlsx(Excel_df, file.path(output_dir, paste0("Context_influence_5mC_", file_number, ".xlsx")))
    #print("Context 5mC influence DOne")
    write.csv(expanded_result_df_filt, file.path(output_dir, paste0("Beta_coef_p_values_filt_5mC_2", file_number, ".csv")), row.names = FALSE)
    write.xlsx(expanded_result_df_filt, file.path(output_dir, paste0("Beta_coef_p_values_filt_5mC_2", file_number, ".xlsx")))
    print("Enzymatic_influence 5mC DOne")
  }
}
main(input_dir, output_dir, enzyme_df_path, dataframe_path, MTase_presence)
