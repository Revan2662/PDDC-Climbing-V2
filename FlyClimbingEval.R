# Helper functions
split_concatenated_genotypes <- function(genotype_string) {
  #' Split concatenated genotypes into separate genotypes
  #' e.g., "G4C2 8 OEG4C2 29 OE" -> c("G4C2 8 OE", "G4C2 29 OE")
  
  if (is.na(genotype_string) || nchar(trimws(genotype_string)) == 0) {
    return(character(0))
  }
  
  # Check for concatenated pattern: "OE" followed by "G4C2" or number
  if (grepl("OE\\s*[G\\d]", genotype_string, ignore.case = TRUE)) {
    # Split on pattern where OE is followed by G4C2 or a digit
    parts <- strsplit(genotype_string, "(?<=OE)(?=\\s*[G\\d])", perl = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- parts[nchar(parts) > 0]
    
    if (length(parts) > 1) {
      return(parts)
    }
  }
  
  return(genotype_string)
}


standardize_genotype <- function(genotype) {
  #' Standardize genotype names that may have been shortened
  
  if (is.na(genotype) || nchar(trimws(genotype)) == 0) {
    return(NA)
  }
  
  genotype <- trimws(genotype)
  
  # Check for DNF - return NA to exclude
  if (grepl("DNF|did not finish", genotype, ignore.case = TRUE)) {
    return(NA)
  }
  
  # W1118 variants
  if (grepl("^W1{0,3}$|^W11$|^W111$|^W1118", genotype, ignore.case = TRUE)) {
    return("W1118")
  }
  
  # TBPH RNAi variants - catches TBPH, TPBH, TBHP, etc.
  if (grepl("^T[BPH][BPH][BPH]", genotype, ignore.case = TRUE)) {
    return("TBPH RNAi")
  }
  
  # G4C2 3 OE
  if (genotype == "3" || grepl("^G4C2\\s*3\\s*O?E?$", genotype, ignore.case = TRUE)) {
    return("G4C2 3 OE")
  }
  
  # G4C2 8 OE
  if (genotype == "8" || grepl("^G4C2\\s*8\\s*O?E?$", genotype, ignore.case = TRUE)) {
    return("G4C2 8 OE")
  }
  
  # G4C2 29 OE
  if (genotype == "29" || grepl("^G4C2\\s*29\\s*O?E?$", genotype, ignore.case = TRUE)) {
    return("G4C2 29 OE")
  }
  
  # G4C2 36 OE
  if (genotype == "36" || grepl("^G4C2\\s*36\\s*O?E?$", genotype, ignore.case = TRUE)) {
    return("G4C2 36 OE")
  }
  
  # G4C2 49 OE
  if (genotype == "49" || grepl("^G4C2\\s*49\\s*O?E?$", genotype, ignore.case = TRUE)) {
    return("G4C2 49 OE")
  }
  
  return(genotype)
}

find_all_trajectories <- function(root_directory, starting_point = NULL, show_progress = TRUE) {
  #' Find all trajectory sessions (with trajectories.csv, Animals_0.txt, OR just session folder)
  #' 
  #' @param root_directory The base directory to search
  #' @param starting_point Optional folder name to start extraction from (e.g., "Set 1")
  #' @param show_progress Show progress indicator (default: TRUE)
  #' @return Data frame with session information
  
  # Internal function to parse session folder and separate genotypes by sex
  parse_session_by_sex <- function(session_folder) {
    if (is.na(session_folder) || nchar(trimws(session_folder)) == 0) {
      return(list())
    }
    
    # Remove "session_" prefix
    cleaned <- sub("^session_", "", session_folder, ignore.case = TRUE)
    
    # Remove trailing "_[number]" at the end
    cleaned <- sub("_\\d+$", "", cleaned)
    
    # Split by underscore
    parts <- strsplit(cleaned, "_")[[1]]
    
    # Remove empty parts and trim whitespace
    parts <- trimws(parts)
    parts <- parts[nchar(parts) > 0]
    
    # Remove first part (experiment ID)
    if (length(parts) > 1) {
      parts <- parts[-1]
    } else {
      return(list())
    }
    
    # Handle cases where sex keywords are embedded in parts with spaces
    expanded_parts <- c()
    for (part in parts) {
      if (grepl("^(males?|females?)\\s+", part, ignore.case = TRUE)) {
        sex_match <- regmatches(part, regexpr("^(males?|females?)", part, ignore.case = TRUE))
        rest <- sub("^(males?|females?)\\s+", "", part, ignore.case = TRUE)
        expanded_parts <- c(expanded_parts, sex_match, rest)
      } else {
        expanded_parts <- c(expanded_parts, part)
      }
    }
    
    parts <- expanded_parts
    
    # Define sex keywords
    is_sex_keyword <- function(part) {
      tolower(trimws(part)) %in% c("males", "male", "females", "female")
    }
    
    # Normalize sex to standard form
    normalize_sex <- function(sex_raw) {
      sex_lower <- tolower(trimws(sex_raw))
      if (sex_lower %in% c("males", "male")) {
        return("Males")
      } else if (sex_lower %in% c("females", "female")) {
        return("Females")
      }
      return(NA)
    }
    
    # Find indices where sex keywords appear
    sex_indices <- which(sapply(parts, is_sex_keyword))
    
    if (length(sex_indices) == 0) {
      return(list())
    }
    
    # ============================================================
    # FIXED: Collect ALL genotypes for each sex, combining duplicates
    # ============================================================
    sex_genotypes_map <- list(
      "Males" = c(),
      "Females" = c()
    )
    
    for (sex_idx_pos in seq_along(sex_indices)) {
      sex_idx <- sex_indices[sex_idx_pos]
      sex_raw <- parts[sex_idx]
      current_sex <- normalize_sex(sex_raw)
      
      if (is.na(current_sex)) {
        next
      }
      
      # Determine range for genotypes following this sex keyword
      start_idx <- sex_idx + 1
      if (sex_idx_pos < length(sex_indices)) {
        end_idx <- sex_indices[sex_idx_pos + 1] - 1
      } else {
        end_idx <- length(parts)
      }
      
      # Extract genotypes
      if (start_idx <= end_idx && start_idx <= length(parts)) {
        genotype_parts <- parts[start_idx:min(end_idx, length(parts))]
        genotype_parts <- genotype_parts[nchar(genotype_parts) > 0]
        
        if (length(genotype_parts) > 0) {
          for (geno_part in genotype_parts) {
            genotype <- trimws(geno_part)
            if (nchar(genotype) > 0 && !is_sex_keyword(genotype)) {
              # Append to the appropriate sex group
              sex_genotypes_map[[current_sex]] <- c(sex_genotypes_map[[current_sex]], genotype)
            }
          }
        }
      }
    }
    
    # Build output list from combined genotypes
    sex_groups <- list()
    
    for (sex_name in names(sex_genotypes_map)) {
      genotypes <- sex_genotypes_map[[sex_name]]
      if (length(genotypes) > 0) {
        sex_groups[[length(sex_groups) + 1]] <- list(
          sex = sex_name,
          genotypes = paste(genotypes, collapse = ", ")
        )
      }
    }
    
    return(sex_groups)
  }
  
  # Main function body
  if (show_progress) {
    cat("Searching for trajectory sessions...\n")
  }
  
  # Find all session directories (look for folders starting with "session")
  all_dirs <- list.dirs(root_directory, recursive = TRUE, full.names = TRUE)
  session_dirs <- all_dirs[grepl("/session[^/]*$", all_dirs)]
  
  if (length(session_dirs) == 0) {
    warning("No session directories found")
    return(data.frame(
      file_path = character(0),
      set = character(0),
      date = character(0),
      line = character(0),
      session = character(0),
      sex = character(0),
      genotypes = character(0),
      has_trajectories = logical(0),
      has_animals_0 = logical(0),
      is_empty_session = logical(0),
      stringsAsFactors = FALSE
    ))
  }
  
  if (show_progress) {
    cat(sprintf("Found %d session directories. Processing...\n", length(session_dirs)))
  }
  
  all_records <- list()
  
  for (idx in seq_along(session_dirs)) {
    session_dir <- session_dirs[idx]
    
    if (show_progress && (idx %% 100 == 0 || idx == length(session_dirs))) {
      cat(sprintf("  Processed %d / %d sessions (%.1f%%)\n", 
                  idx, length(session_dirs), (idx/length(session_dirs))*100))
    }
    
    # Check for trajectories.csv
    traj_csv_path <- file.path(session_dir, "trajectories", "trajectories_csv", "trajectories.csv")
    has_trajectories <- file.exists(traj_csv_path)
    
    # Check for Animals_0.txt
    animals_0_path <- file.path(session_dir, "Animals_0.txt")
    has_animals_0 <- file.exists(animals_0_path)
    
    # Include ALL session folders, even if they lack both files
    is_empty_session <- !has_trajectories && !has_animals_0
    
    # Extract path components
    path_parts <- strsplit(session_dir, "/")[[1]]
    path_parts <- path_parts[nchar(path_parts) > 0]
    
    # Find starting point if specified
    start_idx <- 1
    if (!is.null(starting_point)) {
      start_indices <- which(path_parts == starting_point)
      if (length(start_indices) > 0) {
        start_idx <- start_indices[1]
      }
    }
    
    # Get relevant parts
    relevant_parts <- path_parts[start_idx:length(path_parts)]
    
    # Find session folder
    session_idx <- which(grepl("^session", relevant_parts, ignore.case = TRUE))
    
    if (length(session_idx) == 0) {
      next
    }
    
    # Extract components
    set_name <- NA
    date <- NA
    line <- NA
    session_folder <- relevant_parts[session_idx[1]]
    
    # Get path before session
    path_before_session <- relevant_parts[1:(session_idx[1] - 1)]
    
    # Extract Set, Date, Line
    if (length(path_before_session) >= 1) {
      set_indices <- grep("^Set", path_before_session, ignore.case = TRUE)
      if (length(set_indices) > 0) {
        set_name <- path_before_session[set_indices[1]]
        if (set_indices[1] + 1 <= length(path_before_session)) {
          date <- path_before_session[set_indices[1] + 1]
        }
        if (set_indices[1] + 2 <= length(path_before_session)) {
          line <- path_before_session[set_indices[1] + 2]
        }
      }
    }
    
    # Parse session folder
    session_data <- parse_session_by_sex(session_folder)
    
    # Create records for each sex group
    if (length(session_data) > 0) {
      for (sex_group in session_data) {
        all_records[[length(all_records) + 1]] <- data.frame(
          file_path = traj_csv_path,
          set = set_name,
          date = date,
          line = line,
          session = session_folder,
          sex = sex_group$sex,
          genotypes = sex_group$genotypes,
          has_trajectories = has_trajectories,
          has_animals_0 = has_animals_0,
          is_empty_session = is_empty_session,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine all records
  if (length(all_records) > 0) {
    final_df <- do.call(rbind, all_records)
    final_df <- unique(final_df)
    
    if (show_progress) {
      cat(sprintf("\nComplete! Found %d unique trajectory records from %d sessions.\n", 
                  nrow(final_df), length(session_dirs)))
      cat(sprintf("  Sessions with trajectories.csv: %d\n", 
                  sum(final_df$has_trajectories)))
      cat(sprintf("  Sessions with Animals_0.txt: %d\n", 
                  sum(final_df$has_animals_0)))
      cat(sprintf("  Empty sessions (no data files): %d\n", 
                  sum(final_df$is_empty_session)))
      
      # Debug: Show sample of parsed genotypes
      cat("\n=== SAMPLE OF PARSED SESSIONS ===\n")
      sample_rows <- head(final_df[, c("session", "sex", "genotypes")], 20)
      print(sample_rows)
    }
    
    return(final_df)
  } else {
    return(data.frame(
      file_path = character(0),
      set = character(0),
      date = character(0),
      line = character(0),
      session = character(0),
      sex = character(0),
      genotypes = character(0),
      has_trajectories = logical(0),
      has_animals_0 = logical(0),
      is_empty_session = logical(0),
      stringsAsFactors = FALSE
    ))
  }
}

load_and_assign_rois <- function(trajectories_df, roi_definitions = NULL, time_threshold = 18, y_threshold = 470) {
  #' Load trajectory CSV files and assign flies to ROIs with genotypes
  #' Only includes data up to time_threshold seconds
  #' Handles sessions with Animals_0.txt OR empty sessions (no trajectories.csv) as 0 flies
  #' Counts flies that cross the Y threshold within the time threshold
  #' Calculates means across replicates
  #' 
  #' @param trajectories_df Data frame from find_all_trajectories()
  #' @param roi_definitions List of ROI ranges (x min/max pairs)
  #' @param time_threshold Maximum time in seconds to include (default: 18)
  #' @param y_threshold Y pixel value threshold for climbing (default: 470, lower Y = higher position)
  #' @return Data frame with all trajectory data
  
  # Default ROI definitions
  if (is.null(roi_definitions)) {
    roi_definitions <- list(
      roi_1 = c(470, 725),
      roi_2 = c(780, 1030),
      roi_3 = c(1075, 1350)
    )
  }
  
  # Helper function to determine which ROI a single x value falls into
  get_roi_for_x <- function(x_val, roi_defs) {
    for (roi_name in names(roi_defs)) {
      roi_range <- roi_defs[[roi_name]]
      if (!is.na(x_val) && x_val >= roi_range[1] && x_val <= roi_range[2]) {
        return(roi_name)
      }
    }
    return(NA)
  }
  
  # Helper function to check if a fly crosses the Y threshold
  fly_crossed_threshold <- function(y_data, threshold) {
    valid_y <- y_data[!is.na(y_data) & !is.nan(y_data)]
    if (length(valid_y) == 0) {
      return(FALSE)
    }
    return(any(valid_y <= threshold))
  }
  
  # Add base_session column (without replicate number)
  trajectories_df$base_session <- sub("_\\d+$", "", trajectories_df$session)
  
  cat(sprintf("Processing %d trajectory records...\n", nrow(trajectories_df)))
  cat(sprintf("Time threshold: %.1f seconds\n", time_threshold))
  cat(sprintf("Y threshold for climbing: %d pixels (lower Y = higher position)\n", y_threshold))
  cat("Note: Y threshold crossing is evaluated only within the time threshold.\n")
  
  # ============================================================
  # PRE-PROCESS: Build a map of session -> ALL genotypes (across all sexes)
  # ============================================================
  session_all_genotypes <- list()
  
  for (row_idx in 1:nrow(trajectories_df)) {
    session_name <- trajectories_df$session[row_idx]
    genotypes_string <- trajectories_df$genotypes[row_idx]
    
    genotype_list <- c()
    if (!is.na(genotypes_string) && nchar(genotypes_string) > 0) {
      genotype_list_raw <- strsplit(genotypes_string, ",\\s*")[[1]]
      
      for (geno in genotype_list_raw) {
        split_parts <- split_concatenated_genotypes(geno)
        genotype_list <- c(genotype_list, split_parts)
      }
      
      genotype_list <- sapply(genotype_list, standardize_genotype, USE.NAMES = FALSE)
      genotype_list <- genotype_list[!is.na(genotype_list)]
    }
    
    if (length(genotype_list) > 0) {
      if (is.null(session_all_genotypes[[session_name]])) {
        session_all_genotypes[[session_name]] <- genotype_list
      } else {
        existing <- session_all_genotypes[[session_name]]
        new_genos <- genotype_list[!genotype_list %in% existing]
        session_all_genotypes[[session_name]] <- c(existing, new_genos)
      }
    }
  }
  
  # Debug: Show sessions and their combined genotypes
  cat("\n=== SESSION GENOTYPE MAPPING (showing first 20) ===\n")
  session_names <- names(session_all_genotypes)
  for (i in seq_len(min(20, length(session_names)))) {
    sn <- session_names[i]
    genos <- session_all_genotypes[[sn]]
    cat(sprintf("  %s: %s\n", sn, paste(genos, collapse = ", ")))
  }
  cat("\n")
  
  all_fly_data <- list()
  all_session_summary <- list()
  dnf_count <- 0
  time_filtered_count <- 0
  zero_fly_session_count <- 0
  
  processed_sessions <- list()
  
  # Process each row in trajectories_df
  for (row_idx in 1:nrow(trajectories_df)) {
    
    if (row_idx %% 50 == 0 || row_idx == nrow(trajectories_df)) {
      cat(sprintf("  Processed %d / %d records (%.1f%%)\n", 
                  row_idx, nrow(trajectories_df), (row_idx/nrow(trajectories_df))*100))
    }
    
    file_path <- trajectories_df$file_path[row_idx]
    current_sex <- trajectories_df$sex[row_idx]
    genotypes_string <- trajectories_df$genotypes[row_idx]
    has_trajectories <- trajectories_df$has_trajectories[row_idx]
    has_animals_0 <- trajectories_df$has_animals_0[row_idx]
    is_empty_session <- trajectories_df$is_empty_session[row_idx]
    session_name <- trajectories_df$session[row_idx]
    
    # Parse genotypes for THIS ROW (current sex only)
    genotype_list <- c()
    if (!is.na(genotypes_string) && nchar(genotypes_string) > 0) {
      genotype_list_raw <- strsplit(genotypes_string, ",\\s*")[[1]]
      
      for (geno in genotype_list_raw) {
        split_parts <- split_concatenated_genotypes(geno)
        genotype_list <- c(genotype_list, split_parts)
      }
      
      genotype_list <- sapply(genotype_list, standardize_genotype, USE.NAMES = FALSE)
      genotype_list <- genotype_list[!is.na(genotype_list)]
    }
    
    if (length(genotype_list) == 0) {
      next
    }
    
    # Use ALL genotypes from this session for ROI mapping
    all_session_genos <- session_all_genotypes[[session_name]]
    if (is.null(all_session_genos)) {
      all_session_genos <- genotype_list
    }
    
    # Build ROI to genotype mapping
    roi_genotype_map <- list()
    roi_names <- names(roi_definitions)
    
    if (length(all_session_genos) == 1) {
      single_genotype <- all_session_genos[1]
      for (roi_name in roi_names) {
        roi_genotype_map[[roi_name]] <- single_genotype
      }
      cat(sprintf("  Note: Single genotype '%s' assigned to all ROIs for session %s\n", 
                  single_genotype, session_name))
    } else {
      for (i in seq_along(all_session_genos)) {
        if (i <= length(roi_names)) {
          roi_genotype_map[[roi_names[i]]] <- all_session_genos[i]
        }
      }
    }
    
    # Initialize counts for genotypes in THIS ROW (current sex)
    for (geno in genotype_list) {
      all_session_summary[[length(all_session_summary) + 1]] <- data.frame(
        set = trajectories_df$set[row_idx],
        date = trajectories_df$date[row_idx],
        line = trajectories_df$line[row_idx],
        sex = current_sex,
        genotype = geno,
        session = session_name,
        base_session = trajectories_df$base_session[row_idx],
        file_path = file_path,
        n_flies = 0,
        n_crossed = 0,
        has_data = TRUE,
        stringsAsFactors = FALSE
      )
    }
    
    # Handle zero-fly sessions
    if (has_animals_0 && !has_trajectories) {
      zero_fly_session_count <- zero_fly_session_count + 1
      next
    }
    
    if (is_empty_session) {
      zero_fly_session_count <- zero_fly_session_count + 1
      cat(sprintf("  Note: Empty session (no trajectories.csv) counted as 0 flies: %s\n", 
                  session_name))
      next
    }
    
    if (!file.exists(file_path)) {
      zero_fly_session_count <- zero_fly_session_count + 1
      next
    }
    
    # Check if we've already processed this file
    session_sex_key <- paste(session_name, current_sex, sep = "_")
    
    if (!is.null(processed_sessions[[session_sex_key]])) {
      cached_flies <- processed_sessions[[session_sex_key]]
      
      for (fly_info in cached_flies) {
        if (fly_info$genotype %in% genotype_list) {
          geno_idx <- which(genotype_list == fly_info$genotype)[1]
          summary_idx <- length(all_session_summary) - length(genotype_list) + geno_idx
          
          if (summary_idx > 0 && summary_idx <= length(all_session_summary)) {
            all_session_summary[[summary_idx]]$n_flies <- 
              all_session_summary[[summary_idx]]$n_flies + 1
            
            if (fly_info$crossed_threshold) {
              all_session_summary[[summary_idx]]$n_crossed <- 
                all_session_summary[[summary_idx]]$n_crossed + 1
            }
          }
        }
      }
      next
    }
    
    # Read the CSV
    tryCatch({
      suppressWarnings({
        traj_data <- read.csv(file_path, stringsAsFactors = FALSE, row.names = NULL)
      })
      
      # Mark this session as having data
      for (j in seq_along(genotype_list)) {
        summary_idx <- length(all_session_summary) - length(genotype_list) + j
        if (summary_idx > 0 && summary_idx <= length(all_session_summary)) {
          all_session_summary[[summary_idx]]$has_data <- TRUE
        }
      }
      
      # Filter data by time threshold FIRST
      if ("time" %in% names(traj_data)) {
        time_filtered_indices <- which(traj_data$time <= time_threshold)
        if (length(time_filtered_indices) == 0) {
          processed_sessions[[session_sex_key]] <- list()
          next
        }
        traj_data <- traj_data[time_filtered_indices, ]
        time_filtered_count <- time_filtered_count + 1
      } else {
        warning(sprintf("No 'time' column found in file: %s", file_path))
        processed_sessions[[session_sex_key]] <- list()
        next
      }
      
      cols <- names(traj_data)
      x_cols <- grep("^x\\d+$", cols, value = TRUE)
      fly_numbers <- as.numeric(sub("^x", "", x_cols))
      
      fly_cache <- list()
      
      # Process each fly
      for (fly_num in fly_numbers) {
        x_col <- paste0("x", fly_num)
        y_col <- paste0("y", fly_num)
        
        if (!(x_col %in% cols && y_col %in% cols)) {
          next
        }
        
        x_data <- traj_data[[x_col]]
        y_data <- traj_data[[y_col]]
        
        valid_indices <- which(!is.na(x_data) & !is.nan(x_data))
        valid_x <- x_data[valid_indices]
        
        if (length(valid_x) == 0) {
          next
        }
        
        # Determine ROI for each valid time point
        roi_per_timepoint <- sapply(valid_x, function(x) get_roi_for_x(x, roi_definitions))
        roi_per_timepoint <- roi_per_timepoint[!is.na(roi_per_timepoint)]
        
        if (length(roi_per_timepoint) == 0) {
          dnf_count <- dnf_count + 1
          next
        }
        
        # Find the ROI where the fly spends the most time
        roi_counts <- table(roi_per_timepoint)
        assigned_roi <- names(roi_counts)[which.max(roi_counts)]
        assigned_genotype <- roi_genotype_map[[assigned_roi]]
        mean_x <- mean(valid_x)
        
        if (is.na(assigned_roi) || is.na(assigned_genotype)) {
          dnf_count <- dnf_count + 1
          next
        }
        
        # Check if fly crossed the Y threshold WITHIN TIME THRESHOLD
        crossed <- fly_crossed_threshold(y_data, y_threshold)
        
        # Calculate min Y reached within time threshold
        valid_y <- y_data[!is.na(y_data) & !is.nan(y_data)]
        min_y <- if (length(valid_y) > 0) min(valid_y) else NA
        
        # Add to cache
        fly_cache[[length(fly_cache) + 1]] <- list(
          fly_id = fly_num,
          roi = assigned_roi,
          genotype = assigned_genotype,
          mean_x = mean_x,
          crossed_threshold = crossed,
          min_y = min_y
        )
        
        # Only count and add data if this genotype belongs to current sex
        if (assigned_genotype %in% genotype_list) {
          geno_idx <- which(genotype_list == assigned_genotype)[1]
          summary_idx <- length(all_session_summary) - length(genotype_list) + geno_idx
          
          if (summary_idx > 0 && summary_idx <= length(all_session_summary)) {
            all_session_summary[[summary_idx]]$n_flies <- 
              all_session_summary[[summary_idx]]$n_flies + 1
            
            if (crossed) {
              all_session_summary[[summary_idx]]$n_crossed <- 
                all_session_summary[[summary_idx]]$n_crossed + 1
            }
          }
          
          # Add fly trajectory data with crossing info
          fly_data <- data.frame(
            time = traj_data$time,
            x = x_data,
            y = y_data,
            fly_id = fly_num,
            roi = assigned_roi,
            genotype = assigned_genotype,
            mean_x = mean_x,
            crossed_threshold = crossed,
            min_y = min_y,
            file_path = file_path,
            set = trajectories_df$set[row_idx],
            date = trajectories_df$date[row_idx],
            line = trajectories_df$line[row_idx],
            session = session_name,
            base_session = trajectories_df$base_session[row_idx],
            sex = current_sex,
            stringsAsFactors = FALSE
          )
          rownames(fly_data) <- NULL
          
          all_fly_data[[length(all_fly_data) + 1]] <- fly_data
        }
      }
      
      processed_sessions[[session_sex_key]] <- fly_cache
      
    }, error = function(e) {
      warning(sprintf("Error reading file %s: %s", file_path, e$message))
      processed_sessions[[session_sex_key]] <- list()
    })
  }
  
  # Combine all data
  if (length(all_fly_data) > 0) {
    final_df <- do.call(rbind, all_fly_data)
    rownames(final_df) <- NULL
  } else {
    final_df <- data.frame()
  }
  
  # Combine session summary
  if (length(all_session_summary) > 0) {
    session_summary_df <- do.call(rbind, all_session_summary)
    session_summary_df <- session_summary_df[session_summary_df$has_data, ]
    
    # Calculate percentage crossed per replicate
    session_summary_df$pct_crossed <- ifelse(
      session_summary_df$n_flies > 0,
      round(100 * session_summary_df$n_crossed / session_summary_df$n_flies, 1),
      NA  # Use NA instead of 0 for sessions with no flies
    )
  } else {
    session_summary_df <- data.frame()
  }
  
  # ============================================================
  # OUTPUT SUMMARY
  # ============================================================
  
  cat(sprintf("\nComplete! Loaded %d fly trajectories (filtered to %.1f seconds).\n", 
              sum(session_summary_df$n_flies), time_threshold))
  cat(sprintf("Files with trajectory data: %d\n", time_filtered_count))
  cat(sprintf("Sessions counted as 0 flies (Animals_0.txt or missing trajectories.csv): %d\n", 
              zero_fly_session_count))
  cat(sprintf("Sessions with 0 flies for specific genotypes: %d\n", 
              sum(session_summary_df$n_flies == 0)))
  cat(sprintf("Excluded %d flies (DNF or no ROI match).\n", dnf_count))
  
  # Y threshold crossing summary
  total_flies <- sum(session_summary_df$n_flies)
  total_crossed <- sum(session_summary_df$n_crossed)
  cat(sprintf("\n=== Y THRESHOLD CROSSING SUMMARY (Y <= %d within %.1f sec) ===\n", 
              y_threshold, time_threshold))
  cat(sprintf("Total flies: %d\n", total_flies))
  cat(sprintf("Flies that crossed threshold: %d (%.1f%%)\n", 
              total_crossed, 
              if(total_flies > 0) 100 * total_crossed / total_flies else 0))
  
  # ============================================================
  # REPLICATE-LEVEL DATA (per session)
  # ============================================================
  cat("\n=== REPLICATE-LEVEL DATA (per session) ===\n")
  session_summary_df <- session_summary_df[order(session_summary_df$set, 
                                                 session_summary_df$date,
                                                 session_summary_df$line,
                                                 session_summary_df$sex,
                                                 session_summary_df$genotype,
                                                 session_summary_df$session), ]
  
  display_summary <- session_summary_df[, !names(session_summary_df) %in% "has_data"]
  print(head(display_summary, 40))
  
  # ============================================================
  # MEAN ACROSS REPLICATES (by base_session/genotype/sex)
  # ============================================================
  cat(sprintf("\n=== MEAN ± SD ACROSS REPLICATES (Y <= %d within %.1f sec) ===\n", 
              y_threshold, time_threshold))
  
  # Only include replicates with flies for percentage calculations
  replicate_means <- aggregate(
    cbind(n_flies, n_crossed) ~ set + date + line + sex + genotype + base_session, 
    data = session_summary_df, 
    FUN = sum
  )
  
  # Calculate stats across replicates
  # For n_flies
  flies_stats <- aggregate(
    n_flies ~ set + date + line + sex + genotype + base_session, 
    data = session_summary_df, 
    FUN = function(x) c(
      mean = mean(x), 
      sd = if(length(x) > 1) sd(x) else NA,
      n = length(x)
    )
  )
  
  # For n_crossed
  crossed_stats <- aggregate(
    n_crossed ~ set + date + line + sex + genotype + base_session, 
    data = session_summary_df, 
    FUN = function(x) c(
      mean = mean(x), 
      sd = if(length(x) > 1) sd(x) else NA
    )
  )
  
  # For pct_crossed (only replicates with flies)
  pct_stats <- aggregate(
    pct_crossed ~ set + date + line + sex + genotype + base_session, 
    data = session_summary_df[session_summary_df$n_flies > 0, ], 
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE), 
      sd = if(length(x) > 1) sd(x, na.rm = TRUE) else NA,
      n = sum(!is.na(x))
    )
  )
  
  # Build the summary dataframe
  mean_sd_df <- data.frame(
    set = flies_stats$set,
    date = flies_stats$date,
    line = flies_stats$line,
    sex = flies_stats$sex,
    genotype = flies_stats$genotype,
    base_session = flies_stats$base_session,
    n_replicates = flies_stats$n_flies[, "n"],
    mean_flies = round(flies_stats$n_flies[, "mean"], 2),
    sd_flies = round(flies_stats$n_flies[, "sd"], 2),
    mean_crossed = round(crossed_stats$n_crossed[, "mean"], 2),
    sd_crossed = round(crossed_stats$n_crossed[, "sd"], 2),
    stringsAsFactors = FALSE
  )
  
  # Merge pct_crossed stats (may have fewer rows if some base_sessions have all 0 flies)
  pct_df <- data.frame(
    set = pct_stats$set,
    date = pct_stats$date,
    line = pct_stats$line,
    sex = pct_stats$sex,
    genotype = pct_stats$genotype,
    base_session = pct_stats$base_session,
    mean_pct_crossed = round(pct_stats$pct_crossed[, "mean"], 1),
    sd_pct_crossed = round(pct_stats$pct_crossed[, "sd"], 1),
    n_replicates_with_flies = pct_stats$pct_crossed[, "n"],
    stringsAsFactors = FALSE
  )
  
  mean_sd_df <- merge(mean_sd_df, pct_df, 
                      by = c("set", "date", "line", "sex", "genotype", "base_session"),
                      all.x = TRUE)
  
  mean_sd_df <- mean_sd_df[order(mean_sd_df$set, mean_sd_df$date, 
                                 mean_sd_df$line, mean_sd_df$sex, 
                                 mean_sd_df$genotype), ]
  
  print(mean_sd_df)
  
  # ============================================================
  # CROSSING SUMMARY BY GENOTYPE AND SEX (aggregated across all replicates)
  # ============================================================
  cat(sprintf("\n=== CROSSING SUMMARY BY GENOTYPE AND SEX (within %.1f sec) ===\n", time_threshold))
  cat("(Mean of replicate percentages)\n\n")
  
  # Calculate mean percentage across replicates for each genotype/sex
  crossing_by_geno_sex <- aggregate(
    cbind(n_flies, n_crossed, pct_crossed) ~ genotype + sex,
    data = session_summary_df[session_summary_df$n_flies > 0, ],
    FUN = function(x) c(
      total = sum(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      sd = if(length(x) > 1) sd(x, na.rm = TRUE) else NA,
      n = length(x)
    )
  )
  
  crossing_geno_sex_df <- data.frame(
    genotype = crossing_by_geno_sex$genotype,
    sex = crossing_by_geno_sex$sex,
    total_flies = crossing_by_geno_sex$n_flies[, "total"],
    total_crossed = crossing_by_geno_sex$n_crossed[, "total"],
    n_replicates = crossing_by_geno_sex$pct_crossed[, "n"],
    mean_pct_crossed = round(crossing_by_geno_sex$pct_crossed[, "mean"], 1),
    sd_pct_crossed = round(crossing_by_geno_sex$pct_crossed[, "sd"], 1),
    stringsAsFactors = FALSE
  )
  
  # Add overall percentage (total crossed / total flies)
  crossing_geno_sex_df$overall_pct <- round(
    100 * crossing_geno_sex_df$total_crossed / crossing_geno_sex_df$total_flies, 1
  )
  
  crossing_geno_sex_df <- crossing_geno_sex_df[order(crossing_geno_sex_df$genotype, 
                                                     crossing_geno_sex_df$sex), ]
  print(crossing_geno_sex_df)
  
  # ============================================================
  # CROSSING SUMMARY BY GENOTYPE (ALL SEXES)
  # ============================================================
  cat(sprintf("\n=== CROSSING SUMMARY BY GENOTYPE - ALL SEXES (within %.1f sec) ===\n", time_threshold))
  cat("(Mean of replicate percentages)\n\n")
  
  crossing_by_geno <- aggregate(
    cbind(n_flies, n_crossed, pct_crossed) ~ genotype,
    data = session_summary_df[session_summary_df$n_flies > 0, ],
    FUN = function(x) c(
      total = sum(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      sd = if(length(x) > 1) sd(x, na.rm = TRUE) else NA,
      n = length(x)
    )
  )
  
  crossing_geno_df <- data.frame(
    genotype = crossing_by_geno$genotype,
    total_flies = crossing_by_geno$n_flies[, "total"],
    total_crossed = crossing_by_geno$n_crossed[, "total"],
    n_replicates = crossing_by_geno$pct_crossed[, "n"],
    mean_pct_crossed = round(crossing_by_geno$pct_crossed[, "mean"], 1),
    sd_pct_crossed = round(crossing_by_geno$pct_crossed[, "sd"], 1),
    stringsAsFactors = FALSE
  )
  
  crossing_geno_df$overall_pct <- round(
    100 * crossing_geno_df$total_crossed / crossing_geno_df$total_flies, 1
  )
  
  crossing_geno_df <- crossing_geno_df[order(-crossing_geno_df$mean_pct_crossed), ]
  print(crossing_geno_df)
  
  # ============================================================
  # DETAILED REPLICATE TABLE FOR EXPORT
  # ============================================================
  cat("\n=== DETAILED REPLICATE TABLE (first 50 rows) ===\n")
  
  replicate_detail <- session_summary_df[, c("set", "date", "line", "sex", "genotype", 
                                             "base_session", "session", 
                                             "n_flies", "n_crossed", "pct_crossed")]
  replicate_detail <- replicate_detail[order(replicate_detail$set,
                                             replicate_detail$date,
                                             replicate_detail$line,
                                             replicate_detail$genotype,
                                             replicate_detail$sex,
                                             replicate_detail$session), ]
  print(head(replicate_detail, 50))
  
  # Check for unusual replicate counts
  cat("\n=== GROUPS WITH UNUSUAL REPLICATE COUNTS ===\n")
  unusual_reps <- mean_sd_df[mean_sd_df$n_replicates < 3 | mean_sd_df$n_replicates > 4, ]
  
  if (nrow(unusual_reps) > 0) {
    cat(sprintf("Found %d groups with unusual replicate counts\n", nrow(unusual_reps)))
    print(head(unusual_reps, 10))
  } else {
    cat("All groups have 3-4 replicates!\n")
  }
  
  if (nrow(final_df) > 0) {
    unique_flies <- unique(final_df[, c("file_path", "fly_id", "genotype", "sex", 
                                        "set", "date", "line", "crossed_threshold", "min_y")])
    
    cat("\n=== GENOTYPE DISTRIBUTION (Unique Flies) ===\n")
    print(table(unique_flies$genotype))
    
    cat(sprintf("\nTotal unique flies: %d\n", nrow(unique_flies)))
    cat(sprintf("Total trajectory observations (time points): %d\n", nrow(final_df)))
    cat(sprintf("Average time points per fly: %.1f\n", nrow(final_df) / nrow(unique_flies)))
    cat(sprintf("Time range in data: %.2f to %.2f seconds\n", 
                min(final_df$time, na.rm = TRUE), 
                max(final_df$time, na.rm = TRUE)))
    
    # Y position statistics
    cat(sprintf("\n=== Y POSITION STATISTICS (within %.1f sec) ===\n", time_threshold))
    cat(sprintf("Y threshold: %d pixels\n", y_threshold))
    cat(sprintf("Min Y reached (highest climb): %.1f pixels\n", 
                min(unique_flies$min_y, na.rm = TRUE)))
    cat(sprintf("Max min_Y (lowest climb): %.1f pixels\n", 
                max(unique_flies$min_y, na.rm = TRUE)))
    cat(sprintf("Mean min_Y reached: %.1f pixels\n", 
                mean(unique_flies$min_y, na.rm = TRUE)))
    
    expected_genotypes <- c("G4C2 3 OE", "G4C2 8 OE", "G4C2 29 OE", 
                            "G4C2 36 OE", "G4C2 49 OE", "W1118", "TBPH RNAi")
    problem_genotypes <- unique(unique_flies$genotype[!unique_flies$genotype %in% expected_genotypes])
    problem_genotypes <- problem_genotypes[!is.na(problem_genotypes)]
    
    if (length(problem_genotypes) > 0) {
      cat("\n=== PROBLEMATIC GENOTYPES ===\n")
      print(problem_genotypes)
    } else {
      cat("\nNo problematic genotypes found!\n")
    }
  } else {
    cat("\n=== NO TRAJECTORY DATA ===\n")
    cat("All sessions had Animals_0.txt or no valid flies\n")
  }
  
  # Store summaries
  attr(final_df, "session_summary") <- session_summary_df
  attr(final_df, "replicate_detail") <- replicate_detail
  attr(final_df, "mean_sd_summary") <- mean_sd_df
  attr(final_df, "crossing_by_genotype_sex") <- crossing_geno_sex_df
  attr(final_df, "crossing_by_genotype") <- crossing_geno_df
  attr(final_df, "dnf_excluded") <- dnf_count
  attr(final_df, "zero_fly_sessions") <- zero_fly_session_count
  attr(final_df, "time_threshold") <- time_threshold
  attr(final_df, "y_threshold") <- y_threshold
  
  return(final_df)
}


# Usage:
trajectories_df <- find_all_trajectories("/home/ubuntuomics/PDDC/PDDC Final")
# Run analysis
all_fly_data <- load_and_assign_rois(trajectories_df, time_threshold = 18, y_threshold = 470)

# Access summaries
replicate_detail <- attr(all_fly_data, "replicate_detail")
mean_summary <- attr(all_fly_data, "mean_sd_summary")
crossing_by_geno_sex <- attr(all_fly_data, "crossing_by_genotype_sex")
crossing_by_geno <- attr(all_fly_data, "crossing_by_genotype")

# Export
write.csv(replicate_detail, "replicate_crossing_data.csv", row.names = FALSE)
write.csv(mean_summary, "mean_across_replicates.csv", row.names = FALSE)
write.csv(crossing_by_geno_sex, "crossing_by_genotype_sex.csv", row.names = FALSE)
write.csv(crossing_by_geno, "crossing_by_genotype.csv", row.names = FALSE)