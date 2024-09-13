library(tidyverse)

# Function to generate column names from the data
generate_column_names <- function(df) {
  # Rename the first two columns
  colnames(df)[1:2] <- c("chrom", "pos")
  
  # Extract new column names from the third column onwards
  new_colnames <- sapply(df[1, 3:ncol(df)], function(x) strsplit(x, "=")[[1]][1])
  
  # Assign new column names to the data frame
  colnames(df)[3:ncol(df)] <- new_colnames
  
  return(df)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

wi_gt_file <- args[1]
seq_gt_file <- args[2]
samples_file <- args[3]
#samples_file <- "test_elegans_samples.txt"

out_dir <- args[4]

#check if the output directory exists, if not create it
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}

#out_dir <- "data/proc/test_elegans"

wi_gt <- data.table::fread(
    wi_gt_file
    )%>%
    #create a unique id for each row by combining chrom and pos and then remove chrom and pos
    generate_column_names()%>%
    mutate(id = paste0(chrom, "_", pos)) %>%
    select(-chrom, -pos)%>%
    #remove the everything before the "=" in the gt column
    mutate_at(vars(-id), ~gsub(".*=", "", .))


head(wi_gt)

seq_gt <- data.table::fread(
    seq_gt_file
    )%>%
    #create a unique id for each row by combining chrom and pos and then remove chrom and pos
    generate_column_names()%>%
    mutate(id = paste0(chrom, "_", pos)) %>%
    select(-chrom, -pos)%>%
    #remove the everything before the "=" in the gt column
    mutate_at(vars(-id), ~gsub(".*=", "", .))


calculate_concordance <- function(wi_df, seq_df, wi_sample, seq_sample){
    # Restrict to the columns of interest
    wi_df <- wi_df %>% 
        dplyr::select(id, wi_sample) %>%
        rename(wi_gt = !!sym(wi_sample))
    
    seq_df <- seq_df %>% 
        dplyr::select(id, seq_sample) %>%
        rename(seq_gt = !!sym(seq_sample))
    
    # Left join the two dataframes
    both <- dplyr::left_join(wi_df, seq_df, by = "id")
    
    # Remove the rows where seq_gt is NA or is missing for the sample ./.
    both_fil <- both %>%
        filter(!is.na(seq_gt)) %>%
        filter(seq_gt != "./." & wi_gt != "./.") %>%
        filter(wi_gt != "0/0")
    
    # Create a new column that is TRUE if the WI_gt and seq_gt are the same and FALSE if they are different
    both_fil_con <- both_fil %>%
        mutate(concordant = wi_gt == seq_gt)
    
    # Calculate the total sites and the concordance rate
    concordance_rate <- both_fil_con %>%
        summarise(
            total_sites = n(),
            concordant = sum(concordant),
            discordant = total_sites - sum(concordant),
            concordance_rate = sum(concordant) / total_sites
        )
    
    return(concordance_rate)
}

# calculate_concordance(wi_gt, seq_gt, wi_sample = "N2", seq_sample = "N2")
# calculate_concordance(wi_gt, seq_gt, wi_sample ="JU360", seq_sample = "JU360")
# calculate_concordance(wi_gt, seq_gt, wi_sample = "NIC1107", seq_sample = "NIC1107")


# Run the function where all sample are compared to themselves and compared to the other samples
samples <- read_lines(samples_file)

concordance_results <- tibble::tibble(
    wi_sample = character(),
    seq_sample = character(),
    total_sites = numeric(),
    concordant = numeric(),
    discordant = numeric(),
    concordance_rate = numeric()
)

for (wi_sample in samples) {
    for (seq_sample in samples) {
        concordance_rate <- calculate_concordance(wi_gt, seq_gt, wi_sample = wi_sample, seq_sample = seq_sample)
        concordance_results <- concordance_results %>%
            add_row(
                wi_sample = wi_sample,
                seq_sample = seq_sample,
                total_sites = concordance_rate$total_sites,
                concordant = concordance_rate$concordant,
                discordant = concordance_rate$discordant,
                concordance_rate = concordance_rate$concordance_rate
            )
    }
}

concordance_results

# Write the results to a file
data.table::fwrite(
    concordance_results,
    glue::glue("{out_dir}/concordance_results.tsv"),
    sep = "\t"
)

# Plot heatmap of concordance rates
concordance_matrix <- concordance_results %>%
    dplyr::select(wi_sample, seq_sample, concordance_rate) %>%
    tidyr::spread(seq_sample, concordance_rate) %>%
    column_to_rownames("wi_sample") %>%
    as.matrix()


# Plot the heatmap
hm <- pheatmap::pheatmap(
    concordance_matrix,
    #color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    main = "Concordance Rates"
)

# Save the heatmap
ggsave(
    glue::glue("{out_dir}/concordance_heatmap.png"),
    plot = hm
)