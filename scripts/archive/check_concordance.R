library(tidyverse)

# Load WI GT data
wi_gt <- data.table::fread(
    "SNP_alt_gts.txt",
    col.names = c("chrom", "pos", "WI_gt")
    )%>%
    #create a unique id for each row by combining chrom and pos and then remove chrom and pos
    mutate(id = paste0(chrom, "_", pos)) %>%
    select(-chrom, -pos)%>%
    #remove the everything before the "=" in the gt column
    mutate(WI_gt = gsub(".*=", "", WI_gt))


# Load the seq GT data 
seq_gt <- data.table::fread(
    "seq_WI_fil_gts.txt",
    col.names = c("chrom", "pos", "seq_gt")
    )%>%
    #create a unique id for each row by combining chrom and pos and then remove chrom and pos
    mutate(id = paste0(chrom, "_", pos)) %>%
    select(-chrom, -pos)%>%
    #remove the everything before the "=" in the gt column
    mutate(seq_gt = gsub(".*=", "", seq_gt))

# left join the two dataframes

concordance <- left_join(wi_gt, seq_gt, by = "id")


# Remove the rows where seq_gt is NA or is missing for the sample ./.
concordance_fil <- concordance %>%
    filter(!is.na(seq_gt)) %>%
    filter(seq_gt != "./.")

# create a new column that is TRUE if the WI_gt and seq_gt are the same and FALSE if they are different
concordance_fil <- concordance_fil %>%
    mutate(concordant = WI_gt == seq_gt)

# calculate the total sites and the concordance rate
concordance_rate <- concordance_fil %>%
    summarise(
        total_sites = n(),
        concordant = sum(concordant),
        discordant = total_sites - sum(concordant),
        concordance_rate = sum(concordant) / total_sites
    )

head(concordance_rate)
