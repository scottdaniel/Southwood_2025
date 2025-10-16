#for sunbeam sbx_gene_cluster extension results


read_gene_aln_results <- function(base_dir, pattern = "_1.txt") {

  tibble(FileName = list.files(
    base_dir, pattern=paste0("*", pattern))) %>%
    group_by(FileName) %>%
    do(read_diamond_file(file.path(base_dir, .$FileName), stringsAsFactors = F)) %>%
    ungroup() %>%
    select(-FileName) %>%

    #right_join(select(s_seq, SampleID), by="SampleID") %>%
    filter(!is.na(geneID)) %>%

    pivot_wider(id_cols = SampleID,
                names_from = c(geneID,taxon),
                names_sep = "__",
                values_from = count,
                values_fill = list(count = 0)) %>%
    pivot_longer(-SampleID, names_to = c("geneID","taxon"), names_sep = "__", values_to = "count") %>%

    mutate(database = basename(base_dir))

}

read_gene_results_no_taxon <- function(base_dir, pattern = "_1.txt") {

  tibble(FileName = list.files(
    base_dir, pattern=paste0("*", pattern))) %>%
    group_by(FileName) %>%
    do(read_diamond_file(path_to_file= file.path(base_dir, .$FileName), stringsAsFactors = F, pattern = pattern)) %>%
    ungroup() %>%
    select(-FileName) %>%

    group_by(geneID, SampleID) %>%
    summarize(sum_by_kegg_term = sum(count)) %>%
    #right_join(select(s_seq, SampleID), by="SampleID") %>%
    filter(!is.na(geneID)) %>%

    #much better than using complete(SampleID, nesting(geneID), fill = list(count=0)) %>%
    pivot_wider(id_cols = SampleID,
                names_from = geneID,
                values_from = sum_by_kegg_term,
                values_fill = list(sum_by_kegg_term = 0)) %>%
    pivot_longer(-SampleID, names_to = "geneID", values_to = "count") %>%

    mutate(database = basename(base_dir))

}

#read a single file of diamond results
read_diamond_file <- function(path_to_file, pattern = "_1.txt", ...) {

  if (!str_detect(basename(path_to_file), pattern)) {
    warning(paste("The pattern", pattern, "did not match the file name!",
                  "Your SampleID's will be the same as the filenames"))
  }

  my_df <- read.delim(path_to_file, colClasses = c("geneID" = "character", "taxon" = "character", "count" = "numeric"), ...) %>%
    mutate(SampleID = str_replace(basename(path_to_file), pattern, ""))

  my_df

}

read_table_if_no_rds <- function(name_of_table, rds_file_name, base_dir, pattern = "_1.txt") {

  if (!file.exists(here("Data",rds_file_name))) {
    name_of_table <- read_gene_results_no_taxon(base_dir, pattern)
    write_rds(name_of_table, here("Data",rds_file_name))
    return(name_of_table)
  }else{
    read_rds(here("Data",rds_file_name))
  }

}
