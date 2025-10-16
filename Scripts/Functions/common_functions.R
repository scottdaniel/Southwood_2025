
resetSampleSheet <- function(s) {
  s <- s %>%
    filter(Keep) %>%
    filter(!isControl) %>%
    mutate(SubjectID = factor(SubjectID))
  return(s)
}

se <- function(x) {
  sd(x) / sqrt(length(x))
}

filter_low_coverage <- function(props, frac_cutoff=0.6, min_ab=0, na.rm = T){
  frac_nonzero <- function (x) sum(x > min_ab, na.rm = na.rm) / length(x)
  apply(props, 1, frac_nonzero) >= frac_cutoff
}

filter_low_counts <- function(count_matrix, min_count=5) {
  keep <- count_matrix > min_count
  return(count_matrix[keep])
}

as_proportion <- function (x) x / sum(x)

logit <- function(x) {
  y = log(x / (1-x))

  return(y)

}

spread_to_numeric_matrix <- function (data, row_key, col_key, value) {
  data <- dplyr::select_(data, row_key, col_key, value)
  data_wide <- tidyr::spread(data, col_key, value, fill=0)
  data_wide <- tibble::column_to_rownames(data_wide, row_key)
  as.matrix(as.data.frame(data_wide))
}
#linear mixed-effect model functions
# tidy lmer for the lmTest object
tidy_lmer <- function(lmer_test) {
  mod <- summary(lmer_test)
  data.frame(term  = rownames(mod$tTable), mod$tTable, row.names=NULL)
}

# count based lmer (for taxa abundances) with random effects
run_lmer <- function(props_toTest, s_toTest, form1, rep_mes_label, p_cutoff) {
  rep_mes_form <- paste("~ 1 |", rep_mes_label)
  props_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / Read_Counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    mutate(props_logit_scaled = scale(props_logit)[,1]) %>%
    do(tidy_lmer(nlme::lme(as.formula(form1), random = as.formula(rep_mes_form), data=., na.action=na.omit))) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value<p_cutoff)
}

#just get the lmer
get_lmer <- function(cts_toTest, s_toTest, form1, rep_mes_label, p_cutoff) {
  rep_mes_form <- paste("~ 1 |", rep_mes_label)
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / Read_Counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    nlme::lme(as.formula(form1), random = as.formula(rep_mes_form), data=., na.action=na.omit) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value<p_cutoff)
}

#simple linear model functions
# tidy lm for the lmTest object
tidy_lm <- function(lm_test) {
  mod <- summary(lm_test)
  data.frame(term  = rownames(mod$coefficients), mod$coefficients, row.names=NULL)
}

# count based lm (for taxa abundances)
run_lm <- function(cts_toTest, s_toTest, form1, p_cutoff) {
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / Read_Counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    do(tidy_lm(lm(as.formula(form1), data=., na.action=na.omit))) %>%
    setNames(c("Taxa","term","Estimate","Std.Error","t.value","p.value")) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value < p_cutoff)
}

# just get the lm itself
get_lm <- function(cts_toTest, s_toTest, form1) {
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / Read_Counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    lm(as.formula(form1), data=., na.action=na.omit)
}

#when you already have in long format

#have to do something like this beforehand:
# summed_props_long <- summed_props %>%
#   melt() %>%
#   setNames(c("Taxa", "SampleID", "props")) %>%
#   fix_zeros_v2() %>%
#   mutate(props_logit = log(props/(1-props))) %>%
#   select(SampleID, Taxa, props, props_logit)

#when you already have in long format
run_lm_on_props <- function(factor_to_test = "geneID", props_toTest, s_toTest, form1, p_cutoff) {
  props_toTest %>%
    inner_join(s_toTest, by = "SampleID") %>% #this has the added effect of dropping all the samples that we do not want to test
    group_by(!!sym(factor_to_test)) %>%
    do(tidy_lm(lm(as.formula(form1), data=., na.action=na.omit))) %>%
    setNames(c(factor_to_test,"term","Estimate","Std.Error","t.value","p.value")) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value < p_cutoff)
}

#lmer version

run_lmer_on_props <- function(factor_to_test = "geneID", props_toTest, s_toTest, form1, rep_mes_label, p_cutoff) {

  rep_mes_form <- paste("~ 1 |", rep_mes_label)

  props_toTest %>%
    inner_join(s_toTest, by = "SampleID") %>% #this has the added effect of dropping all the samples that we do not want to test
    group_by(!!sym(factor_to_test)) %>%
    do(tidy_lmer(nlme::lme(as.formula(form1), random = as.formula(rep_mes_form), data=., na.action=na.omit))) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value<p_cutoff)
}

# run glm on dataframe
run_glm_on_props <- function(factor_to_test = "geneID", props_toTest, s_toTest, form1, p_cutoff) {
  props_toTest %>%
    inner_join(s_toTest, by = "SampleID") %>% #this has the added effect of dropping all the samples that we do not want to test
    group_by(!!sym(factor_to_test)) %>%
    do(tidy_lm(glm(as.formula(form1), data=., na.action=na.omit, family=binomial))) %>%
    setNames(c(factor_to_test,"term","Estimate","Std.Error","z.value","p.value")) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value < p_cutoff)
}

# get a glm (as kyle says: It's like lm() but with a g in front of it!
get_glm <- function(cts_toTest, s_toTest, form1) {
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / read_counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    glm(as.formula(form1), data=., na.action=na.omit, family=binomial)
}



se <- function(x) sd(x)/sqrt(length(x))

top_table <- function(summed_props, s_toPlot, thre=0.8, option=1, prop_cut=0.01){

  s_props <- summed_props[,s_toPlot$SampleID]

  if (option == 1) {
    rows_to_keep <- filter_low_coverage(s_props, frac_cutoff=thre)
  } else if (option == 2) {
    rows_to_keep <- apply(s_props,1,max) >= prop_cut
  }

  s_props <- as.data.frame(s_props[rows_to_keep,])

}

### clustered taxa heatmap
heatmap_cluster_taxa <- function(summed_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1, prop_cut=0.01, satu_limit=0.4){

  #color = saturated_rainbow(101)
  color = saturated_rainbow(101, saturation_limit=satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))

  heatmap_props <- summed_props[,heatmap_s$SampleID]

  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, frac_cutoff=thre)
  } else if (option == 2) {
    rows_to_keep <- apply(heatmap_props,1,max) >= prop_cut
  }
  heatmap_props <- heatmap_props[rows_to_keep,]

  ## group the SampleIDs
  heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]

  ## update the annotation
  anno <- heatmap_s[,grps] %>% as.data.frame()
  rownames(anno) <- heatmap_s$SampleID
  colnames(anno) <- grps

  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = anno, color = color, breaks = breaks, filename = fname,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = TRUE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = anno, color = color, breaks = breaks,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = TRUE,cellheight = 8, cellwidth = 8)
}

change_data_format <- function(d) {
  #for dates like 20181121
  #changes to 11-21-2018 (but should be sorted on the 20181121 version)
  #we have 2019-03-29 and we want 03-29-2019
  paste(substr(d,6,7), substr(d,9,10), substr(d,1,4), sep="-")
}

#
#  make_pcoa_plot <- function(uu, s, shape_by, color_by, title)
#  uu: distance, s: mapping file, shape_by: variable used for shape, color_by: variable used for color
#

make_pcoa_plot <- function(dm, s, shape_by, color_by) {
  dm <- usedist::dist_subset(dm, s$SampleID)
  pc <- pcoa(dm)
  pc_df <- merge(s, pc$vectors[, 1:3], by.x="SampleID", by.y="row.names")
  pc_pct <- round(pc$values$Relative_eig * 100)

  pcoa_plot = ggplot(pc_df, aes(x=Axis.1, y=Axis.2)) +
    theme_bw() +
    scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_colour_brewer(name=sub("_", " ", color_by), palette = "Set2") +
    labs(
      x=paste0("PCoA axis 1 (", pc_pct[1], "%)"),
      y=paste0("PCoA axis 2 (", pc_pct[2], "%)")
    )

  if (is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else {
    pcoa_plot <- pcoa_plot + geom_point()
  }
  return(pcoa_plot)
}

heatmap_grouped <- function(summed_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1, prop_cut=0.01, satu_limit=0.4){

  #color = saturated_rainbow(101)
  color = saturated_rainbow(101, saturation_limit=satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))

  heatmap_props <- summed_props[,heatmap_s$SampleID]

  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, frac_cutoff=thre)
  } else if (option == 2) {
    rows_to_keep <- apply(heatmap_props,1,max) >= prop_cut
  }
  heatmap_props <- heatmap_props[rows_to_keep,]

  ## group the SampleIDs
  grps_symbol <- syms(grps)
  heatmap_s %<>% arrange(!!!grps_symbol)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]

  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps

  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
}

fix_zeros_v2 <- function(props_long, props_label = props) {

  props_label <- enquo(props_label)

  props_long %>%
    mutate(!!props_label := !!props_label + min(filter(props_long, !!props_label > 0) %>% pull(!!props_label)) / 10)

}


### These functions to be used like this:
# vars <- c("shannon","richness","faith_pd") %>%
#   set_names(c("shannon","richness","faith_pd"))
# models <- vars %>%
#   map(diversity_lmers, s_toTest = s_toTest, variable_toTest = "Simple_Genotype")
# models <- vars %>%
#   map(diversity_lmers, s_toTest = s_toTest, variable_toTest = "original_study_group * Timepoint", random_effects = "MouseID / Cage")
# better <- bind_rows(models, .id = "Diversity measure")

diversity_lms <- function(response, s_toTest, variable_toTest) {

  form1 = paste0(response, " ~ ", variable_toTest)

  df <- tidy_lm(lm(as.formula(form1), data = s_toTest, na.action = na.omit))

  return(df)

}

diversity_lmers <- function(response, s_toTest, variable_toTest, random_effects) {

  form1 = paste0(response, " ~ ", variable_toTest)

  randform <- paste0(" ~ 1|", random_effects)

  df <- tidy_lmer(nlme::lme(as.formula(form1), random = as.formula(randform), data = s_toTest, na.action = na.omit))

  return(df)

}

tidy_permanova <- function(anov){
  data.frame(Term = rownames(anov$aov.tab), anov$aov.tab, row.names = NULL) %>%
    rename(p.value = Pr..F.)
}
permanova_test <- function(dist_matrix, s_toTest, form1, perm, strata=NULL){
  set.seed(42)
  if (!grepl("~", form1)) {
    form1 <- paste0("dist_matrix ~ ", form1)
  }
  dist_matrix <- usedist::dist_subset(dist_matrix, s_toTest$SampleID)
  form1 <- as.formula(form1)
  if(is.null(strata)) {
    tidy_permanova(adonis(form1, data=s_toTest, permutations=perm))
  } else {
    tidy_permanova(adonis(form1, data=s_toTest, permutations=perm, strata=s_toTest[,strata]))
  }
}

permanova_posthoc <- function(dist_matrix, s_toTest, form1, perm, strata=NULL, group_label, p_cutoff=0.05){
  if (!grepl("~", form1)) {
    form1 <- paste0("dist_matrix ~ ", form1)
  }
  a_ixn <- permanova_test(dist_matrix, s_toTest, form1, perm, strata) %>%
    mutate(comparison = "all") %>%
    mutate(group1 = NA) %>%
    mutate(group2 = NA)

  combs <- combn(as.character(unique(s_toTest[[group_label]])), 2)
  num_tests <- dim(combs)[2]

  # do post hoc tests
  if (filter(a_ixn, Term == group_label)$p.value < p_cutoff) {
    for (i in 1:num_tests){
      s_temp <- filter(s_toTest, .data[[group_label]] %in% combs[,i])
      #dist_toTest = dist_subset(dist_matrix, s_temp$SampleID)
      a_ixn <- rbind(a_ixn,
                     permanova_test(dist_matrix, s_temp, form1, perm, strata) %>%
                       mutate(comparison = paste(combs[,i], collapse=' - ')) %>%
                       mutate(group1 = combs[,i][1]) %>%
                       mutate(group2 = combs[,i][2])
      )
    }
  }
  a_ixn
}

# from Ceylan's script: https://github.research.chop.edu/tanesc/igram/blob/master/merge_kraken.R

merge_kraken <- function(p1, p2){

  if (is.character(p1)){
    o_r1 <- read.delim(p1, skip=1, stringsAsFactors = F)
    colnames(o_r1) <- sub(".taxa", "", colnames(o_r1))
  } else {
    o_r1 <- p1
  }

  o_r2 <- read.delim(p2, skip=1, stringsAsFactors = F)
  colnames(o_r2) <- sub(".taxa", "", colnames(o_r2))

  o <- merge(o_r1, o_r2, by=c("X.OTU.ID"), all=T) %>%
    mutate(Consensus.Lineage.x = ifelse(is.na(Consensus.Lineage.x), Consensus.Lineage.y, Consensus.Lineage.x)) %>%
    select(-Consensus.Lineage.y) %>%
    rename(Consensus.Lineage = Consensus.Lineage.x)

  o[is.na(o)] <- 0
  o
}

#use like:
# o <- merge_kraken(p1, p2)
# o <- o %>%
#   select(-Consensus.Lineage, everything()) %>%
#   rename(`#OTU.ID`=X.OTU.ID)
# fp_out <- file.path(data_dir, "20191210_kraken_merged.tsv")
# cat("# Constructed from biom file",file=fp_out,sep="\n")
# write.table(o, fp_out, sep='\t', quote=F, row.names=F, append=T)

# This is so you can output your tables to csv to share with collaborators
# e.g. table_export(my_data_frame, "report_table")
# NOTE: table_iterator is not to be used directly
table_iterator <- function() {
  table_count <- 0
  function(x, title = "table_export", ...) {
    table_count <<- table_count + 1
    table_name <- paste0(title, table_count, ".csv")
    table_path <- here("Output", table_name)
    write_csv(x = x, file = table_path, ...)
    print(paste0("This table is saved as ", table_name, "."))
  }
}
table_export <- table_iterator()

# to split up a graph that has too many facets (Taxa, genes, etc.)
chunk_it <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))