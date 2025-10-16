make_pcoa_plot_with_closed_and_open_circles_for_points <- function(dm, s, shape_by, color_by) {
  dm <- usedist::dist_subset(dm, s$SampleID)
  pc <- pcoa(dm)
  pc_df <- merge(s, pc$vectors[, 1:3], by.x="SampleID", by.y="row.names")
  pc_pct <- round(pc$values$Relative_eig * 100)

  pcoa_plot = ggplot(pc_df, aes(x=Axis.1, y=Axis.2)) +
    theme_bw() +
    scale_shape_manual(name=sub("_", " ", shape_by), values=c(1,16)) +
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

make_pcoa_plot_with_ellipses_around_groups <- function(dm, s, shape_by, color_by, label_points = "", fname = NULL, ellipse_by = NULL) {
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

  if (label_points != "") {
    pcoa_plot = pcoa_plot + geom_text(aes(label=get(label_points)), hjust=0, vjust=0) }

  if (is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else {
    pcoa_plot <- pcoa_plot + geom_point()
  }

  if (!is.null(fname)) {
    ggsave(pcoa_plot,filename = fname, device = "pdf")
  }

  if (!is.null(ellipse_by)) {
    pcoa_plot <- pcoa_plot + stat_ellipse(aes(colour=factor(get(ellipse_by))))
  }

  return(pcoa_plot)
}


make_pcoa_plot_with_ellipses_around_groups <- function(dm, s, shape_by, color_by, label_points = "", fname = NULL, ellipse_by = NULL) {
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

  if (label_points != "") {
    pcoa_plot = pcoa_plot + geom_text(aes(label=get(label_points)), hjust=0, vjust=0) }

  if (is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else {
    pcoa_plot <- pcoa_plot + geom_point()
  }

  if (!is.null(fname)) {
    ggsave(pcoa_plot,filename = fname, device = "pdf")
  }

  if (!is.null(ellipse_by)) {
    pcoa_plot <- pcoa_plot + stat_ellipse(aes(colour=factor(get(ellipse_by))))
  }

  return(pcoa_plot)
}
