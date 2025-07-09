dotplot_gsea <- function(fgsea_results, title = "GSEA results", top_n = 5, bottom_n=5, comparisons = NULL, category = NULL, FDR_threshold = 0.05) {
  
  if(is.null(comparisons)) comparisons = unique(fgsea_results$comparison)
  if(is.null(category)) category = unique(fgsea_results$GeneSetCollection)
  
  data = as.data.table(subset(fgsea_results, GeneSetCollection == category & FDR_global < FDR_threshold & comparison %in% comparisons))
  data$comparison = factor(data$comparison, levels=comparisons)

  data$direction = ifelse(data$NES>0, "Up","Down")
  data[, rank_bottom:=frank(NES), by=c("comparison")]
  data[, rank_top:=frank(-NES), by=c("comparison")]
  data = subset(data, ((rank_top<=top_n & direction=="Up") | (rank_bottom <=bottom_n & direction=="Down") ) )

  if(nrow(data)<1) return()

  pw_counts = as.data.table(data)[, .(pw_count=.N, comps = paste(sort(comparison),collapse=",")), by=c("pathway")]
  # shorten pw names
  pw_counts$pathway_shorter = ifelse(nchar(pw_counts$pathway)<60, pw_counts$pathway, paste0(substr(pw_counts$pathway, 1, 60), "..."))
  dup_short_pw = duplicated(pw_counts$pathway_shorter)
  pw_counts$pathway_shorter[dup_short_pw] <- paste0(pw_counts$pathway_shorter, 1:sum(dup_short_pw))
  setkey(pw_counts, pathway)
  
  data$pathway_orig = data$pathway
  data$pathway = pw_counts[data$pathway_orig, "pathway_shorter"]
  setorder(pw_counts, pw_count, comps)
  data$pathway = factor(data$pathway, levels = pw_counts$pathway_shorter)
  data$leading_edge_len = unlist(lapply(data$leadingEdge, length))
  data$leading_edge_proportion = data$leading_edge_len / data$size
  data$circle_size = data$leading_edge_proportion * 10
  
  y_label_fun = function(x) levels(data$pathway)[x]
  
  p = ggplot(data, aes(x=comparison, y= pathway)) + 
    geom_point(aes(size=circle_size, color=NES)) + 
    scale_colour_gradient2(low="blue", high="red", breaks=seq(-3, 3, by=1), limits=c(-3,3), oob=scales::squish) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=20 ), axis.text.y = element_text(size=max(10,floor(20-nrow(pw_counts)/10) ) ), 
          plot.title = element_text(size=20)) + 
    labs(title=title) + xlab("Comparison") + ylab("Gene set") + 
    scale_size_continuous(
      limits=c(0,10),
      breaks = seq(0,10,2.5),
      labels = paste(c(0,25,50,75,100), "%")) + 
    guides(size = guide_legend(title = "% in Leading Edge")) + 
    scale_x_discrete(drop = FALSE) # maintain comparisons without any significant gene sets
  
  
  return(p)
}

