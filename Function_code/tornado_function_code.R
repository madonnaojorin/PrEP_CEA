############ Final code after checking ##############
############ Dec 27

theme_set(theme_classic(base_size = 20,base_family = "Helvetica"))

ggplot_tornado <- function(dat,
                           baseline_output = NA,
                           annotate_scale = 0,
                           ORDER = TRUE){
  
  if (all(class(dat) != "tornado")) stop("Input data must be tornado class data frame.")
  if (length(baseline_output) != 1) stop("Input baseline_output must be length one.")
  
  output_name <- attr(dat, "output_name")
  
  if (is.na(baseline_output)) {
    
    baseline_output <- attr(dat, "baseline")[output_name]
  }
  
  dat$baseline <- unlist(baseline_output, use.names = FALSE)
  
  # don't strictly need this
  # order output columns as decending and ascending
  datplot <-
    dat[ ,c(output_name, "baseline")] %>%
    dplyr::mutate("min" = apply(., 1, min),
                  "max" = apply(., 1, max)) %>%
    dplyr::select(min, max)
  
  datplot <- cbind(dat, datplot)
  NAMES <- as.character(datplot$names)
  
  # order by length of bars
  ##TODO## assumes symmetrical; what about otherwise?
  if (ORDER) {
    
    datplot$names = factor(as.character(datplot$names),
                           levels = rev(unique(datplot$names[order(dat$size, decreasing = TRUE)])))
  }
  
  # check if parameter values are provided
  if (all(NAMES %in% names(datplot))) {
    barLabels <- datplot[, NAMES] %>%
      OpenMx::diag2vec()
  }else{barLabels <- ""}
  
  # shift annotation left or right
  nudge <- (with(datplot, eval(parse(text = output_name)) > baseline) - 0.5) * annotate_scale
  
  ggplot2::ggplot(datplot,
                  aes(names, ymin = min, ymax = max, colour = val)) +
    geom_linerange(size = 10) +
    coord_flip() +
    ylab("Incremantal cost per quality-adjusted life-year gained (USD)") +
    xlab("")+
    #geom_hline(yintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = dat$baseline, linetype = "dashed") +
    geom_hline(yintercept = 45454.5, linetype = "solid") +
    scale_y_continuous(breaks = c(25000,50000,75000), labels= c(25000,50000,75000))+
    scale_color_manual(values = c("#00468B","#ED2200"),labels = c("Low value","High value"))+
    #theme_bw() +
    theme(axis.text = element_text(face = "bold",colour = "black"),
          axis.title = element_text(face = "bold",size=20),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.position = "right",
          legend.title = element_blank()) +
    annotate("text", x = datplot$names, y = datplot[, output_name] + nudge, label = barLabels)
}
