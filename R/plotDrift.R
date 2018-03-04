plotPulseChase <- function(x, time_x)
{
  time_sim        = time_sim + tau
  fractions.split = nrow(sim)
  col_names       = c("Age", paste(1:fractions.split, ".", fractions.split, sep = ""))
  
  ## Sim process  
  model.preds           = data.frame(Age = time_sim, t(sim))
  colnames(model.preds) = col_names
  model.melt            = melt(model.preds, id = "Age")

  if(!is.null(x)){
    ## Data process
    clone.quarters           = data.frame(Age = time_x, t(normalise.counts(x)))
    colnames(clone.quarters) = col_names
    quarters.melt            = melt(clone.quarters, id = "Age")
    # SE calculations
    se.clone.quarters           = data.frame(Age = time_x, t(se.counts.multinom(x)))  
    colnames(se.clone.quarters) = col_names
    se.melt                     = melt(se.clone.quarters, id = "Age")
    quarters.melt               = cbind(quarters.melt, se = se.melt[,4])
    plot_data = geom_point(data = quarters.melt, aes(x=Age, y=value), size = 2.2) + 
    geom_errorbar(data = quarters.melt, aes(ymax = value + se, ymin=value - se), width=3) 
  }else{plot_data=NULL}
#   browser()  
    pp <- ggplot(data = model.melt, aes(x=Age, y=value)) + geom_line() + plot_data +
    labs(x = "Time post labelling (days)", y = "Clone fraction") +
    theme_bw() + facet_wrap(~variable, nrow = 1) + ggtitle("") # facet_grid(Location~variable) 
  pp
}



