#' Calculate proportions with 95\% credible interval (confidence interval).
#' If there are replicates they will be pooled(!). It will work for both count table format and list format.
#' For it to work for table format give time vector
format_exp_data = function(data_x, time_points)
{
  colnames(data_x) = paste0("T", time_points)
  data_x= data_x %>% mutate(CloneSize = paste0(1:n(), ".", n())) %>% gather(Day, Counts, -CloneSize) %>%   
    group_by(Day) %>% mutate(TotalCounts = sum(Counts), 
                             value       = qbeta(  0.5, Counts, TotalCounts-Counts), 
                             low_lim     = qbeta(0.025, Counts, TotalCounts-Counts), 
                             hi_lim      = qbeta(0.975, Counts, TotalCounts-Counts)) %>%  
    mutate(Age = as.numeric(str_replace(Day, pattern = "T", replacement = ""))) %>% ungroup() %>% 
    select(CloneSize, Age, value, low_lim, hi_lim)
  data_x  
}

#' Adjust theory drift data for plotting
format_th_data = function(th_x, time_points)
{
  fractions_split = nrow(th_x)
  col_names       = c("Age", paste0(1:fractions_split, ".", fractions_split))
  ## Sim process  
  model_preds           = data.frame(Age = time_points, t(th_x))
  colnames(model_preds) = col_names
  model_list            = model_preds %>% gather(CloneSize, value, -Age)
  model_list
}


