agg_meta <- function(df, groupvar){
  out_meta <- df %>% dplyr::filter(my_varnam==eval(parse_character(groupvar))) %>% 
    metafor::rma( logr, logr_var, method = "DL", data = . )
  
  # transform back
  out_meta_scaled <- predict( out_meta, transf=exp )
  
  tibble(
    my_varnam=groupvar, 
    middle = out_meta$b[1,1], 
    ymin   = out_meta$ci.lb, 
    ymax   = out_meta$ci.ub,
    
    middle_scaled = out_meta_scaled$pred, 
    ymin_scaled   = out_meta_scaled$ci.lb, 
    ymax_scaled   = out_meta_scaled$ci.ub
  )
}