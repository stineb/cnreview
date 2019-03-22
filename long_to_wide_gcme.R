long_to_wide_gcme <- function(df_long, keyvar){
  
  joinvars <- names(df_long)[-which(names(df_long) %in% c("level", "mean", "sd", "se", keyvar))]
  factors_all <- df_long %>%
    filter(!(factors %in% c("x","x99"))) %>%      # Something weird with this
    select(factors) %>% 
    unique() %>% 
    mutate( factors_sep = strsplit(factors, NULL) ) %>% 
    select(factors_sep) %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  factorvars <- factors_all[-which(factors_all==keyvar)]
  
  # take all data where factor 'c' is TRUE
  df_wide <- dplyr::filter(df_long, eval(parse(text=keyvar)) ) %>% 
    
    # call this 'elevated'
    dplyr::rename(elevated=mean, elevated_sd=sd, elevated_se=se) %>% 
    
    # remove column 'c', is no longer used
    dplyr::select(-keyvar, -level) %>% 
    
    # merge this with the corresponding row where all other factors are the same, while 'c' is FALSE
    left_join( 
      
      # take all data where factor 'c' is FALSE
      dplyr::filter(df_long, !(eval(parse(text=keyvar)))) %>% 
        
        # call this 'ambient'
        dplyr::rename(ambient=mean, ambient_sd=sd, ambient_se=se) %>% 
        
        # remove column 'c', is no longer used
        dplyr::select(-keyvar, -level),
      
      ## merge by all other columns
      by = joinvars ) %>% 
    
    # order columns
    select( c(joinvars, ambient, elevated, ambient_se, elevated_se, ambient_sd, elevated_sd)  )
  
  return(df_wide)
}