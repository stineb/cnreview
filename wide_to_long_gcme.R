wide_to_long_gcme <- function( df_wide, keyvars ){
  
  require(stringr)
  
  ## if level is 'elevated', turn all factor levels that are part of 'treatment' to TRUE
  switch_factor <- function(df){
    if (df$level=="elevated"){
      for (ifactor in strsplit(df$treatment, NULL) %>% unlist()){
        df[[ifactor]] <- TRUE
      }
    }
    return(df)
  }
  
  ## First get available factors
  factors_avl <- df_wide %>%
    filter(!(factors %in% c("x","x99"))) %>%      # Something weird with this
    select(factors) %>% 
    unique() %>% 
    mutate( factors_sep = strsplit(factors, NULL) ) %>% 
    select(factors_sep) %>% 
    unlist() %>% 
    unname() %>% 
    unique()  

  df_long <- df_wide %>%
    
    ## gather 'mean' based on columns: 'ambient', 'elevated'
    select( -ambient_Sd, -elevated_Sd, -ambient_Se, -elevated_Se ) %>%
    tidyr::gather(level, mean, c(ambient, elevated)) %>%
    
    ## gather 'sd' based on columns: 'ambient_Sd', 'elevated_Sd'
    left_join( 
      select( df_wide, -ambient, -elevated, -ambient_Se, -elevated_Se ) %>%
        tidyr::gather(level, sd, c(ambient_Sd, elevated_Sd)) %>%
        mutate( level = stringr::str_split(.$level, "_") %>% purrr::map(., 1) %>% unlist() ),
      by = c(keyvars) ) %>% 

    ## gather 'se' based on columns: 'ambient_Se', 'elevated_Se'
    left_join( 
      select( df_wide, -ambient, -elevated, -ambient_Sd, -elevated_Sd ) %>%
        tidyr::gather(level, se, c(ambient_Se, elevated_Se)) %>%
        mutate( level = stringr::str_split(.$level, "_") %>% purrr::map(., 1) %>% unlist() ),
      by = c(keyvars) ) %>% 

    ## magically add columns corresponding to all available levels of 'factors'
    `is.na<-`(factors_avl) %>% 
    
    ## fill all new factor columns with FALSE
    dplyr::mutate_at( factors_avl, ~FALSE ) 
  
  ## Determine new factor columns based on information in 'treatment' and 'level'
  df_long <- purrr::map_dfr(as.list(1:nrow(df_long)), ~switch_factor(df_long[.,])) %>% 
    
    ## abandon columns, now obsolete with boolean columns for factors 
    select(-treatment, -factors, -level, -id) %>%
    
    # use only distinct columns (remember that the "absolute ambient" was a three-times repeated entry in the wide table)
    distinct()
  
  return(df_long)
}