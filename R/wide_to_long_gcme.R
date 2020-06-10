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
    dplyr::filter(!(factors %in% c("x","x99"))) %>%      # Something weird with this
    dplyr::select(factors) %>% 
    unique() %>% 
    mutate( factors_sep = strsplit(factors, NULL) ) %>% 
    dplyr::select(factors_sep) %>% 
    unlist() %>% 
    unname() %>% 
    unique()  
  
  remove_suffix_sd <- function(string){stringr::str_remove(string, "_Sd")}
  remove_suffix_se <- function(string){stringr::str_remove(string, "_Se")}
  
  if (!("level" %in% keyvars)) keyvars <- c(keyvars, "level")
  
  df_long <- df_wide %>%
    
    ## gather 'mean' based on columns: 'ambient', 'elevated'
    dplyr::select( -ambient_Sd, -elevated_Sd, -ambient_Se, -elevated_Se ) %>%
    tidyr::pivot_longer(c(ambient, elevated), names_to = "level", values_to = "mean" ) %>%
    
    ## gather 'sd' based on columns: 'ambient_Sd', 'elevated_Sd'
    left_join( 
      select( df_wide, -ambient, -elevated, -ambient_Se, -elevated_Se ) %>%
        tidyr::pivot_longer(c(ambient_Sd, elevated_Sd), names_to = "level", values_to = "sd", names_transform = list(level = remove_suffix_sd)),
        # mutate( level = stringr::str_split(.$level, "_") %>% purrr::map(., 1) %>% unlist() ),
      by = keyvars ) %>% 
    
    ## gather 'se' based on columns: 'ambient_Se', 'elevated_Se'
    left_join( 
      dplyr::select( df_wide, -ambient, -elevated, -ambient_Sd, -elevated_Sd ) %>%
        tidyr::pivot_longer(c(ambient_Se, elevated_Se), names_to = "level", values_to = "se", names_transform = list(level = remove_suffix_se)),
      by = keyvars ) %>% 
    
    ## magically add columns corresponding to all available levels of 'factors'
    `is.na<-`(factors_avl) %>% 
    
    ## fill all new factor columns with FALSE
    dplyr::mutate_at( factors_avl, ~FALSE ) 
  
  ## Determine new factor columns based on information in 'treatment' and 'level'
  df_long <- purrr::map_dfr(as.list(1:nrow(df_long)), ~switch_factor(df_long[.,]))
  
  ## Now we have duplicate ambients. Sort it out: only difference is in "treatment"
  kevars_red <- c(keyvars[-which(keyvars %in% c("treatment", "id"))], factors_avl)
  df_long <- df_long %>% 
    distinct(across(kevars_red), .keep_all = TRUE) %>% 
    
    ## abandon columns, now obsolete with boolean columns for factors 
    dplyr::select(-treatment, -level, -id)
  
  return(df_long)
}