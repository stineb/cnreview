check_cf <- function(df_wide, df_long, verbose=FALSE){
  
  ## check cf experiments
  df_cf <- df_wide %>% filter(treatment %in% c("f","c","cf","fc"))
  
  out <- tibble()
  
  all_exp <- df_cf$exp_nam %>% unique
  
  for (iexp in all_exp){
    
    df_sub <- df_cf %>% filter(exp_nam==iexp)
    
    all_var <- df_sub$Data_type %>% unique
    
    for (ivar in all_var){
      
      df_subsub <- df_sub %>% filter(Data_type==ivar)
      
      all_var_unit <- df_subsub$Unit %>% unique
      
      for (iunit in all_var_unit){
        
        df_subsubsub <- df_subsub %>% filter(Unit==iunit)
        
        # elevated-f treatment
        correct_f <- df_subsubsub %>% 
          filter(treatment=="f") %>% 
          select(exp_nam, Data_type, Unit, Year, mean=elevated, se=elevated_Se, sd=elevated_Sd) %>%
          all_equal( ., df_long %>% filter(f & !c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        # elevated-c treatment
        correct_c <- df_subsubsub %>% 
          filter(treatment=="c") %>% 
          select(exp_nam, Data_type, Unit, Year, mean=elevated, se=elevated_Se, sd=elevated_Sd) %>%
          all_equal( ., df_long %>% filter(!f & c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        # elevated-cf treatment
        correct_cf <- df_subsubsub %>% 
          filter(treatment=="cf") %>% 
          select(exp_nam, Data_type, Unit, Year, mean=elevated, se=elevated_Se, sd=elevated_Sd) %>%
          all_equal( ., df_long %>% filter(f & c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        # (absolute) control: WARNING: DUPLICATED ROWS, DIFFERING ONLY BY 'ALIAS'
        correct_0 <- df_subsubsub %>% 
          distinct(exp_nam, Data_type, Unit, Year, ambient, ambient_Se, ambient_Sd, .keep_all=TRUE) %>% 
          select(exp_nam, Data_type, Unit, Year, mean=ambient, se=ambient_Se, sd=ambient_Sd) %>%
          all_equal( ., df_long %>% filter(!f & !c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        if (verbose){
          if (!identical(correct_0, TRUE)){
            rlang::warn(correct_0)
            print("WIDE:")
            df_subsubsub %>%
              distinct(exp_nam, Data_type, Unit, Year, ambient, ambient_Se, ambient_Sd, .keep_all=TRUE) %>%
              select(exp_nam, Data_type, Unit, Year, mean=ambient, se=ambient_Se, sd=ambient_Sd) %>%
              print()
            print("LONG:")
            df_long %>%
              filter(!f & !c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>%
              select(exp_nam, Data_type, Unit, Year, mean, se, sd) %>%
              print()
            print("---------")
          }
        }
        
        out <- bind_rows(out, c(factors="cf", exp_nam=iexp, Data_type=ivar, Unit=iunit, correct_0=correct_0, correct_c=correct_c, correct_f=correct_f, correct_cf=correct_cf))
        
      }
    }
  }
  
  return(out)
}