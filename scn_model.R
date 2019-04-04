scn_model <- function( ctot0, csoil0, ppfd, lue, n_in, par, settings, method="scn", accelerate=FALSE ){

  f_supply <- function( cbg, n0, f_unavoidable=0.0, kr ){
    (1.0 - f_unavoidable) * n0 * cbg / (cbg + kr )
  }
  prod <- function(leafarea, ppfd, lue, kl ){
    ppfd * lue * leafarea / (leafarea + kl )
  }
  f_ndemand <- function( ..., r_cton_plant ){
    (1/r_cton_plant) * prod(...)
  }
  
  # leaf area and belowground C as a function of fraction in shoots (alpha)
  calc_leafarea_alpha <- function(alpha, ctot, sla){
    alpha * sla * ctot
  }
  calc_cbg_alpha <- function(alpha, ctot){
    (1-alpha) * ctot
  }
  
  setzero_alpha <- function(alpha, ctot, n0, ppfd, lue, r_cton_plant, sla, kr, kl){
    nsupply <- f_supply(  calc_cbg_alpha(      alpha, ctot ), n0=n0, kr=kr )
    ndemand <- f_ndemand( calc_leafarea_alpha( alpha, ctot, sla ), ppfd=ppfd, lue=lue, r_cton_plant=r_cton_plant, kl=kl )
    out <- nsupply - ndemand
    return(out)
  }

  ## to simplify life
  r_ntoc_plant <- 1/par$r_cton_plant
  r_ntoc_soil  <- 1/par$r_cton_soil

  ## initialise output variables
  ntout <- ifelse( settings$out_spinup, settings$spinupyears + settings$nyeartrend, settings$nyeartrend)
  out_cplant_ag   <- c()
  out_nplant_ag   <- c()
  out_cplant_bg   <- c()
  out_nplant_bg   <- c()
  out_csoil       <- c()
  out_nsoil       <- c()
  out_clabl       <- c()
  out_nlabl       <- c()
  out_nloss       <- c()
  out_netmin      <- c()
  out_nup         <- c()
  out_npp         <- c()
  out_clitterfall <- c()
  out_nlitterfall <- c()
  
  ## plant, aboveground, start with root:shoot ratio of 0.5
  cplant_ag <- ctot0 * 0.5
  nplant_ag <- cplant_ag * r_ntoc_plant

  ## plant, belowground
  cplant_bg <- ctot0 * 0.5
  nplant_bg <- cplant_bg * r_ntoc_plant

  ## soil
  csoil <- csoil0
  nsoil <- csoil * r_ntoc_plant * par$tau_soil_n/par$tau_soil_c

  ## plant labile pools
  nlabl <- 0
  clabl <- 0

  ## to avoid numerical oscillation
  calc_turnover <- function( c0, tau, dt = 1.0 ){
    c0 * (1.0 - exp(-1/tau * dt ) )
    # c0/tau
  }
  
  if (method=="conly"){
    if (par$tau_plant^-1*par$kl>par$alpha_fix*par$sla*ppfd*lue){
      print("WARNING: PLANT WILL DIE.")
    }
  }

  spinup <- TRUE

  ##----------------------------------------------BEGIN OF LOOP
  itout <- 0
  for (it in 1:(settings$spinupyears + settings$nyeartrend)){
    
    if (it==settings$spinupyears+1) spinup <- FALSE
    if (settings$out_spinup){
      itout <- itout + 1
    } else {
      if (!spinup) itout <- itout + 1
    }

    ## MANIPULATION ----------------
    # if (it>(settings$spinupyears+100)) lue <- 1.1
    ##------------------------------    
    
    itin <- max(it-settings$spinupyears, 1)
    my_ppfd <- ifelse( length(ppfd) == settings$nyeartrend, ppfd[itin], ppfd[1] )
    my_lue  <- ifelse( length(lue ) == settings$nyeartrend, lue[itin],  lue[1]  )
    my_n_in <- ifelse( length(n_in) == settings$nyeartrend, n_in[itin], n_in[1] )
    
    ## Soil turnover
    csoil_turnover <- calc_turnover( csoil, par$tau_soil_c )
    csoil          <- csoil - csoil_turnover
    
    nsoil_turnover <- calc_turnover( nsoil, par$tau_soil_n )
    nsoil          <- nsoil - nsoil_turnover
    
    ## Net mineralisation
    netmin         <- nsoil_turnover + my_n_in
    
    ## Plant turnover, needs to be after acquisition and before new balance evaluation
    ctot <- cplant_ag + cplant_bg + clabl
    ntot <- nplant_ag + nplant_bg + nlabl
    # print( paste( "C:N ratio of plant before turnover:", ctot/ntot ) )
    
    cturnover_ag <- calc_turnover( cplant_ag, par$tau_plant )
    cplant_ag    <- cplant_ag - cturnover_ag
    
    nturnover_ag <- calc_turnover( nplant_ag, par$tau_plant)
    nplant_ag    <- nplant_ag - nturnover_ag
    
    cturnover_bg <- calc_turnover( cplant_bg, par$tau_plant )
    cplant_bg    <- cplant_bg - cturnover_bg
    
    nturnover_bg <- calc_turnover( nplant_bg, par$tau_plant)
    nplant_bg    <- nplant_bg - nturnover_bg
    
    c_litterfall <- cturnover_ag + cturnover_bg
    n_litterfall <- nturnover_ag + nturnover_bg
    
    # print(paste("C:N of litterfall", c_litterfall/n_litterfall))
    
    csoil <- csoil + c_litterfall
    nsoil <- nsoil + n_litterfall

    # if (spinup & accelerate & it == (settings$spinupyears-2700)){
    #   csoil <- c_litterfall * par$tau_soil_c
    #   nsoil <- n_litterfall * par$tau_soil_n
    # }
    
    ## update total biomass with (allocatable) labile C and N 
    ctot <- cplant_ag + cplant_bg + clabl
    ntot <- nplant_ag + nplant_bg + nlabl
    # print( paste( "C:N ratio of plant after turnover:", ctot/ntot ) )
    
    ## Get balanced allocation
    if (method=="conly"){
      
      root <- par$alpha_fix

      ## redesign the plant (immediate par$effect assumption)
      clabl <- 0
      nlabl <- 0
      cplant_ag <- root * ctot
      aleaf <- par$sla * cplant_ag
      cplant_bg <- (1 - root) * ctot
      
      nplant_ag <- cplant_ag * r_ntoc_plant
      nplant_bg <- cplant_bg * r_ntoc_plant
      
      clabl <- prod( aleaf, ppfd=my_ppfd, lue=my_lue, kl=par$kl ) #+clabl
      nlabl <- r_ntoc_plant * clabl
      # print( paste( "C:N ratio of labile:", clabl/nlabl ) )
      
    } else if (method=="scn") {

      if (spinup & accelerate & it < (settings$spinupyears-2000) ){
        ## short-cut: by-passing soil
        #nlabl <- f_noloss( cplant_bg ) * r_ntoc_plant * (cplant_bg + cplant_ag) / tau_plant  #+ nlabl
        netmin <- my_n_in + (nplant_ag + nplant_bg) / par$tau_plant
      }
      
      root <- uniroot( 
        function(x) 
          setzero_alpha( x, ctot, netmin, my_ppfd, my_lue, par$r_cton_plant, par$sla, par$kr, par$kl ), 
        interval=c(0,1) 
      )$root
    
      ## redesign the plant (immediate par$effect assumption)
      clabl <- 0
      nlabl <- 0
      cplant_ag <- root * ctot
      aleaf <- par$sla * cplant_ag
      cplant_bg <- (1 - root) * ctot
      
      nplant_ag <- cplant_ag * r_ntoc_plant
      nplant_bg <- cplant_bg * r_ntoc_plant
      
      nlabl <- f_supply( cplant_bg, n0=netmin, kr=par$kr ) #+ nlabl
      clabl <- prod( aleaf, ppfd=my_ppfd, lue=my_lue, kl=par$kl ) #+clabl
      # print( paste( "C:N ratio of labile:", clabl/nlabl ) )
      
    }
    
    ## gather output variables
    # print(paste("itout: ", itout))
    if ( itout > 0 ){
      out_cplant_ag[itout]   <- cplant_ag
      out_nplant_ag[itout]   <- nplant_ag
      out_cplant_bg[itout]   <- cplant_bg
      out_nplant_bg[itout]   <- nplant_bg
      out_csoil[itout]       <- csoil
      out_nsoil[itout]       <- nsoil
      out_clabl[itout]       <- clabl
      out_nlabl[itout]       <- nlabl
      out_nloss[itout]       <- netmin - nlabl
      out_netmin[itout]      <- netmin
      out_nup[itout]         <- nlabl
      out_npp[itout]         <- clabl
      out_clitterfall[itout] <- c_litterfall
      out_nlitterfall[itout] <- n_litterfall
    }
    
  }
  ##----------------------------------------------END OF LOOP

  df_out <- tibble(
    simyear      = 1:ntout,
    cplant_ag    = out_cplant_ag,
    nplant_ag    = out_nplant_ag,
    cplant_bg    = out_cplant_bg,
    nplant_bg    = out_nplant_bg,
    csoil        = out_csoil,
    nsoil        = out_nsoil,
    clabl        = out_clabl,
    nlabl        = out_nlabl,
    nloss        = out_nloss,
    netmin       = out_netmin,
    nup          = out_nup,
    npp          = out_npp,
    c_litterfall = out_clitterfall,
    n_litterfall = out_nlitterfall
  )

  return(df_out)
}
