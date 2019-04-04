scn_model <- function( ctot0, csoil0, ppfd, lue, n_in, par, settings, method="scn" ){

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
  out_cplant_ag <- rep(NA, settings$ntsteps)
  out_nplant_ag <- rep(NA, settings$ntsteps)
  out_cplant_bg <- rep(NA, settings$ntsteps)
  out_nplant_bg <- rep(NA, settings$ntsteps)
  out_csoil     <- rep(NA, settings$ntsteps)
  out_nsoil     <- rep(NA, settings$ntsteps)
  out_clabl     <- rep(NA, settings$ntsteps)
  out_nlabl     <- rep(NA, settings$ntsteps)
  out_nloss     <- rep(NA, settings$ntsteps)
  out_netmin    <- rep(NA, settings$ntsteps)
  out_clitterfall <- rep(NA, settings$ntsteps)
  out_nlitterfall <- rep(NA, settings$ntsteps)
  
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

  for (it in 1:settings$ntsteps){
    
    # ## step increase in light use par$efficiency ~ CO2 fertilisation
    # if (it==1000) lue <- 1.2 * lue
    
    ## Soil turnover
    csoil_turnover <- calc_turnover( csoil, par$tau_soil_c )
    csoil          <- csoil - csoil_turnover
    
    nsoil_turnover <- calc_turnover( nsoil, par$tau_soil_n )
    nsoil          <- nsoil - nsoil_turnover
    
    ## Net mineralisation
    netmin         <- nsoil_turnover + n_in
    
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
    
    print(paste("C:N of litterfall", c_litterfall/n_litterfall))
    
    csoil <- csoil + c_litterfall
    nsoil <- nsoil + n_litterfall
    
    ## update total biomass with (allocatable) labile C and N 
    ctot <- cplant_ag + cplant_bg + clabl
    ntot <- nplant_ag + nplant_bg + nlabl
    # print( paste( "C:N ratio of plant after turnover:", ctot/ntot ) )
    
    ## Get balanced allocation
    if (method=="conly"){
      
      root <- 0.5

      ## redesign the plant (immediate par$effect assumption)
      clabl <- 0
      nlabl <- 0
      cplant_ag <- root * ctot
      aleaf <- par$sla * cplant_ag
      cplant_bg <- (1 - root) * ctot
      
      nplant_ag <- cplant_ag * r_ntoc_plant
      nplant_bg <- cplant_bg * r_ntoc_plant
      
      clabl <- prod( aleaf, ppfd=ppfd, lue=lue, kl=par$kl ) #+clabl
      nlabl <- r_ntoc_plant * clabl
      # print( paste( "C:N ratio of labile:", clabl/nlabl ) )
      
    } else if (method=="scn") {
      
      root <- uniroot( 
        function(x) 
          setzero_alpha( x, ctot, netmin, ppfd, lue, par$r_cton_plant, par$sla, par$kr, par$kl ), 
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
      clabl <- prod( aleaf, ppfd=ppfd, lue=lue, kl=par$kl ) #+clabl
      # print( paste( "C:N ratio of labile:", clabl/nlabl ) )
      
    }
    
    ## gater output variables
    out_cplant_ag[it] <- cplant_ag
    out_nplant_ag[it] <- nplant_ag
    out_cplant_bg[it] <- cplant_bg
    out_nplant_bg[it] <- nplant_bg
    out_csoil[it]     <- csoil
    out_nsoil[it]     <- nsoil
    out_clabl[it]     <- clabl
    out_nlabl[it]     <- nlabl
    out_nloss[it]     <- netmin - nlabl
    out_netmin[it]    <- netmin
    out_clitterfall[it] <- c_litterfall
    out_nlitterfall[it] <- n_litterfall
    
  }
  ##----------------------------------------------END OF LOOP

  df_out <- tibble(
    simyear   = 1:settings$ntsteps,
    cplant_ag = out_cplant_ag,
    nplant_ag = out_nplant_ag,
    cplant_bg = out_cplant_bg,
    nplant_bg = out_nplant_bg,
    csoil     = out_csoil,
    nsoil     = out_nsoil,
    clabl     = out_clabl,
    nlabl     = out_nlabl,
    nloss     = out_nloss,
    netmin    = out_netmin,
    c_litterfall = out_clitterfall,
    n_litterfall = out_nlitterfall
  )

  return(df_out)
}

## Simulation settings
settings <- list(
  ntsteps = 5000
  )

## model parameters
par <- list(
  r_cton_plant = 30,
  r_cton_soil  = 10,
  tau_plant    = 10,
  tau_soil_c   = 50,
  tau_soil_n   = 150,
  tau_labl     = 0.5,
  sla          = 0.1,
  eff          = 1.0,
  kl           = 50,
  kr           = 50
)

## Environmental conditions
ppfd <- 100
lue  <- 1
n_in <- 0.3

## Run the model
df_scn <- scn_model( ctot0=10, csoil0=100, 
                     ppfd=ppfd, lue=lue, n_in=n_in, 
                     par=par, settings=settings, method="conly" 
                     )

library(ggplot2)
df_scn %>%
  tidyr::gather(varnam, value, c(cplant_ag, cplant_bg)) %>% 
  ggplot( aes(x=simyear, y=value, color=varnam)) +
  geom_line() +
  labs(title="Plant C", x="Simulation Year", y=expression(paste("C pool (g C m"^{-2}, ")")))

df_scn %>% 
  ggplot( aes(x=simyear, y=csoil) ) +
  geom_line() +
  labs(title="Soil C", x="Simulation Year", y=expression(paste("C pool (g C m"^{-2}, ")")))

df_scn %>% 
  ggplot( aes(x=simyear, y=nsoil) ) +
  geom_line() +
  labs(title="Soil N", x="Simulation Year", y=expression(paste("N pool (g N m"^{-2}, ")")))

df_scn %>% 
  ggplot( aes(x=simyear, y=clabl) ) +
  geom_line() +
  labs(title="Labile C", x="Simulation Year", y=expression(paste("C pool (g C m"^{-2}, ")")))

df_scn %>% 
  ggplot( aes(x=simyear, y=clabl/nlabl) ) +
  geom_line() +
  labs(title="Labile C:N", x="Simulation Year", y=expression(paste("C:N ratio (g C g N"^{-1}, ")"))) +
  geom_hline(yintercept=par$r_cton_plant, linetype="dotted")

soil_cton_expected <- par$r_cton_plant * par$tau_soil_c / par$tau_soil_n
df_scn %>% 
  ggplot( aes(x=simyear, y=csoil/nsoil) ) +
  geom_line() +
  labs(title="Soil C:N", x="Simulation Year", y=expression(paste("C:N ratio (g C g N"^{-1}, ")"))) +
  geom_hline(yintercept=soil_cton_expected, linetype="dotted")

df_scn %>% 
  ggplot( aes(x=simyear, y=c_litterfall/n_litterfall) ) +
  geom_line() +
  labs(title="Litterfall C:N", x="Simulation Year", y=expression(paste("C:N ratio (g C g N"^{-1}, ")"))) +
  geom_hline(yintercept=par$r_cton_plant, linetype="dotted")

df_scn %>% 
  ggplot( aes(x=simyear, y=(cplant_ag+cplant_bg)/(nplant_ag+nplant_bg)) ) +
  geom_line() +
  labs(title="Plant C:N", x="Simulation Year", y=expression(paste("C:N ratio (g C g N"^{-1}, ")")))

df_scn %>% 
  ggplot( aes(x=simyear, y=nloss) ) +
  geom_line() +
  labs(title="N losses", x="Simulation Year", y=expression(paste("N flux (g N m"^{-2}, " yr"^{-1}, ")"))) +
  geom_hline(yintercept = n_in, linetype="dashed")

# plot( 1:ntsteps, out_clabl, type = "l" )
# plot( 1:ntsteps, out_nlabl, type = "l" )
# plot( 1:ntsteps, out_cplant_ag/(out_cplant_bg+out_cplant_ag), type = "l" )

# print("Steady state f_ag (fraction of aboveground to total plant C):")
# print( out_cplant_ag[ntsteps]/(cplant_bg[ntsteps]+cplant_ag[ntsteps]) )

# print("Steady state total plant C:")
# print( cplant_bg[ntsteps]+cplant_ag[ntsteps] )

# print("Steady state net mineralisation:")
# print( netmin )

# print("Expected steady state net mineralisation:")
# print( n_in + (nplant_ag + nplant_bg) / par$tau_plant )
# # print( calc_turnover( nsoil, par$tau_soil_n ))
# print( (cplant_ag + cplant_bg) / (par$r_cton_plant * par$tau_plant) + n_in )

