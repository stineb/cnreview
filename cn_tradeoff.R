prod <- function( cleaf, ppfd, lue ){
  kbeer         = 0.5
  sla           = 0.1

  prod <- ppfd * lue * ( 1.0 - exp( - kbeer * sla * cleaf ) )

  return(prod)
}

f_noloss <- function( croot ){
  f_unavoidable = 0.1
  k_noloss      = 100

  # f_avl <- (1.0 - f_unavoidable) * croot / (croot + 100 )
  # f_avl <- min( 1.0, croot / 30 )
  out   <- (1.0 - f_unavoidable) * (1.0 - exp( - 0.01 * croot ) )
  
  return(out)
}

eval_cnbalance <- function( f_ag, ctot, navl, lue, ppfd ){

  eff  = 0.7
  resp = 0.1
  r_cton_plant = 30

  out <- (( eff  * ( prod( f_leaf_of_c_ag_tot( f_ag * ctot ) * f_ag * ctot, ppfd, lue ) - resp * ctot ) ) /
         ( navl * f_noloss( (1.0 - f_ag) * ctot ) ) ) - r_cton_plant

  return( out )
}

c_wood <- function( ctot ){
  par = 50
  sqrt( ctot^2 + par ) - sqrt(par)
}

c_leaf <- function( ctot ){
  kl = 0.3
  c_leaf_max = 20
  c_leaf_max * (1 - exp(-kl*ctot))
}

f_leaf_of_c_ag_tot <- function( c_ag_tot ){
  ifelse( c_ag_tot==0, 1.0, c_leaf( c_ag_tot ) / (c_leaf( c_ag_tot ) + c_wood( c_ag_tot )) )
}

## Static CN-model for one state-------------------------------
lue  = 1
ppfd = 100
navl = 0.8
ctot = 350

curve( c_wood(x), from = 0, to = 50 )
curve( c_leaf(x), from = 0, to = 50, add = TRUE )
curve( f_leaf_of_c_ag_tot(x), from = 0, to = ctot )

curve( prod( f_leaf_of_c_ag_tot( x ) * x, ppfd, lue ), from = 0, to = ctot )
curve( f_noloss, from = 0, to = ctot*2 )

## Visualise imbalance curve
par(mfrow=c(1,1))
curve( eval_cnbalance( x, ctot, navl, lue, ppfd ), from = 0, to = 1.0, col="red", ylim=c(-20,100) )
abline( h=0, lty=3 )
out_root <- uniroot( function(x) eval_cnbalance( x, ctot, navl, lue, ppfd ), interval=c(0,1.0) )
abline( v = out_root$root, lty = 2 )

## production curve
curve( prod( f_leaf_of_c_ag_tot( x ) * x, ppfd, lue ), from = 0, to = ctot ) # x is total aboveground
points( out_root$root * ctot, prod( f_leaf_of_c_ag_tot( out_root$root * ctot ) * out_root$root * ctot, ppfd, lue ), pch=16, col="red" )

## N uptake curve
curve( navl * f_noloss( x ), from = 0, to = ctot ) # x is total belowground
points( (1.0 - out_root$root) * ctot, navl * f_noloss( (1.0 - out_root$root) * ctot ), pch=16, col="red" )


## Dynamic CN-model--------------------------------------------
lue          = 1       # change this to two and the whole thing collapses. works best for 1. 
ppfd         = 100
r_ntoc_plant = 1/30
r_ntoc_soil  = 1/10
eff          = 0.7
resp         = 0.1
tau_plant    = 10
tau_soil_c   = 50
tau_soil_n   = 150
tau_labl     = 0.5
n_input      = 1.0

## Simulation settings
ntsteps <- 3000

## initialise output variables
out_cplant_ag <- rep(NA, ntsteps)
out_nplant_ag <- rep(NA, ntsteps)
out_cplant_bg <- rep(NA, ntsteps)
out_nplant_bg <- rep(NA, ntsteps)
out_csoil     <- rep(NA, ntsteps)
out_nsoil     <- rep(NA, ntsteps)
out_clabl     <- rep(NA, ntsteps)
out_nlabl     <- rep(NA, ntsteps)

## plant, aboveground
cplant_ag <- 150
nplant_ag <- cplant_ag * r_ntoc_plant

## plant, belowground
cplant_bg <- 150
nplant_bg <- cplant_bg * r_ntoc_plant

## soil
csoil <- 1000
nsoil <- csoil * r_ntoc_plant * tau_soil_n/tau_soil_c

## plant labile pools
nlabl <- 0
clabl <- 0

for (it in 1:ntsteps){

  ## Soil turnover and net mineralisation
  csoil          <- csoil - csoil / tau_soil_c
  nsoil_turnover <- nsoil / tau_soil_n 
  nsoil          <- nsoil - nsoil_turnover
  netmin         <- nsoil_turnover + n_input
  
  ## actual C and N acquisition with current plant
  clabl <- prod( f_leaf_of_c_ag_tot( cplant_ag ) * cplant_ag, ppfd, lue ) #+ clabl
  nlabl <- f_noloss( cplant_bg ) * netmin #+ nlabl
  
  # ## short-cut: by-passing soil
  # nlabl <- f_noloss( cplant_bg ) * r_ntoc_plant * (cplant_bg + cplant_ag) / tau_plant  #+ nlabl

  ## limit allocatable C depending on labile N
  c_alloc <- min( nlabl / (eff * r_ntoc_plant), clabl )

  # ##----------
  # netmin <- 0.8 # xxx test
  # ##----------

  ## Get optimal fraction to leaves
  out_root <- uniroot( function(x) eval_cnbalance( x, ctot, netmin, lue, ppfd ), interval=c(0,1.0) )
  fleaf <- out_root$root
  print(fleaf)

  # ##----------
  # ## Visualise imbalance curve
  # ctot <- (cplant_ag + cplant_bg + eff * c_alloc)
  # par(mfrow=c(1,1))
  # curve( eval_cnbalance( x, ctot, netmin, lue, ppfd ), from = 0, to = 1.0, col="red", ylim=c(-20,100) )
  # abline( h=0, lty=3 )
  # abline( v = fleaf, lty = 2 )

  # ## production curve
  # curve( prod( f_leaf_of_c_ag_tot( x ) * x, ppfd, lue ), from = 0, to = ctot ) # x is total aboveground
  # points( fleaf * ctot, prod( f_leaf_of_c_ag_tot( fleaf * ctot ) * fleaf * ctot, ppfd, lue ), pch=16, col="red" )

  # ## N uptake curve
  # curve( netmin * f_noloss( x ), from = 0, to = ctot ) # x is total belowground
  # points( (1.0 - fleaf) * ctot, netmin * f_noloss( (1.0 - fleaf) * ctot ), pch=16, col="red" )
  # ##----------

  ## Additional check only necessary for legacy allocation setup (Version A) 
  # if ( eval_cnbalance( 0, ctot, netmin, lue, ppfd ) > 0.0 ){
  #   ## all to roots
  #   fleaf <- 0.0
  # } else if ( eval_cnbalance( 1, ctot, netmin, lue, ppfd ) < 0.0 ){
  #   ## all to leaves
  #   fleaf <- 1.0
  # } else {
  #   ## find allocation fraction to leaves so that return next year matches (eff * r_cton_plant)
  #   out_root <- uniroot( function(x) eval_cnbalance( x, (cplant_ag + cplant_bg + eff * c_alloc), netmin, lue, ppfd ), interval=c(1e-12,1.0) )
  #   fleaf <- out_root$root
  # }

  ## update pools
  clabl <- clabl - c_alloc
  nlabl <- nlabl - eff * c_alloc * r_ntoc_plant

  ## labile pool decay
  clabl <- clabl - tau_labl * clabl
  nlabl <- nlabl - tau_labl * nlabl

  ##---------------------------------------------------------------------  
  # ## VERSION A: fleaf only acts on new growth
  # ## allocate C, N follows from r_ntoc_plant
  # cturnover_ag <- cplant_ag / tau_plant
  # c_alloc_ag   <- fleaf * c_alloc * eff
  # cplant_ag    <- cplant_ag + c_alloc_ag - cturnover_ag

  # nturnover_ag <- nplant_ag / tau_plant
  # n_alloc_ag   <- c_alloc_ag * r_ntoc_plant
  # nplant_ag    <- nplant_ag + n_alloc_ag - nturnover_ag

  # cturnover_bg <- cplant_bg / tau_plant
  # c_alloc_bg   <- (1.0 - fleaf) * c_alloc * eff
  # cplant_bg    <- cplant_bg + c_alloc_bg - cturnover_bg

  # nturnover_bg <- nplant_bg / tau_plant
  # n_alloc_bg   <- c_alloc_bg * r_ntoc_plant
  # nplant_bg    <- nplant_bg + n_alloc_bg - nturnover_bg
  ##---------------------------------------------------------------------  


  ##---------------------------------------------------------------------  
  ## VERSION B: fleaf redistributes entire plant (to avoid legacy in mal-allocation)
  ## allocate C, N follows from r_ntoc_plant
  cturnover_ag <- cplant_ag / tau_plant
  cplant_ag    <- cplant_ag - cturnover_ag
  cplant_ag    <- fleaf * (c_alloc * eff + cplant_ag + cplant_bg)
  
  nturnover_ag <- nplant_ag / tau_plant
  nplant_ag    <- nplant_ag - nturnover_ag
  nplant_ag    <- cplant_ag * r_ntoc_plant
  
  cturnover_bg <- cplant_bg / tau_plant
  cplant_bg    <- cplant_bg - cturnover_bg
  cplant_bg    <- (1.0 - fleaf) * (c_alloc * eff + cplant_ag + cplant_bg)
  
  nturnover_bg <- nplant_ag / tau_plant
  nplant_ag    <- nplant_ag - nturnover_bg
  nplant_ag    <- cplant_ag * r_ntoc_plant
  ##---------------------------------------------------------------------  

  csoil <- csoil + cturnover_ag + cturnover_bg
  nsoil <- nsoil + nturnover_ag + nturnover_bg

  ## gater output variables
  out_cplant_ag[it] <- cplant_ag 
  out_nplant_ag[it] <- nplant_ag 
  out_cplant_bg[it] <- cplant_bg 
  out_nplant_bg[it] <- nplant_bg 
  out_csoil[it]     <- csoil     
  out_nsoil[it]     <- nsoil     
  out_clabl[it]     <- clabl     
  out_nlabl[it]     <- nlabl     

}
##---------------------------------------------------------
plot( 1:ntsteps, out_cplant_ag )
plot( 1:ntsteps, out_cplant_bg )
plot( 1:ntsteps, out_csoil )
plot( 1:ntsteps, out_nsoil )
plot( 1:ntsteps, out_clabl )
plot( 1:ntsteps, out_nlabl )
plot( 1:ntsteps, out_cplant_ag/(out_cplant_bg+out_cplant_ag) )
