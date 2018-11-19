## Parameters
lue   <<- 1
ppfd   <<- 100
kbeer <<- 0.5
sla   <<- 1.0
fleaf_fixed <<- 0.5
r_ntoc_plant <<- 1/30
r_ntoc_soil <<- 1/10
eff <<- 0.7
resp <<- 0.1
f_unavoidable <<- 0.0
k_noloss <<- 0.5
tau_plant <<- 10
tau_soil_c <<- 50
tau_soil_n <<- 150
tau_labl <<- 0.5
n_input <<- 1.0

## Functions
prod <- function( cleaf ){
  prod <- ppfd * lue * ( 1.0 - exp( - kbeer * sla * cleaf ) )
  return(prod)
}
curve( prod, from = 0, to = 10 )

npp <- function( ctot ){
  npp <- eff * prod( ctot * fleaf_fixed ) - resp * ctot
  return(npp)
}

n_u_demand <- function( ctot ){
  n_u <- r_ntoc_plant * npp(ctot)
  return(n_u)
}

turnover <- function( ctot, tau ){
  turnover <- ctot / tau
  return(turnover)
}

f_noloss <- function( croot ){
  # out <- (1.0 - f_unavoidable) / (croot + k_noloss )
  out <- (1.0 - f_unavoidable) * (1.0 - exp( - k_noloss * croot ))
  return(out)
}
curve( f_noloss, from = 0, to = 10 )

n_u_supply <- function( ctot ){
  n_u <- f_noloss( ctot ) * (r_ntoc_plant * turnover( ctot ) + n_input)
  return(n_u)
}


## C-only model--------------------------------------------
## Simulation settings
ntsteps <- 1000

## starting value
cplant <- 20
csoil  <- 20

## output initialisation
out_cplant  <- rep(NA, ntsteps)
out_csoil   <- rep(NA, ntsteps)
out_nsoil   <- rep(NA, ntsteps)
out_netmin  <- rep(NA, ntsteps)

for (it in 1:ntsteps){
  
  ## Soil turnover and net mineralisation
  csoil_turnover <- csoil / tau_soil_c
  csoil          <- csoil - csoil_turnover
  nsoil_turnover <- nsoil / tau_soil_n 
  nsoil          <- nsoil - nsoil_turnover
  netmin         <- nsoil_turnover + n_input

  cplant <- cplant + npp( cplant ) - turnover( cplant, tau_plant )

  ## plant turnover
  cplant_turnover <- cplant / tau_plant
  nplant_turnover <- cplant_turnover * r_ntoc_plant

  csoil <- csoil + cplant_turnover
  nsoil <- nsoil + nplant_turnover

  ## output collection
  out_cplant[it] <- cplant
  out_csoil[it]  <- csoil
  out_nsoil[it]  <- nsoil
  out_netmin[it] <- netmin

}
plot( 1:ntsteps, out_cplant )
plot( 1:ntsteps, out_csoil )
plot( 1:ntsteps, out_nsoil )
plot( 1:ntsteps, out_netmin )
##---------------------------------------------------------




##---------------------------------------------------------
## CN MODEL BELOW
##---------------------------------------------------------





ntoc_balance <- function( fleaf, ctot ){

  ## net C gain (c_acq) and N acquisition
  c_acq <- prod( fleaf * ctot ) - resp * ctot
  n_acq <- f_noloss( (1.0 - fleaf) * ctot ) * (r_ntoc_plant * ctot / tau_plant + n_input)
  ntoc_acq <- n_acq / c_acq

  balance <- ntoc_acq - r_ntoc_plant * eff

  return(balance)
}

## doesn't have a root - always positive
curve( ntoc_balance(x, 10), from = 0, to = 1.0 )


eval_alloc <- function( fleaf, c_alloc, cplant_ag, cplant_bg, netmin ){

  # ## VERSION A: fleaf only acts on new growth
  # ## allocate C, N follows from r_ntoc_plant
  # c_alloc_ag <- fleaf * c_alloc
  # cplant_ag  <- cplant_ag + c_alloc_ag * eff - cplant_ag / tau_plant
  # 
  # c_alloc_bg <- (1.0 - fleaf) * c_alloc
  # cplant_bg  <- cplant_bg + c_alloc_bg * eff - cplant_bg / tau_plant

  ## VERSION B: fleaf redistributes entire plant (to avoid legacy in mal-allocation)
  ## total c to re-distribute:
  ctot <- (cplant_ag + cplant_bg) * (1 - 1/tau_plant) + eff * c_alloc

  ## create new plant
  cplant_ag <- fleaf * ctot
  cplant_bg <- (1.0 - fleaf) * ctot
  
  ## C and N acquisition after allocation
  c_acq <- prod( cplant_ag ) - resp * (cplant_ag + cplant_bg)
  n_acq <- f_noloss( cplant_bg ) * netmin

  eval <- eff * c_acq / n_acq - (1.0 / r_ntoc_plant)

  return(eval)

}

calc_r_cton_acq <- function( fleaf, c_alloc, cplant_ag, cplant_bg, netmin ){
    
  ## total c to re-distribute:
  ctot <- (cplant_ag + cplant_bg) * (1 - 1/tau_plant) + eff * c_alloc

  ## create new plant
  cplant_ag <- fleaf * ctot
  cplant_bg <- (1.0 - fleaf) * ctot
  
  ## C and N acquisition after allocation
  c_acq <- prod( cplant_ag ) - resp * (cplant_ag + cplant_bg)
  n_acq <- f_noloss( cplant_bg ) * netmin
  
  r_cton_acq <- eff * c_acq / n_acq
  
  return(r_cton_acq)
  
}

cplant_ag <- 175
cplant_bg <- 175
netmin <- 2.165898
c_alloc <- 2

out_root <- uniroot( function(x) eval_alloc( x, c_alloc, cplant_ag, cplant_bg, netmin), interval=c(0.0,1.0) )

curve( calc_r_cton_acq( x, c_alloc, cplant_ag, cplant_bg, netmin), from = 0, to = 1.0 )
abline(h=(1/r_ntoc_plant), lty=3)
abline(v=out_root$root, col="red")

curve(      eval_alloc( x, c_alloc, cplant_ag, cplant_bg, netmin), from = 0, to = 1.0 )
abline(h=0, lty=3)
abline(v=out_root$root, col="red")

calc_r_cton_acq( 0.99, 5, 20, 20, 0.5 )
calc_r_cton_acq( 0.01, 5, 20, 20, 0.5 )
print("moin he")

## CN-model--------------------------------------------

## Simulation settings
ntsteps <- 1000

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
cplant_ag <- 10
nplant_ag <- cplant_ag * r_ntoc_plant

## plant, belowground
cplant_bg <- 10
nplant_bg <- cplant_bg * r_ntoc_plant

## soil
csoil <- 100
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
  clabl <- prod( cplant_ag ) - resp * (cplant_ag + cplant_bg) #+ clabl
  nlabl <- f_noloss( cplant_bg ) * netmin #+ nlabl
  
  # ## short-cut: by-passing soil
  # nlabl <- f_noloss( cplant_bg ) * r_ntoc_plant * (cplant_bg + cplant_ag) / tau_plant  #+ nlabl

  ## limit allocatable C depending on labile N
  c_alloc <- min( nlabl / (eff * r_ntoc_plant), clabl )

  # curve( eval_alloc( x, c_alloc, cplant_ag, cplant_bg, netmin ), from = 0, to = 1.0 )

  curve( calc_r_cton_acq( x, c_alloc, cplant_ag, cplant_bg, netmin), from = 0, to = 1.0 )

  if ( eval_alloc( 0.0, c_alloc, cplant_ag, cplant_bg, netmin ) < 0.0 ){
    ## all to roots
    fleaf <- 0.0
  } else if ( eval_alloc( 1.0, c_alloc, cplant_ag, cplant_bg, netmin ) > 0.0 ){
    ## all to leaves
    fleaf <- 1.0
  } else {
    ## find allocation fraction to leaves so that return next year matches (eff * r_cton_plant)
    out_root <- uniroot( function(x) eval_alloc( x, c_alloc, cplant_ag, cplant_bg, netmin), interval=c(1e-12,1.0) )
    fleaf <- out_root$root
  }

  print(fleaf)
  
  # fleaf <- 0.5

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
  
  nturnover_ag <- nplant_ag / tau_plant
  nplant_ag    <- nplant_ag - nturnover_ag
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

