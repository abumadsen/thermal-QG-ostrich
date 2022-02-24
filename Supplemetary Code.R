#Title:
#Authors: Mads F. Schou
#Date: 17/06/2021

#------- 1. Thermal plasticity part I (Supplementary table 1)

#------- 2. Estimating non-linear selection gradient (Supplementary tables 2-3)

#------- 3. Tests of stabilizing selection - non-animal models (Supplementary tables 4-5)

#------- 4. Tests of stabilizing selection - animal models (Supplementary tables 6-7)

#------- 5. Thermal plasticity part II (Supplementary tables 8-9)

#------- 6. Phenotypic components of thermal plasticity (Supplementary tables 10-12)

#------- 7. Thermal plasticity of populations and their hybrids (Supplementary tables 13-14)

#------- 8. Additional tests of stabilizing selection - animal models (Supplementary tables 15-16)

#------- 9. Additional tests of stabilizing selection - non-animal models (Supplementary tables 17-18)

#------- 10. Additional estimation non-linear selection gradient (Supplementary tables 19-20)

#------- 11. Estimating evolutionary parameters

pacman::p_load(MCMCglmm,parallel,reshape)

Myburn = 100000
Mythin = 5000
Mynitt = 9000000 + Myburn
Nsamples = (Mynitt-Myburn)/Mythin

#For selection tests:
# Myburn = 100000
# Mythin = 3000
# Mynitt = 3000000 + Myburn
# Nsamples = (Mynitt-Myburn)/Mythin

#---- Data descriptions:
# cbind(Neggs, NoEggsHalf) = number of two-day intervals with and without an egg (within the specified temperature range, T.Con or T.State)
# T.Con = Continous change in temperature from thermal optimum (20C), used for random regression models
# T.Type = Type of temperatures change (cold or heat), used for random regression models
# Pop.fem = Population of female (SAB, ZB, KR or Hybrid) {note: not all populations are included in all models, see details in methods}
# Pop.male = Population of male (SAB, ZB, KR or Hybrid) {note: not all populations are included in all models, see details in methods}
# Age.fem_z = age of female (continous, scaled and centered)
# enclosure = enclosure id of the enclosure where female and male resided (factorial) 
# year = breeding year where observations were made (factorial) 
# animal = reserved term for the id of female as given in the pedigree (MyPed), genetic term
# damid = id of female (factorial), among individual term
# year_damid = female id split by breeding year, paste(year, damid, sep = "_") (factorial), within individual term
# units = reserved name for the residual term

# Extra variables in character state models:
# T.State = Thermal state (cold, benign and hot)

# Extra variables in selection models:
# cbind(NeggsTot, NoEggsTot) = number of two-day intervals with and without an egg (across the entire breeding year), in MS refered to as reproductive success
# relfit = Relative Fitness (see details in methods)
# HeatTol_z = Relative heat tolerance (continous, scaled and centered), see details in methods.
# HeatTol_z2 = HeatTol_z^2
# ColdTol_z = Relative cold tolerance (continous, scaled and centered), see details in methods.
# ColdTol_z2 = ColdTol_z^2
# trait = reserved name for the different response variables
# HeatTolFC_z = Relative heat toleranceFC (continous, scaled and centered), see details in methods.
# HeatTolFC_z2 = HeatTolFC_z^2
# ColdTolFC_z = Relative cold toleranceFC (continous, scaled and centered), see details in methods.
# ColdTolFC_z2 = ColdTolFC_z^2

# All models were run 3 times like this
# m1_3 <- mclapply(1:3, function(i) {
#   MCMCglmm()
# }, mc.cores=3)

#############################################################
# 1. Thermal plasticity part I (Supplementary table 1)
#############################################################

### Supplementary Table 1: Random slope animal model of thermal plasticity  

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002)))

m1 <- MCMCglmm(cbind(Neggs, NoEggsHalf) ~ T.Con:T.Type*Pop.fem + Pop.male + Age.fem_z,
         random = ~ us(1 + T.Con:T.Type):damid + us(1+T.Con:T.Type):animal + enclosure + year,
         data   = repro,
         family = "multinomial2",
         prior  = prior.m1,
         pedigree = MyPed,
         burnin = Myburn, thin = Mythin,,nitt = Mynitt,
         verbose = TRUE,
         pr = TRUE)



#############################################################
# 2. Estimating non-linear selection gradient (Supplementary tables 2-3)
#############################################################

### Supplementary Table 2: Non-linear selection gradient cold tolerance 

Myburn_selgrad = 100000
Mythin_selgrad = 3000
Mynitt_selgrad = 3000000 + Myburn
Nsamples_selgrad = (Mynitt-Myburn)/Mythin

prior.m2.3 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(1),
                 nu        = 0.002),
         G2=list(V        = diag(1),
                 nu        = 0.002)))

MCMCglmm(relfit ~ ColdTol_z + ColdTol_z2 + Pop.fem + Pop.male + Age.fem_z,
         random = ~ damid + enclosure,
         rcov = ~ units,
         data   = repro_cold,
         family = "gaussian",
         prior  = prior.m2.3,
         burnin = Myburn_selgrad, thin = Mythin_selgrad,nitt = Mynitt_selgrad,
         #burnin = 0, thin = 1,nitt = 100,
         verbose = TRUE,
         pr = TRUE)

### Supplementary Table 3: Non-linear selection gradient heat tolerance

MCMCglmm(relfit ~ HeatTol_z + HeatTol_z2 + Pop.fem + Pop.male + Age.fem_z,
         random = ~ damid + enclosure,
         rcov = ~ units,
         data   = repro_heat,
         family = "gaussian",
         prior  = prior.m2.3,
         burnin = Myburn_selgrad, thin = Mythin_selgrad,nitt = Mynitt_selgrad,
         #burnin = 0, thin = 1,nitt = 100,
         verbose = TRUE,
         pr = TRUE)


#############################################################
# 3. Tests of stabilizing selection - non-animal models (Supplementary tables 4-5)
#############################################################


### Supplementary Table 4: Test of stabilizing selection for cold resilience

prior.m4.5 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE),
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(3),
                 nu        = 2.002)))

m4 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot), ColdTol_z,ColdTol_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
               random = ~ us(trait):damid + idh(trait):enclosure + idh(trait):year,
               rcov = ~us(trait):units,
               data   = repro_cold,
               family = c("multinomial2",rep("gaussian",2)),
               prior  = prior.m4.5,
               burnin = Myburn, thin = Mythin,nitt = Mynitt,
               verbose = TRUE,
               pr = TRUE)

### Supplementary Table 5: Test of stabilizing selection for heat resilience

m5 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot), HeatTol_z,HeatTol_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
               random = ~ us(trait):damid + idh(trait):enclosure + idh(trait):year,
               rcov = ~us(trait):units,
               data   = repro_heat,
               family = c("multinomial2",rep("gaussian",2)),
               prior  = prior.m4.5,
               burnin = Myburn, thin = Mythin,nitt = Mynitt,
               verbose = TRUE,
               pr = TRUE)


#############################################################
# 4. Tests of stabilizing selection - animal models (Supplementary tables 6-7) 
#############################################################


### Supplementary Table 6: Test of stabilizing selection for cold resilience

prior.m6.7 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(3),
                 nu        = 2.002),
         G4=list(V        = diag(3),   
                 nu        = 2.002)))

m6 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot), ColdTol_z,ColdTol_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
               random = ~ us(trait):damid + us(trait):animal + idh(trait):enclosure + idh(trait):year,
               rcov = ~us(trait):units,
               data   = repro_cold,
               family = c("multinomial2",rep("gaussian",2)),
               prior  = prior.m6.7,
               pedigree = MyPed,
               burnin = Myburn, thin = Mythin,nitt = Mynitt,
               verbose = TRUE,
               pr = TRUE)


### Supplementary Table 7: Test of stabilizing selection for heat resilience

m7 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot), HeatTol_z,HeatTol_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
               random = ~ us(trait):damid + us(trait):animal + idh(trait):enclosure + idh(trait):year,
               rcov = ~us(trait):units,
               data   = repro_heat,
               family = c("multinomial2",rep("gaussian",2)),
               prior  = prior.m6.7,
               pedigree = MyPed,
               burnin = Myburn, thin = Mythin,nitt = Mynitt,
               verbose = TRUE,
               pr = TRUE)


#############################################################
# 5. Thermal plasticity part II (Supplementary tables 8-9)
#############################################################

### Supplementary Table 8: Random slope animal model with year-specific slopes of thermal plasticity

prior.m8 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(3),
                 nu        = 2.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002),
         G5=list(V        = diag(1),   
                 nu        = 0.002)))

m2 <- MCMCglmm(cbind(Neggs, NoEggsHalf) ~ T.Con:T.Type*Pop.fem + Pop.male + Age.fem_z,
         random = ~ us(1 + T.Con:T.Type):damid + us(1 + T.Con:T.Type):year_damid + us(1+T.Con:T.Type):animal + enclosure + year,
         data   = repro,
         family = "multinomial2",
         prior  = prior.m8,
         pedigree = MyPed,
         burnin = Myburn, thin = Mythin,,nitt = Mynitt,
         verbose = TRUE,
         pr = TRUE)

### Supplementary Table 9: Character-state animal model of thermal plasticity

prior.m9 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 n        = 2.002),
         G2=list(V        = diag(3),
                 n        = 2.002),
         G3=list(V        = diag(1),   
                 n        = 0.002),
         G4=list(V        = diag(1),   
                 n        = 0.002)))

m3 <- MCMCglmm(cbind(Neggs, NoEggsHalf) ~ 1 + T.State*Pop.fem + Pop.male + Age.fem_z,
         random = ~ us(T.State):damid + us(T.State):animal + enclosure + year,
         rcov = ~idh(T.State):units,
         data   = repro,
         family = "multinomial2",
         prior  = prior.m9,
         pedigree = MyPed,
         burnin = Myburn, thin = Mythin,,nitt = Mynitt,
         verbose = TRUE,
         pr = FALSE)


#############################################################
# 6. Phenotypic components of thermal plasticity (Supplementary tables 10-12)
#############################################################


### Supplementary Table 10: Random slope model of thermal plasticity

prior.m8 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(1),   
                 nu        = 0.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002)))

m8 <-  MCMCglmm(cbind(Neggs, NoEggsHalf) ~ T.Con:T.Type*Pop.fem + Pop.male + Age.fem_z,
                random = ~ us(1 + T.Con:T.Type):damid + enclosure + year,
                data   = repro,
                family = "multinomial2",
                prior  = prior.m8,
                burnin = Myburn, thin = Mythin,,nitt = Mynitt,
                verbose = TRUE,
                pr = TRUE)


### Supplementary Table 11: Random slope model with year-specific slopes of thermal plasticity

prior.m9 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002)))

m9 <- MCMCglmm(cbind(Neggs, NoEggsHalf) ~ T.Con:T.Type*Pop.fem + Pop.male + Age.fem_z,
               random = ~ us(1 + T.Con:T.Type):damid + us(1 + T.Con:T.Type):year_damid + enclosure + year,
               data   = repro,
               family = "multinomial2",
               prior  = prior.m9,
               burnin = Myburn, thin = Mythin,,nitt = Mynitt,
               verbose = TRUE,
               pr = TRUE)

### Supplementary Table 12: Character-state model of thermal plasticity

prior.m10 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 3.002 because 3 in diag!
  G=list(G1=list(V        = diag(3),
                 n        = 2.002),
         G2=list(V        = diag(1),   
                 n        = 0.002),
         G3=list(V        = diag(1),   
                 n        = 0.002)))

m10 <- MCMCglmm(cbind(Neggs, NoEggsHalf) ~ 1 + T.State*Pop.fem + Pop.male + Age.fem_z,
                random = ~ us(T.State):damid + enclosure + year,
                rcov = ~idh(T.State):units,
                data   = repro,
                family = "multinomial2",
                prior  = prior.m10,
                burnin = Myburn, thin = Mythin,,nitt = Mynitt,
                verbose = TRUE,
                pr = FALSE)


#############################################################
# 7. Thermal plasticity of populations and their hybrids (Supplementary tables 13-14)
#############################################################

### Supplementary Table 13: Population specific random slope model of thermal plasticity

prior.m13 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(3),
                 nu        = 2.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002),
         G5=list(V        = diag(1),   
                 nu        = 0.002)))


m11 <- MCMCglmm(cbind(Neggs, NoEggsHalf) ~ Pop.fem-1 + T.Con:T.Type:Pop.fem + Pop.male + Age.fem_z,
                random = ~
                  us(1 + T.Con:T.Type:at.level(Pop.fem,"ZB")):damid +
                  us(1 + T.Con:T.Type:at.level(Pop.fem,"KR")):damid +
                  us(1 + T.Con:T.Type:at.level(Pop.fem,"SAB")):damid +
                  enclosure + year,
                data   = repro,
                family = "multinomial2",
                prior  = prior.m13,
                burnin = Myburn, thin = Mythin,nitt = Mynitt,
                verbose = TRUE,
                pr = TRUE)


### Supplementary Table 14: Population hybrids specific random slope model of thermal plasticity

prior.m14 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(3),
                 nu        = 2.002),
         G4=list(V        = diag(3),   
                 nu        = 2.002),
         G5=list(V        = diag(3),   
                 nu        = 2.002),
         G6=list(V        = diag(1),   
                 nu        = 0.002),
         G7=list(V        = diag(1),   
                 nu        = 0.002)))

m12 <- MCMCglmm(cbind(Neggs, NoEggsHalf) ~ Pop.fem-1 + T.Con:T.Type:Pop.fem + Pop.male + Age.fem_z,
         random = ~
           us(1 + T.Con:T.Type:at.level(Pop.fem,"ZB")):damid +
           us(1 + T.Con:T.Type:at.level(Pop.fem,"KR")):damid +
           us(1 + T.Con:T.Type:at.level(Pop.fem,"SAB")):damid +
           us(1 + T.Con:T.Type:at.level(Pop.fem,"SABxKR")):damid +
           us(1 + T.Con:T.Type:at.level(Pop.fem,"SABxZB")):damid +
           enclosure + year,
         data   = repro,
         family = "multinomial2",
         prior  = prior.m14,
         burnin = Myburn, thin = Mythin,nitt = Mynitt,
         #burnin = 0, thin = 1,nitt = 100,
         verbose = TRUE,
         pr = TRUE)


#############################################################
# 8. Additional tests of stabilizing selection - animal models (Supplementary tables 15-16)
#############################################################


### Supplementary Table 15: Additional test of stabilizing selection for cold tolerance

prior.m15.16 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(3),
                 nu        = 2.002),
         G4=list(V        = diag(3),   
                 nu        = 2.002)))

m22 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot),ColdTolFC_z,ColdTolFC_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
                random = ~ us(trait):damid + us(trait):animal + idh(trait):enclosure + idh(trait):year,
                rcov = ~us(trait):units,
                data   = repro_cold,
                family = c("multinomial2",rep("gaussian",2)),
                prior  = prior.m15.16,
                pedigree = MyPed,
                burnin = Myburn, thin = Mythin,nitt = Mynitt,
                verbose = TRUE,
                pr = TRUE)

### Supplementary Table 16: Additional test of stabilizing selection for heat tolerance

m23 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot), HeatTolFC_z,HeatTolFC_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
                random = ~ us(trait):damid + us(trait):animal + idh(trait):enclosure + idh(trait):year,
                rcov = ~us(trait):units,
                data   = repro_heat,
                family = c("multinomial2",rep("gaussian",2)),
                prior  = prior.m15.16,
                pedigree = MyPed,
                burnin = Myburn, thin = Mythin,nitt = Mynitt,
                verbose = TRUE,
                pr = TRUE)


#############################################################
# 9. Additional tests of stabilizing selection - non-animal models (Supplementary tables 17-18)
#############################################################


### Supplementary Table 17: Additional test of stabilizing selection for cold tolerance

prior.m17.18 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(3),
                 nu        = 2.002)))

m24 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot), ColdTolFC_z,ColdTolFC_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
                random = ~ us(trait):damid + idh(trait):enclosure + idh(trait):year,
                rcov = ~us(trait):units,
                data   = repro_cold,
                family = c("multinomial2",rep("gaussian",2)),
                prior  = prior.m17.18,
                burnin = Myburn, thin = Mythin,nitt = Mynitt,
                verbose = TRUE,
                pr = TRUE)

### Supplementary Table 18: Additional test of stabilizing selection for heat tolerance


m25 <- MCMCglmm(cbind(cbind(NeggsTot, NoEggsTot), HeatTolFC_z,HeatTolFC_z2) ~ trait-1 + Pop.fem*trait + Pop.male + Age.fem_z,
                random = ~ us(trait):damid + idh(trait):enclosure + idh(trait):year,
                rcov = ~us(trait):units,
                data   = repro_heat,
                family = c("multinomial2",rep("gaussian",2)),
                prior  = prior.m17.18,
                burnin = Myburn, thin = Mythin,nitt = Mynitt,
                verbose = TRUE,
                pr = TRUE)

#############################################################
# 10. Additional estimation non-linear selection gradient (Supplementary tables 19-20)
#############################################################

### Supplementary Table 19: Non-linear selection gradient cold tolerance

prior.m19.20 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  #Not fixed at one as binomial model!, inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(1),
                 nu        = 0.002),
         G2=list(V        = diag(1),
                 nu        = 0.002)))

MCMCglmm(relfit ~ ColdTolFC_z + ColdTolFC_z2 + Pop.fem + Pop.male + Age.fem_z,
         random = ~ damid + enclosure,
         rcov = ~ units,
         data   = repro_cold,
         family = "gaussian",
         prior  = prior.m19.20,
         burnin = Myburn_selgrad, thin = Mythin_selgrad,nitt = Mynitt_selgrad,
         #burnin = 0, thin = 1,nitt = 100,
         verbose = TRUE,
         pr = TRUE)

### Supplementary Table 20: Non-linear selection gradient heat tolerance

MCMCglmm(relfit ~ HeatTolFC_z + HeatTolFC_z2 + Pop.fem + Pop.male + Age.fem_z,
         random = ~ damid + enclosure,
         rcov = ~ units,
         data   = repro_heat,
         family = "gaussian",
         prior  = prior.m19.20,
         burnin = Myburn_selgrad, thin = Mythin_selgrad,nitt = Mynitt_selgrad,
         #burnin = 0, thin = 1,nitt = 100,
         verbose = TRUE,
         pr = TRUE)


#############################################################
# 11. Estimating evolutionary parameters
#############################################################


#------------ Quadratic selection gradients -----------------#


Nonlinear <- cbind(posterior.mode(m1$Sol[,"ColdTol_z2"]*2), HPDinterval(m1$Sol[,"ColdTol_z2"]*2))


#------------ Heritability of slopes, i.e. heat resilience (random regression models) -----------------#


h2_slope_heat<-m2$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.animal']/
  (
    m2$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.animal']+
      m2$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.damid']+
      m2$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.year_damid'])

h2_slope_heatPM <- data.frame(t(c(posterior.mode(h2_slope_heat),HPDinterval(h2_slope_heat))))


#------------ Change in Va, Vp and H2 with increasing thermal T.Type -----------------#


covariate_values_hot <- seq(min(repro$T.Con[repro$T.Type %in% "Heat"]),max(repro$T.Con[repro$T.Type %in% "Heat"]),0.01)
covariate_values_cold <- seq(min(repro$T.Con[repro$T.Type %in% "Cold"]),max(repro$T.Con[repro$T.Type %in% "Cold"]),0.01)

# #----- Va function
get.va.across.slope.yrid = function(myX = NULL, Va.int = NULL, Va.slope = NULL, Va.cov = NULL,Ves.int = NULL, Ves.slope = NULL, Ves.cov = NULL,Vperm.int = NULL, Vperm.slope = NULL, Vperm.cov = NULL, enclosure = NULL, year = NULL, units = NULL) {
  Va.all =    Va.int + 2*myX*Va.cov + Va.slope*(myX^2)
  return(Va.all)
}

# #----- h2 function
get.h2.across.slope.yrid = function(myX = NULL, Va.int = NULL, Va.slope = NULL, Va.cov = NULL,Ves.int = NULL, Ves.slope = NULL, Ves.cov = NULL,Vperm.int = NULL, Vperm.slope = NULL, Vperm.cov = NULL, enclosure = NULL, year = NULL, units = NULL) {
  Ves.all =    Ves.int + 2*myX*Ves.cov + Ves.slope*(myX^2)
  Va.all =    Va.int + 2*myX*Va.cov + Va.slope*(myX^2)
  Vperm.all = Vperm.int + 2*myX*Vperm.cov + Vperm.slope*(myX^2)
  h2 = Va.all / (Va.all + Ves.all+ Vperm.all + enclosure + year + units)
  return(h2)
}

#----- Vp function
get.Vp.across.slope.yrid = function(myX = NULL, Va.int = NULL, Va.slope = NULL, Va.cov = NULL,Ves.int = NULL, Ves.slope = NULL, Ves.cov = NULL,Vperm.int = NULL, Vperm.slope = NULL, Vperm.cov = NULL, enclosure = NULL, year = NULL, units = NULL) {
  Ves.all =    Ves.int + 2*myX*Ves.cov + Ves.slope*(myX^2)
  Va.all =    Va.int + 2*myX*Va.cov + Va.slope*(myX^2)
  Vperm.all = Vperm.int + 2*myX*Vperm.cov + Vperm.slope*(myX^2)
  Vp = (Ves.all + Va.all + Vperm.all + enclosure + year + units)
  return(Vp)
}

#Each of these functions are then applied like this (here get.va.across.slope.yrid at heat T.Type)
Vaout <- as.numeric()
for(i in covariate_values_hot){
  REPpost <- get.va.across.slope.yrid(myX = i, 
                                      Va.int = m2$VCV[,'(Intercept):(Intercept).animal'], Va.slope = m2$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.animal'], Va.cov = m2$VCV[,'(Intercept):T.Con:T.TypeHeat.animal'],
                                      Ves.int = m2$VCV[,'(Intercept):(Intercept).year_damid'], Ves.slope = m2$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.year_damid'], Ves.cov = m2$VCV[,'(Intercept):T.Con:T.TypeHeat.year_damid'],
                                      Vperm.int = m2$VCV[,'(Intercept):(Intercept).damid'], Vperm.slope = m2$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.damid'], Vperm.cov = m2$VCV[,'(Intercept):T.Con:T.TypeHeat.damid'], 
                                      enclosure =  m2$VCV[,'enclosure'], year = m2$VCV[,'year'], units = m2$VCV[,'units'])
  
  REP <- data.frame(t(c(posterior.mode(REPpost),HPDinterval(REPpost))))
  
  Vaout <- rbind(Vaout,REP)
}


#------------ Correlations -------------------------#

#Correlations were estimated with this function which spits out all random effect correlations within us matrices in the model (adjusted from Charlie K. Cornwallis)

ReportCorrelationsMCMC = function(x, roundto = 2){
  #Get vars
  Vars <- ReportRandomVarianceMCMC(x)[c("Random Effects: Variances","Level")]
  colnames(Vars)[1] <- "Var"
  #"units" was automatically renamed to "residuals"; we reverse it
  Vars$Level[Vars$Level %in% "residuals"] <- "units"
  #Vars <- Vars[Vars$Level != "units",]
  
  #Get covars within each level
  Covars <- as.numeric()
  for(mylevel in levels(factor(Vars$Level))){
    tempCovars <- expand.grid(Vars[Vars$Level %in% mylevel ,c(1,1)])
    tempCovars$Level <- mylevel
    Covars <- rbind(Covars,tempCovars)
  }
  Covars <- Covars[Covars$Var != Covars$Var.1,]
  Covars <- Covars[!duplicated(Covars),]
  Covars$CovarNames <- paste(paste(Covars$Var,Covars$Var.1,sep = ":"),Covars$Level, sep = ".") #Gives two of each covar, just as in summary$VCV
  Covars$VarNames1 <- paste(paste(Covars$Var,Covars$Var,sep = ":"),Covars$Level, sep = ".")
  Covars$VarNames2 <- paste(paste(Covars$Var.1,Covars$Var.1,sep = ":"),Covars$Level, sep = ".")
  #Estimate correlations
  corrs <- NULL
  corrs=matrix(nrow = nrow(x$VCV), ncol = 0)
  for(i in 1:nrow(Covars)){
    covar.post <- x$VCV[,colnames(x$VCV) %in% Covars$CovarNames[i]]
    var1.post <- x$VCV[,colnames(x$VCV) %in% Covars$VarNames1[i]]
    var2.post <- x$VCV[,colnames(x$VCV) %in% Covars$VarNames2[i]]
    tmpcor<-covar.post/sqrt(var1.post*var2.post)
    
    if(!all(is.na(tmpcor[]))){ #If not empty: covar estimated by model
      corrs<-cbind(corrs,tmpcor)
      colnames(corrs)[ncol(corrs)] <- Covars$CovarNames[i]
    }
  }
  corrs<-as.mcmc(corrs)
  
  #Cor Summaries
  cor1=paste(round(posterior.mode(corrs),roundto)," (",round(HPDinterval(corrs)[,1],roundto), ",",round(HPDinterval(corrs)[,2],roundto),")",sep="")
  ncors<-ifelse(is.null(dim(corrs)), 1,dim(corrs)[2])
  nits<-ifelse(is.null(dim(corrs)),length(corrs), dim(corrs)[1])
  if(ncors >1){
    pCor=pmax(0.5/nits, pmin(colSums(corrs[,1:ncors, drop = FALSE] > 0)/nits, 1 - colSums(corrs[, 1:ncors, drop = FALSE] > 0)/nits))*2
  } else  {
    pCor=pmax(0.5/nits, pmin(sum(corrs[,drop = FALSE] > 0)/nits, 1 - sum(corrs[, drop = FALSE] > 0)/nits))*2
  }
  randomCorr<-data.frame("Random Effects: Correlations"=colnames(corrs),"Posterior Mode (CI)"=cor1,"pMCMC"=round(pCor,3), check.names=FALSE)
  #Remove duplicates (each corr will apear twice as they do so in the VCV)
  randomCorr <- randomCorr[!duplicated(randomCorr$`Posterior Mode (CI)`),]
  randomCorr[,c("Random Effects: Correlations")] <-  gsub("units","residuals",randomCorr[,c("Random Effects: Correlations")])
  return(randomCorr)
}

#The key element of the function is "covar.post/sqrt(var1.post*var2.post)"
#This function was therefore used to get correlations between reproductive success and heat (or cold) resilience^2 in the selection test models.
#It was also used to estimate correlations between cold and heat resilience. 
#Below a few examples of how this can also be done without the complex function above.

#------------ Genetic correlations -----------------#


# Benign vs Resilience to heat (intercept vs heat slope)

rg.bh <- m1$VCV[,'(Intercept):T.Con:T.TypeHeat.animal']/
  sqrt(m1$VCV[,'(Intercept):(Intercept).animal']*m1$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.animal'])
rg.bhPM <- data.frame(t(c(posterior.mode(rg.bh),HPDinterval(rg.bh))))

# Resilience to cold vs Resilience to heat (cold slope vs heat slope)

rgch <-m1$VCV[,'T.Con:T.TypeCold:T.Con:T.TypeHeat.animal']/
  sqrt(m1$VCV[,'T.Con:T.TypeCold:T.Con:T.TypeCold.animal']*m1$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.animal'])
rgchPM <- data.frame(t(c(posterior.mode(rgch),HPDinterval(rgch))))


#------------ Phenotypic correlations -----------------#

# Benign vs Resilience to heat (intercept vs heat slope)

rp.bh <- m9$VCV[,'(Intercept):T.Con:T.TypeHeat.damid']/
  sqrt(m9$VCV[,'(Intercept):(Intercept).damid']*m9$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.damid'])
rp.bhPM <- data.frame(t(c(posterior.mode(rp.bh),HPDinterval(rp.bh))))

# Resilience to cold vs Resilience to heat (cold slope vs heat slope)

rpch <-m9$VCV[,'T.Con:T.TypeCold:T.Con:T.TypeHeat.damid']/
  sqrt(m9$VCV[,'T.Con:T.TypeCold:T.Con:T.TypeCold.damid']*m9$VCV[,'T.Con:T.TypeHeat:T.Con:T.TypeHeat.damid'])
rpchPM <- data.frame(t(c(posterior.mode(rpch),HPDinterval(rpch))))



#---------- differences in hot vs cold resilience across Populations and crosses -----------------#

#Example with KR (Kenyan reds)
RedHotCold <- m12$Sol[,grep("Pop.femKR:T.Con:T.TypeHeat", colnames(m12$Sol))]-m12$Sol[,grep("Pop.femKR:T.Con:T.TypeCold", colnames(m12$Sol))]
RedHotCold <- as.mcmc(RedHotCold)
RedHotCold <- c(posterior.mode(RedHotCold),HPDinterval(RedHotCold))
