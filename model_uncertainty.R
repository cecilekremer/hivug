args <- commandArgs(trailingOnly =  TRUE)
f_main <- as.numeric(args[1])
f_reg <- as.numeric(args[2])
actsf11_main <- as.numeric(args[3])
actsf12_main <- as.numeric(args[4])
actsf13_main <- as.numeric(args[5])
actsf21_main <- as.numeric(args[6])
actsf22_main <- as.numeric(args[7])
actsf23_main <- as.numeric(args[8])
actsf31_main <- as.numeric(args[9])
actsf32_main <- as.numeric(args[10])
actsf33_main <- as.numeric(args[11])
actsf22_reg <- as.numeric(args[12])
actsf23_reg <- as.numeric(args[13])
actsf32_reg <- as.numeric(args[14])
actsf33_reg <- as.numeric(args[15])
actsf33_casual <- as.numeric(args[16])
df3_casual <- as.numeric(args[17])
dm3_casual <- as.numeric(args[18])
df2_reg <- as.numeric(args[19])
dm2_reg <- as.numeric(args[20])
df3_reg <- as.numeric(args[21])
dm3_reg <- as.numeric(args[22])
df1_main <- as.numeric(args[23])
dm1_main <- as.numeric(args[24])
df2_main <- as.numeric(args[25])
dm2_main <- as.numeric(args[26])
df3_main <- as.numeric(args[27])
dm3_main <- as.numeric(args[28])
A_main <- as.numeric(args[29])
A_reg <- as.numeric(args[30])
beta_f <- as.numeric(args[31])
beta_m <- as.numeric(args[32])
prop_circ <- as.numeric(args[33])
sim <- as.numeric(args[34])
pCIRC <- as.numeric(args[35])
pHCT <- as.numeric(args[36])
tau1 <- as.numeric(args[37])
tau2 <- as.numeric(args[38])
tau3 <- as.numeric(args[39])
tau4 <- as.numeric(args[40])
HCT_m <- as.numeric(args[41])
ART_m <- as.numeric(args[42])
const1 <- as.numeric(args[43])
const2 <- as.numeric(args[44])
const3 <- as.numeric(args[45])
f_casual <- as.numeric(args[46])
scaleART1 <- as.numeric(args[47])
scaleART2 <- as.numeric(args[48])
scaleART3 <- as.numeric(args[49])
scaleART4 <- as.numeric(args[50])
scaleCIRC <- as.numeric(args[51])
scaleCOND1 <- as.numeric(args[52])
scaleCOND2 <- as.numeric(args[53])
scaleCOND3 <- as.numeric(args[54])
pHCT2 <- as.numeric(args[55])
pPREP <- as.numeric(args[56])
PhctS <- as.numeric(args[57])

######################################################################################################
### BASELINE
#################
### FUNCTIONS ###

## HCT uptake (monthly rate) / pHCT = desired proportion aware
fun.tauT = function(t,rHCT,pHCT,prop, rHCT2, pHCT2){
  tau_t=0
  if(t>217 & t<=(217+1/rHCT)){ # testing since 1987
    tau_t = ((t-217)*rHCT*pHCT) - prop
  } 

  else if(t>(217+1/rHCT)){ # increase prop aware to 65%
    tau_t = tau_t + ((t-540)*rHCT2*(0.65-pHCT)) - (prop-pHCT)
  }
  
  return(tau_t)
}

## HCT for susceptibles (after 2014)
fun.hct = function(t,pmax){
  tau_t = 0
  return(tau_t)
}


## ART uptake (increasing probability of accepting treatment after positive test)
fun.tauA34 = function(t,pART,scaleART){ # pART = max prob of accepting treatment
  tau_t = 0
  if(t>420 & t<=540){
    h = 2
    tau_t = pART*((t-420)^h/((t-420)^h+10^h))
  }
  return(tau_t)
}

fun.tauA12 = function(t,pART,scaleART){ # pART = max prob of accepting treatment
  tau_t = 0
  if(t>504 & t<=540){
    h = 2
    tau_t = pART*((t-504)^h/((t-504)^h+10^h))
  }
  return(tau_t)
}

## PrEP from 2015 onward
fun.prep = function(t, theta){
  tau_p = 0
  return(tau_p)
}

## Circumcision (no 'extra' circumcision before HCT in 1987)
fun.circ = function(t,pCIRC,scaleCIRC){
  tau_t = 0
  if(t>217 & t<=540){
    h = 4
    tau_t = pCIRC*((t-217)^h/((t-217)^h+100^h))
  }
  return(tau_t)
}

## Condom use
fun.cond = function(t,cond_max,scaleCOND){
  h=4
  cond = cond_max*((t-0)^h/((t-0)^h+200^h))
  return(cond)
}


#####################################
library(deSolve)

### INITIAL STATES

init=0.01

N = 2000

Nf1=0.43*N; Nm1=0.33*N
Nf2=0.02*N; Nm2=0.07*N
Nf3=0.05*N ; Nm3=0.1*N 

H2m1=0; H2m2=0; H2m3=0; H2f1=0; H2f2=0; H2f3=0
H3m1=0; H3m2=0; H3m3=0; H3f1=0; H3f2=0; H3f3=0

Y2m1=0; Y2m2=0; Y2m3=0; Y2f1=0; Y2f2=0; Y2f3=0
Y3m1=0; Y3m2=0; Y3m3=0; Y3f1=0; Y3f2=0; Y3f3=0

Z2m1=0; Z2m2=0; Z2m3=0; Z2f1=0; Z2f2=0; Z2f3=0
Z3m1=0; Z3m2=0; Z3m3=0; Z3f1=0; Z3f2=0; Z3f3=0

A2m1=0; A2m2=0; A2m3=0; A2f1=0; A2f2=0; A2f3=0
A3m1=0; A3m2=0; A3m3=0; A3f1=0; A3f2=0; A3f3=0
# 
H1m1=init*Nm1; Y1m1=0; Z1m1=0; A1m1=0
H1f1=init*Nf1; Y1f1=0; Z1f1=0; A1f1=0
H1m2=init*Nm2; Y1m2=0; Z1m2=0; A1m2=0
H1f2=init*Nf2; Y1f2=0; Z1f2=0; A1f2=0
H1m3=init*Nm3; Y1m3=0; Z1m3=0; A1m3=0
H1f3=init*Nf3; Y1f3=0; Z1f3=0; A1f3=0

S2f1=0; S2f2=0; S2f3=0; S2m1=0; S2m2=0; S2m3=0
S3m1=0; S3m2=0; S3m3=0; S4m1=0; S4m2=0; S4m3=0

S1f1=Nf1-S2f1-H1f1-Y1f1-Z1f1-A1f1-H2f1-Y2f1-Z2f1-A2f1-H3f1-Y3f1-Z3f1-A3f1
S1f2=Nf2-S2f2-H1f2-Y1f2-Z1f2-A1f2-H2f2-Y2f2-Z2f2-A2f2-H3f2-Y3f2-Z3f2-A3f2 
S1f3=Nf3-S2f3-H1f3-Y1f3-Z1f3-A1f3-H2f3-Y2f3-Z2f3-A2f3-H3f3-Y3f3-Z3f3-A3f3 

S1m1=Nm1-S2m1-S3m1-S4m1-H1m1-Y1m1-Z1m1-A1m1-H2m1-Y2m1-Z2m1-A2m1-H3m1-Y3m1-Z3m1-A3m1
S1m2=Nm2-S2m2-S3m2-S4m2-H1m2-Y1m2-Z1m2-A1m2-H2m2-Y2m2-Z2m2-A2m2-H3m2-Y3m2-Z3m2-A3m2 
S1m3=Nm3-S2m3-S3m3-S4m3-H1m3-Y1m3-Z1m3-A1m3-H2m3-Y2m3-Z2m3-A2m3-H3m3-Y3m3-Z3m3-A3m3 

########################
### INPUT PARAMETERS ###
########################

## Demographic
psi = 0.0019
mu = 0.0019 + (1.025^(1/12)-1)

## Epidemiologic
avg.nu1 = 3.25; avg.omega1 = 100; avg.epi1 = 14; avg.delta1 = 16
avg.nu2 = 98; avg.omega2 = 648; avg.epi2 = 30; avg.delta2 = 33


## Intervention parameters
ART_eff = 0.96 # reduced infectiousness
cond_eff = 0.1
HCT_cond = 1.1 # 10% increase in condom use (1=no increase) due to HCT
e_circ = 0.6 # relative susceptibility of circumcised men = 0.4 = reduction of 60%
prop_circ = prop_circ # proportion circumcised when entering the population
HCT_red = 0.9 # 10% reduction in number of sex acts (1=no reduction) due to HCT
AIDS_red = 0.1 # 90% reduction in number of sex acts due to untreated AIDS

## Intervention uptake
rHCT = 1/323 # HCT coverage increases over 323 months 
pHCT = pHCT # specified coverage level 2014 (proportion aware of HIV+ status)
tau1 = (1+(-log(1-tau1)))^(1/12) - 1 # proportion initiating ART on average one year after testing
tau2 = (1+(-log(1-tau2)))^(1/12) - 1
tau3 = (1+(-log(1-tau3)))^(1/12) - 1
tau4 = (1+(-log(1-tau4)))^(1/12) - 1
pCIRC = (1+(-log(1-pCIRC)))^(1/12) - 1 # pCIRC = proportion accepting VMMC on average one year after testing

gamma = (1+(-log(1-0.2)))^(1/12) - 1 # monthly rate of ART dropout (20% of individuals no longer on treatment one year after starting)

ART_m = ART_m # ART uptake relative to females
HCT_m = HCT_m # HCT uptake relative to females

## PrEP
# pPREP = 0.5 # prop initiating prep on average one year after testing
# rPREP = 1/168
pPREP = (1+(-log(1-pPREP)))^(1/12) - 1 # monthly rate of PrEP uptake; theta = proportion initiating PrEP on average one year after testing
om = (1+(-log(1-0.2)))^(1/12) - 1 # monthly rate of PrEP discontinuation (20% of individuals no longer on PrEP one year after starting)
# om = 0
e_p = 0.67 # relative reduction in susceptibility on PrEP
prepF = 0.6 # prep adherence
prepM = 0.6
# prep efficacy = e_p*prepF
e_pF = e_p*prepF
e_pM = e_p*prepM

##### Intervention scale-up
scaleART1 = (1+(-log(1-scaleART1)))^(1/12) - 1
scaleART2 = (1+(-log(1-scaleART2)))^(1/12) - 1
scaleART3 = (1+(-log(1-scaleART3)))^(1/12) - 1
scaleART4 = (1+(-log(1-scaleART4)))^(1/12) - 1
scaleCIRC = (1+(-log(1-scaleCIRC)))^(1/12) - 1
scaleCOND1 = scaleCOND1
scaleCOND2 = scaleCOND2
scaleCOND3 = scaleCOND3
rHCT2 = 1/180
pHCT2 = pHCT2
PhctS = (1+(-log(1-PhctS)))^(1/12) - 1 # yearly rate of accepting HCT for susceptibles

## Transmission
high_v1 = 26
high_v2 = 7
high_v3 = 1

## Behavioral
f_main = f_main # condom use
f_reg = f_reg
f_casual = f_casual

## Longer partnership = more sex acts
actsf11_main = actsf11_main
actsm11_main = actsf11_main
actsf12_main = actsf12_main
actsm21_main = actsf12_main
actsf13_main = actsf13_main
actsm31_main = actsf13_main
actsf21_main = actsf21_main
actsm12_main = actsf21_main
actsf22_main = actsf22_main
actsm22_main = actsf22_main
actsf23_main = actsf23_main
actsm32_main = actsf23_main
actsf31_main = actsf31_main
actsm13_main = actsf31_main
actsf32_main = actsf32_main
actsm23_main = actsf32_main
actsf33_main = actsf33_main
actsm33_main = actsf33_main
actsf22_reg = actsf22_reg
actsm22_reg = actsf22_reg
actsf23_reg = actsf23_reg
actsm32_reg = actsf23_reg
actsf32_reg = actsf32_reg
actsm23_reg = actsf32_reg
actsf33_reg = actsf33_reg
actsm33_reg = actsf33_reg
actsf33_casual = actsf33_casual
actsm33_casual = actsf33_casual

## Rate of partner acquisition (= average number of partners per month)
df3_casual = df3_casual
dm3_casual = dm3_casual
df2_reg = df2_reg
dm2_reg = dm2_reg
df3_reg = df3_reg
dm3_reg = dm3_reg
df1_main = df1_main
dm1_main = dm1_main
df2_main = df2_main
dm2_main = dm2_main
df3_main = df3_main
dm3_main = dm3_main

A_main = A_main # degree of assortativity for main partnerships
A_reg = A_reg


##################################
#### Transmission probabilities D

### Transmission parameters (depending on disease stage, PrEP use, ART, and circumcision status)
beta_f = beta_f # baseline (partner in chronic LV stage) assumed same for male/female infected partner
beta_m = beta_m
# Disease stages (1=initial, 2=preA, 3=AIDS)
beta1_f = high_v1*beta_f
beta1_m = high_v1*beta_m
beta2_f = high_v2*beta_f
beta2_m = high_v2*beta_m
beta3_f = 0 # no transmission in untreated AIDS stage
beta3_m = 0
# Infected partner on ART
beta4_f = beta_f*(1-ART_eff)
beta5_f = beta_f*(1-ART_eff)
beta6_f = beta_f*(1-ART_eff)
beta7_f = beta_f*(1-ART_eff)
beta4_m = beta_m*(1-ART_eff)
beta5_m = beta_m*(1-ART_eff)
beta6_m = beta_m*(1-ART_eff)
beta7_m = beta_m*(1-ART_eff)
# Susceptible on PrEP
# beta8_f = beta_f*(1-e_p)
# beta9_f = beta1_f*(1-e_p)
# beta10_f = beta2_f*(1-e_p)
# beta11_f = beta3_f*(1-e_p)
# beta8_m = beta_m*(1-e_p)
# beta9_m = beta1_m*(1-e_p)
# beta10_m = beta2_m*(1-e_p)
# beta11_m = beta3_m*(1-e_p)
beta8_f = beta_f*(1-e_pF)
beta9_f = beta1_f*(1-e_pF)
beta10_f = beta2_f*(1-e_pF)
beta11_f = beta3_f*(1-e_pF)
beta8_m = beta_m*(1-e_pM)
beta9_m = beta1_m*(1-e_pM)
beta10_m = beta2_m*(1-e_pM)
beta11_m = beta3_m*(1-e_pM)
# Circumcised susceptible males
beta12_m = beta_m*(1-e_circ)
beta13_m = beta1_m*(1-e_circ)
beta14_m = beta2_m*(1-e_circ)
beta15_m = beta3_m*(1-e_circ)
# Circumcised males on PrEP
# beta28_m = beta12_m*(1-e_p)
# beta29_m = beta13_m*(1-e_p)
# beta30_m = beta14_m*(1-e_p)
# beta31_m = beta15_m*(1-e_p)
beta28_m = beta12_m*(1-e_pM)
beta29_m = beta13_m*(1-e_pM)
beta30_m = beta14_m*(1-e_pM)
beta31_m = beta15_m*(1-e_pM)
# Infected partner on ART and susceptible on PrEP
# beta16_f = beta4_f*(1-e_p)
# beta17_f = beta5_f*(1-e_p)
# beta18_f = beta6_f*(1-e_p)
# beta19_f = beta7_f*(1-e_p)
# beta16_m = beta4_m*(1-e_p)
# beta17_m = beta5_m*(1-e_p)
# beta18_m = beta6_m*(1-e_p)
# beta19_m = beta7_m*(1-e_p)
beta16_f = beta4_f*(1-e_pF)
beta17_f = beta5_f*(1-e_pF)
beta18_f = beta6_f*(1-e_pF)
beta19_f = beta7_f*(1-e_pF)
beta16_m = beta4_m*(1-e_pM)
beta17_m = beta5_m*(1-e_pM)
beta18_m = beta6_m*(1-e_pM)
beta19_m = beta7_m*(1-e_pM)
# Infected partner on ART and susceptible male circumcised
beta20_m = beta12_m*(1-ART_eff)
beta21_m = beta12_m*(1-ART_eff)
beta22_m = beta12_m*(1-ART_eff)
beta23_m = beta12_m*(1-ART_eff)
# Infected partner on ART & susceptible male circumcised + PrEP
# beta24_m = beta20_m*(1-e_p)
# beta25_m = beta21_m*(1-e_p)
# beta26_m = beta22_m*(1-e_p)
# beta27_m = beta23_m*(1-e_p)
beta24_m = beta20_m*(1-e_pM)
beta25_m = beta21_m*(1-e_pM)
beta26_m = beta22_m*(1-e_pM)
beta27_m = beta23_m*(1-e_pM)


#######################
### ODE SYSTEM      ###
#######################

SHYZA_FM=function(t, state, parameters)
{
  
  with(as.list(c(state,parameters)),
       {
####################################################
### FORCE OF INFECTION --> FREQUENCY-DEPENDENT ! ###
####################################################
         
Nf1 = S1f1 + S2f1 + H1f1 + H2f1 + Y1f1 + Y2f1 + Z1f1 + Z2f1 + A1f1 + A2f1 + H3f1 + Y3f1 + Z3f1 + A3f1
Nf2 = S1f2 + S2f2 + H1f2 + H2f2 + Y1f2 + Y2f2 + Z1f2 + Z2f2 + A1f2 + A2f2 + H3f2 + Y3f2 + Z3f2 + A3f2
Nf3 = S1f3 + S2f3 + H1f3 + H2f3 + Y1f3 + Y2f3 + Z1f3 + Z2f3 + A1f3 + A2f3 + H3f3 + Y3f3 + Z3f3 + A3f3
Nm1 = S1m1 + S2m1 + S3m1 + S4m1 + H1m1 + H2m1 + Y1m1 + Y2m1 + Z1m1 + Z2m1 + A1m1 + A2m1 + H3m1 + Y3m1 + Z3m1 + A3m1
Nm2 = S1m2 + S2m2 + S3m2 + S4m2 + H1m2 + H2m2 + Y1m2 + Y2m2 + Z1m2 + Z2m2 + A1m2 + A2m2 + H3m2 + Y3m2 + Z3m2 + A3m2
Nm3 = S1m3 + S2m3 + S3m3 + S4m3 + H1m3 + H2m3 + Y1m3 + Y2m3 + Z1m3 + Z2m3 + A1m3 + A2m3 + H3m3 + Y3m3 + Z3m3 + A3m3
N = Nf1 + Nf2 + Nf3 + Nm1 + Nm2 + Nm3  
         
Stot = S1f1 + S2f1 + S1f2 + S2f2 + S1f3 + S2f3 + S1m1 + S2m1 + S3m1 + S4m1 + S1m2 + S2m2 + S3m2 + S4m2 + S1m3 + S2m3 + S3m3 + S4m3
Htot = H1f1 + H2f1 + H3f1 + H1f2 + H2f2 + H3f2 + H1f3 + H2f3 + H3f3 + H1m1 + H2m1 + H3m1 + H1m2 + H2m2 + H3m2 + H1m3 + H2m3 + H3m3
Ytot = Y1f1 + Y2f1 + Y3f1 + Y1f2 + Y2f2 + Y3f2 + Y1f3 + Y2f3 + Y3f3 + Y1m1 + Y2m1 + Y3m1 + Y1m2 + Y2m2 + Y3m2 + Y1m3 + Y2m3 + Y3m3
Ztot = Z1f1 + Z2f1 + Z3f1 + Z1f2 + Z2f2 + Z3f2 + Z1f3 + Z2f3 + Z3f3 + Z1m1 + Z2m1 + Z3m1 + Z1m2 + Z2m2 + Z3m2 + Z1m3 + Z2m3 + Z3m3
Atot = A1f1 + A2f1 + A3f1 + A1f2 + A2f2 + A3f2 + A1f3 + A2f3 + A3f3 + A1m1 + A2m1 + A3m1 + A1m2 + A2m2 + A3m2 + A1m3 + A2m3 + A3m3

S.male = S1m1 + S2m1 + S3m1 + S4m1 + S1m2 + S2m2 + S3m2 + S4m2 + S1m3 + S2m3 + S3m3 + S4m3
S.female = S1f1 + S2f1 + S1f2 + S2f2 + S1f3 + S2f3 

## Proportion aware of HIV+ status
propaware = (H2f1+H3f1+H2f2+H3f2+H2f3+H3f3+H2m1+H3m1+H2m2+H3m2+H2m3+H3m3+
             Y2f1+Y3f1+Y2f2+Y3f2+Y2f3+Y3f3+Y2m1+Y3m1+Y2m2+Y3m2+Y2m3+Y3m3+
             Z2f1+Z3f1+Z2f2+Z3f2+Z2f3+Z3f3+Z2m1+Z3m1+Z2m2+Z3m2+Z2m3+Z3m3+
             A2f1+A3f1+A2f2+A3f2+A2f3+A3f3+A2m1+A3m1+A2m2+A3m2+A2m3+A3m3)/(Htot+Ytot+Ztot+Atot)

# Proportion of susceptibles on PrEP
propPREP = (S2f1+S2f2+S2f3+S2m1+S2m2+S2m3+S4m1+S4m2+S4m3)/Stot

## Prevalence
PRf1 = (H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1)/Nf1
PRf2 = (H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2)/Nf2
PRf3 = (H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3)/Nf3
PRm1 = (H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1)/Nm1
PRm2 = (H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2)/Nm2
PRm3 = (H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/Nm3
         
PRtot = (( H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1
          +H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2
          +H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3
          +H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1
          +H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2
          +H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/N)

## Treatment uptake current time
# start ART in 2004
if(t<=419) {tau3 = 0; tau4 = 0; gamma = 0}
if(t<504) {tau1 = 0; tau2 = 0}
# if(t>552) {tau1=1; tau2=1; tau3=1; tau4=1} # test-and-treat 2015
# Uptake rates
tauC = fun.circ(t,pCIRC,scaleCIRC-pCIRC)*HCT_m*fun.tauT(t,rHCT,pHCT,propaware,rHCT2,pHCT2)
# if(t>540){tauC = fun.circ(t,pCIRC,scaleCIRC-pCIRC)*HCT_m*fun.hct(t,PhctS)}
tauA1 = fun.tauA12(t,tau1,scaleART1-tau1)
tauA2 = fun.tauA12(t,tau2,scaleART2-tau2)
tauA3 = fun.tauA34(t,tau3,scaleART3-tau3)
tauA4 = fun.tauA34(t,tau4,scaleART4-tau4)
tauT = fun.tauT(t,rHCT,pHCT,propaware,rHCT2,pHCT2)
theta = fun.prep(t, pPREP)
tauTS = fun.hct(t,PhctS)

# Baseline
if(t>540){
  tauA1 = fun.tauA12(540,tau1,scaleART1-tau1)
  tauA2 = fun.tauA12(540,tau2,scaleART2-tau2)
  tauA3 = fun.tauA34(540,tau3,scaleART3-tau3)
  tauA4 = fun.tauA34(540,tau4,scaleART4-tau4)
  tauC = fun.circ(540,pCIRC,scaleCIRC-pCIRC)*HCT_m*fun.tauT(t,rHCT,pHCT,propaware,rHCT2,pHCT2)
}

# ART incidence per 100 PY
ARTinc = ((tauA1*H2f1 + tauA2*Y2f1 + tauA3*Z2f1 + tauA4*A2f1 
         + tauA1*H2f2 + tauA2*Y2f2 + tauA3*Z2f2 + tauA4*A2f2 
         + tauA1*H2f3 + tauA2*Y2f3 + tauA3*Z2f3 + tauA4*A2f3
         + ART_m*tauA1*H2m1 + ART_m*tauA2*Y2m1 + ART_m*tauA3*Z2m1 + ART_m*tauA4*A2m1 
         + ART_m*tauA1*H2m2 + ART_m*tauA2*Y2m2 + ART_m*tauA3*Z2m2 + ART_m*tauA4*A2m2 
         + ART_m*tauA1*H2m3 + ART_m*tauA2*Y2m3 + ART_m*tauA3*Z2m3 + ART_m*tauA4*A2m3)/(Htot+Ytot+Ztot+Atot)*100)*12

# ART dropout
ART.drop = (gamma*H3f1 + gamma*H3f2 + gamma*H3f3 +
            gamma*Y3f1 + gamma*Y3f2 + gamma*Y3f3 +
            gamma*Z3f1 + gamma*Z3f2 + gamma*Z3f3 +
            gamma*A3f1 + gamma*A3f2 + gamma*A3f3 +
            gamma*H3m1 + gamma*H3m2 + gamma*H3m3 +
            gamma*Y3m1 + gamma*Y3m2 + gamma*Y3m3 +
            gamma*Z3m1 + gamma*Z3m2 + gamma*Z3m3 +
            gamma*A3m1 + gamma*A3m2 + gamma*A3m3)

# PrEP dropout
PrEP.drop = (om*S2f1 + om*S2f2 + om*S2f3 + om*S2m1 + om*S2m2 + om*S2m3 +
               om*S4m1 + om*S4m2 + om*S4m3)

### MIXING PARAMETERS

## MAIN PARTNERSHIPS
# FEMALES
rho_f1_1_main = (1 - A_main)* ((dm1_main*Nm1)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) + A_main 
rho_f1_2_main = (1 - A_main)* ((dm2_main*Nm2)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
rho_f1_3_main = (1 - A_main)* ((dm3_main*Nm3)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 

rho_f2_1_main = (1 - A_main)* ((dm1_main*Nm1)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
rho_f2_2_main = (1 - A_main)* ((dm2_main*Nm2)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) + A_main
rho_f2_3_main = (1 - A_main)* ((dm3_main*Nm3)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 

rho_f3_1_main = (1 - A_main)* ((dm1_main*Nm1)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
rho_f3_2_main = (1 - A_main)* ((dm2_main*Nm2)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
rho_f3_3_main = (1 - A_main)* ((dm3_main*Nm3)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) + A_main

# MALES
rho_m1_1_main = (1 - A_main)* ((df1_main*Nf1)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) + A_main 
rho_m1_2_main = (1 - A_main)* ((df2_main*Nf2)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
rho_m1_3_main = (1 - A_main)* ((df3_main*Nf3)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 

rho_m2_1_main = (1 - A_main)* ((df1_main*Nf1)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
rho_m2_2_main = (1 - A_main)* ((df2_main*Nf2)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) + A_main
rho_m2_3_main = (1 - A_main)* ((df3_main*Nf3)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 

rho_m3_1_main = (1 - A_main)* ((df1_main*Nf1)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
rho_m3_2_main = (1 - A_main)* ((df2_main*Nf2)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
rho_m3_3_main = (1 - A_main)* ((df3_main*Nf3)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) + A_main

## Balancing the number of sexual partnerships
B11_main = (dm1_main*rho_m1_1_main*Nm1)/(df1_main*rho_f1_1_main*Nf1)
const1 = const1
dm11_main = B11_main^(const1-1) * dm1_main
df11_main = B11_main^const1 * df1_main

B12_main = (dm2_main*rho_m2_1_main*Nm2)/(df1_main*rho_f1_2_main*Nf1)
const1 = const1
dm21_main = B12_main^(const1-1) * dm2_main
df12_main = B12_main^const1 * df1_main

B13_main = (dm3_main*rho_m3_1_main*Nm3)/(df1_main*rho_f1_3_main*Nf1)
const1 = const1
dm31_main = B13_main^(const1-1) * dm3_main
df13_main = B13_main^const1 * df1_main

B21_main = (dm1_main*rho_m1_2_main*Nm1)/(df2_main*rho_f2_1_main*Nf2)
const1 = const1
dm12_main = B21_main^(const1-1) * dm1_main
df21_main = B21_main^const1 * df2_main

B22_main = (dm2_main*rho_m2_2_main*Nm2)/(df2_main*rho_f2_2_main*Nf2)
const1 = const1
dm22_main = B22_main^(const1-1) * dm2_main
df22_main = B22_main^const1 * df2_main

B23_main = (dm3_main*rho_m3_2_main*Nm3)/(df2_main*rho_f2_3_main*Nf2)
const1 = const1
dm32_main = B23_main^(const1-1) * dm3_main
df23_main = B23_main^const1 * df2_main

B31_main = (dm1_main*rho_m1_3_main*Nm1)/(df3_main*rho_f3_1_main*Nf3)
const1 = const1
dm13_main = B31_main^(const1-1) * dm1_main
df31_main = B31_main^const1 * df3_main

B32_main = (dm2_main*rho_m2_3_main*Nm2)/(df3_main*rho_f3_2_main*Nf3)
const1 = const1
dm23_main = B32_main^(const1-1) * dm2_main
df32_main = B32_main^const1 * df3_main

B33_main = (dm3_main*rho_m3_3_main*Nm3)/(df3_main*rho_f3_3_main*Nf3)
const1 = const1
dm33_main = B33_main^(const1-1) * dm3_main
df33_main = B33_main^const1 * df3_main

## REGULAR PARTNERSHIPS
# FEMALES
rho_f2_2_reg = (1 - A_reg)* ((dm2_reg*Nm2)/(dm2_reg*Nm2+dm3_reg*Nm3)) + A_reg 
rho_f2_3_reg = (1 - A_reg)* ((dm3_reg*Nm3)/(dm2_reg*Nm2+dm3_reg*Nm3))

rho_f3_2_reg = (1 - A_reg)* ((dm2_reg*Nm2)/(dm2_reg*Nm2+dm3_reg*Nm3))
rho_f3_3_reg = (1 - A_reg)* ((dm3_reg*Nm3)/(dm2_reg*Nm2+dm3_reg*Nm3)) + A_reg

# MALES
rho_m2_2_reg = (1 - A_reg)* ((df2_reg*Nf2)/(df2_reg*Nf2+df3_reg*Nf3)) + A_reg 
rho_m2_3_reg = (1 - A_reg)* ((df3_reg*Nf3)/(df2_reg*Nf2+df3_reg*Nf3))

rho_m3_2_reg = (1 - A_reg)* ((df2_reg*Nf2)/(df2_reg*Nf2+df3_reg*Nf3))
rho_m3_3_reg = (1 - A_reg)* ((df3_reg*Nf3)/(df2_reg*Nf2+df3_reg*Nf3)) + A_reg


B22_reg = (dm2_reg*rho_m2_2_reg*Nm2)/(df2_reg*rho_f2_2_reg*Nf2)
const2 = const2
dm22_reg = B22_reg^(const2-1) * dm2_reg
df22_reg = B22_reg^const2 * df2_reg

B23_reg = (dm3_reg*rho_m3_2_reg*Nm3)/(df2_reg*rho_f2_3_reg*Nf2)
const2 = const2
dm32_reg = B23_reg^(const2-1) * dm3_reg
df23_reg = B23_reg^const2 * df2_reg

B32_reg = (dm2_reg*rho_m2_3_reg*Nm2)/(df3_reg*rho_f3_2_reg*Nf3)
const2 = const2
dm23_reg = B32_reg^(const2-1) * dm2_reg
df32_reg = B32_reg^const2 * df3_reg

B33_reg = (dm3_reg*rho_m3_3_reg*Nm3)/(df3_reg*rho_f3_3_reg*Nf3)
const2 = const2
dm33_reg = B33_reg^(const2-1) * dm3_reg
df33_reg = B33_reg^const2 * df3_reg


## CASUAL PARTNERSHIPS (always 1)
# FEMALES
rho_f3_3_casual = 1

# MALES
rho_m3_3_casual = 1

B33_casual = (dm3_casual*rho_m3_3_casual*Nm3)/(df3_casual*rho_f3_3_casual*Nf3)
const3 = const3
dm33_casual = B33_casual^(const3-1) * dm3_casual
df33_casual = B33_casual^const3 * df3_casual


##### PROBABILITY OF TRANSMISSION

## Condom use at current time
cond_main = fun.cond(t,f_main,scaleCOND1-f_main)
cond_reg = fun.cond(t,f_reg,scaleCOND2-f_reg)
cond_casual = fun.cond(t,f_casual,scaleCOND3-f_casual)
# Baseline
if(t>540){
  cond_main = fun.cond(540,f_main,scaleCOND1-f_main)
  cond_reg = fun.cond(540,f_reg,scaleCOND2-f_reg)
  cond_casual = fun.cond(540,f_casual,scaleCOND3-f_casual)
}

### Susceptibles not on PrEP/not circumcised ###

## Partner in LV stage, not on HCT/ART

D1G1_f11_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf11_main) * (1-beta_f)^((1-cond_main)*actsf11_main))
D1G1_f12_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf12_main) * (1-beta_f)^((1-cond_main)*actsf12_main))
D1G1_f13_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf13_main) * (1-beta_f)^((1-cond_main)*actsf13_main))
D1G1_f21_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf21_main) * (1-beta_f)^((1-cond_main)*actsf21_main))
D1G1_f22_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf22_main) * (1-beta_f)^((1-cond_main)*actsf22_main))
D1G1_f23_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf23_main) * (1-beta_f)^((1-cond_main)*actsf23_main))
D1G1_f31_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf31_main) * (1-beta_f)^((1-cond_main)*actsf31_main))
D1G1_f32_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf32_main) * (1-beta_f)^((1-cond_main)*actsf32_main))
D1G1_f33_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf33_main) * (1-beta_f)^((1-cond_main)*actsf33_main))
D1G1_m11_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm11_main) * (1-beta_m)^((1-cond_main)*actsm11_main))
D1G1_m12_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm12_main) * (1-beta_m)^((1-cond_main)*actsm12_main))
D1G1_m13_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm13_main) * (1-beta_m)^((1-cond_main)*actsm13_main))
D1G1_m21_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm21_main) * (1-beta_m)^((1-cond_main)*actsm21_main))
D1G1_m22_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm22_main) * (1-beta_m)^((1-cond_main)*actsm22_main))
D1G1_m23_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm23_main) * (1-beta_m)^((1-cond_main)*actsm23_main))
D1G1_m31_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm31_main) * (1-beta_m)^((1-cond_main)*actsm31_main))
D1G1_m32_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm32_main) * (1-beta_m)^((1-cond_main)*actsm32_main))
D1G1_m33_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm33_main) * (1-beta_m)^((1-cond_main)*actsm33_main))
D1G1_f22_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta_f)^((1-cond_reg)*actsf22_reg))
D1G1_f23_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta_f)^((1-cond_reg)*actsf23_reg))
D1G1_f32_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta_f)^((1-cond_reg)*actsf32_reg))
D1G1_f33_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta_f)^((1-cond_reg)*actsf33_reg))
D1G1_m22_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta_m)^((1-cond_reg)*actsm22_reg))
D1G1_m23_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta_m)^((1-cond_reg)*actsm23_reg))
D1G1_m32_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta_m)^((1-cond_reg)*actsm32_reg))
D1G1_m33_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta_m)^((1-cond_reg)*actsm33_reg))
D1G1_f33_casual = 1 - ((1-beta_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta_f)^((1-cond_casual)*actsf33_casual))
D1G1_m33_casual = 1 - ((1-beta_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta_m)^((1-cond_casual)*actsm33_casual))

## Partner in LV stage, on HCT

D1G2_f11_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf11_main))
D1G2_f12_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf12_main))
D1G2_f13_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf13_main))
D1G2_f21_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf21_main))
D1G2_f22_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf22_main))
D1G2_f23_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf23_main))
D1G2_f31_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf31_main))
D1G2_f32_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf32_main))
D1G2_f33_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf33_main))
D1G2_m11_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
D1G2_m12_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
D1G2_m13_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
D1G2_m21_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
D1G2_m22_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
D1G2_m23_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
D1G2_m31_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
D1G2_m32_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
D1G2_m33_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
D1G2_f22_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf22_reg))
D1G2_f23_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf23_reg))
D1G2_f32_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf32_reg))
D1G2_f33_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf33_reg))
D1G2_m22_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
D1G2_m23_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
D1G2_m32_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
D1G2_m33_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
D1G2_f33_casual = 1 - ((1-beta_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta_f)^((1-cond_casual*HCT_cond)*HCT_red*actsf33_casual))
D1G2_m33_casual = 1 - ((1-beta_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))

## Partner in LV stage, on ART

D1G3_f11_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf11_main) * (1-beta4_f)^((1-cond_main)*actsf11_main))
D1G3_f12_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf12_main) * (1-beta4_f)^((1-cond_main)*actsf12_main))
D1G3_f13_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf13_main) * (1-beta4_f)^((1-cond_main)*actsf13_main))
D1G3_f21_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf21_main) * (1-beta4_f)^((1-cond_main)*actsf21_main))
D1G3_f22_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf22_main) * (1-beta4_f)^((1-cond_main)*actsf22_main))
D1G3_f23_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf23_main) * (1-beta4_f)^((1-cond_main)*actsf23_main))
D1G3_f31_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf31_main) * (1-beta4_f)^((1-cond_main)*actsf31_main))
D1G3_f32_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf32_main) * (1-beta4_f)^((1-cond_main)*actsf32_main))
D1G3_f33_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf33_main) * (1-beta4_f)^((1-cond_main)*actsf33_main))
D1G3_m11_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm11_main) * (1-beta4_m)^((1-cond_main)*actsm11_main))
D1G3_m12_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm12_main) * (1-beta4_m)^((1-cond_main)*actsm12_main))
D1G3_m13_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm13_main) * (1-beta4_m)^((1-cond_main)*actsm13_main))
D1G3_m21_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm21_main) * (1-beta4_m)^((1-cond_main)*actsm21_main))
D1G3_m22_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm22_main) * (1-beta4_m)^((1-cond_main)*actsm22_main))
D1G3_m23_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm23_main) * (1-beta4_m)^((1-cond_main)*actsm23_main))
D1G3_m31_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm31_main) * (1-beta4_m)^((1-cond_main)*actsm31_main))
D1G3_m32_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm32_main) * (1-beta4_m)^((1-cond_main)*actsm32_main))
D1G3_m33_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm33_main) * (1-beta4_m)^((1-cond_main)*actsm33_main))
D1G3_f22_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta4_f)^((1-cond_reg)*actsf22_reg))
D1G3_f23_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta4_f)^((1-cond_reg)*actsf23_reg))
D1G3_f32_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta4_f)^((1-cond_reg)*actsf32_reg))
D1G3_f33_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta4_f)^((1-cond_reg)*actsf33_reg))
D1G3_m22_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta4_m)^((1-cond_reg)*actsm22_reg))
D1G3_m23_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta4_m)^((1-cond_reg)*actsm23_reg))
D1G3_m32_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta4_m)^((1-cond_reg)*actsm32_reg))
D1G3_m33_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta4_m)^((1-cond_reg)*actsm33_reg))
D1G3_f33_casual = 1 - ((1-beta4_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta4_f)^((1-cond_casual)*actsf33_casual))
D1G3_m33_casual = 1 - ((1-beta4_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta4_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in initial HV stage, not on HCT/ART

D1P1_f11_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf11_main) * (1-beta1_f)^((1-cond_main)*actsf11_main))
D1P1_f12_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf12_main) * (1-beta1_f)^((1-cond_main)*actsf12_main))
D1P1_f13_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf13_main) * (1-beta1_f)^((1-cond_main)*actsf13_main))
D1P1_f21_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf21_main) * (1-beta1_f)^((1-cond_main)*actsf21_main))
D1P1_f22_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf22_main) * (1-beta1_f)^((1-cond_main)*actsf22_main))
D1P1_f23_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf23_main) * (1-beta1_f)^((1-cond_main)*actsf23_main))
D1P1_f31_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf31_main) * (1-beta1_f)^((1-cond_main)*actsf31_main))
D1P1_f32_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf32_main) * (1-beta1_f)^((1-cond_main)*actsf32_main))
D1P1_f33_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf33_main) * (1-beta1_f)^((1-cond_main)*actsf33_main))
D1P1_m11_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm11_main) * (1-beta1_m)^((1-cond_main)*actsm11_main))
D1P1_m12_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm12_main) * (1-beta1_m)^((1-cond_main)*actsm12_main))
D1P1_m13_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm13_main) * (1-beta1_m)^((1-cond_main)*actsm13_main))
D1P1_m21_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm21_main) * (1-beta1_m)^((1-cond_main)*actsm21_main))
D1P1_m22_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm22_main) * (1-beta1_m)^((1-cond_main)*actsm22_main))
D1P1_m23_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm23_main) * (1-beta1_m)^((1-cond_main)*actsm23_main))
D1P1_m31_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm31_main) * (1-beta1_m)^((1-cond_main)*actsm31_main))
D1P1_m32_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm32_main) * (1-beta1_m)^((1-cond_main)*actsm32_main))
D1P1_m33_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm33_main) * (1-beta1_m)^((1-cond_main)*actsm33_main))
D1P1_f22_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta1_f)^((1-cond_reg)*actsf22_reg))
D1P1_f23_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta1_f)^((1-cond_reg)*actsf23_reg))
D1P1_f32_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta1_f)^((1-cond_reg)*actsf32_reg))
D1P1_f33_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta1_f)^((1-cond_reg)*actsf33_reg))
D1P1_m22_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta1_m)^((1-cond_reg)*actsm22_reg))
D1P1_m23_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta1_m)^((1-cond_reg)*actsm23_reg))
D1P1_m32_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta1_m)^((1-cond_reg)*actsm32_reg))
D1P1_m33_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta1_m)^((1-cond_reg)*actsm33_reg))
D1P1_f33_casual = 1 - ((1-beta1_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta1_f)^((1-cond_casual)*actsf33_casual))
D1P1_m33_casual = 1 - ((1-beta1_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta1_m)^((1-cond_casual)*actsm33_casual))

## Partner in initial HV stage, on HCT

D1P2_f11_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf11_main))
D1P2_f12_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf12_main))
D1P2_f13_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf13_main))
D1P2_f21_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf21_main))
D1P2_f22_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf22_main))
D1P2_f23_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf23_main))
D1P2_f31_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf31_main))
D1P2_f32_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf32_main))
D1P2_f33_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf33_main))
D1P2_m11_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
D1P2_m12_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
D1P2_m13_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
D1P2_m21_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
D1P2_m22_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
D1P2_m23_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
D1P2_m31_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
D1P2_m32_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
D1P2_m33_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
D1P2_f22_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf22_reg))
D1P2_f23_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf23_reg))
D1P2_f32_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf32_reg))
D1P2_f33_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf33_reg))
D1P2_m22_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
D1P2_m23_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
D1P2_m32_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
D1P2_m33_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
D1P2_f33_casual = 1 - ((1-beta1_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta1_f)^((1-cond_casual*HCT_cond)*HCT_red*actsf33_casual))
D1P2_m33_casual = 1 - ((1-beta1_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta1_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))

## Partner in initial HV stage, on ART

D1P3_f11_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf11_main) * (1-beta5_f)^((1-cond_main)*actsf11_main))
D1P3_f12_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf12_main) * (1-beta5_f)^((1-cond_main)*actsf12_main))
D1P3_f13_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf13_main) * (1-beta5_f)^((1-cond_main)*actsf13_main))
D1P3_f21_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf21_main) * (1-beta5_f)^((1-cond_main)*actsf21_main))
D1P3_f22_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf22_main) * (1-beta5_f)^((1-cond_main)*actsf22_main))
D1P3_f23_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf23_main) * (1-beta5_f)^((1-cond_main)*actsf23_main))
D1P3_f31_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf31_main) * (1-beta5_f)^((1-cond_main)*actsf31_main))
D1P3_f32_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf32_main) * (1-beta5_f)^((1-cond_main)*actsf32_main))
D1P3_f33_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf33_main) * (1-beta5_f)^((1-cond_main)*actsf33_main))
D1P3_m11_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm11_main) * (1-beta5_m)^((1-cond_main)*actsm11_main))
D1P3_m12_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm12_main) * (1-beta5_m)^((1-cond_main)*actsm12_main))
D1P3_m13_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm13_main) * (1-beta5_m)^((1-cond_main)*actsm13_main))
D1P3_m21_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm21_main) * (1-beta5_m)^((1-cond_main)*actsm21_main))
D1P3_m22_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm22_main) * (1-beta5_m)^((1-cond_main)*actsm22_main))
D1P3_m23_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm23_main) * (1-beta5_m)^((1-cond_main)*actsm23_main))
D1P3_m31_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm31_main) * (1-beta5_m)^((1-cond_main)*actsm31_main))
D1P3_m32_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm32_main) * (1-beta5_m)^((1-cond_main)*actsm32_main))
D1P3_m33_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm33_main) * (1-beta5_m)^((1-cond_main)*actsm33_main))
D1P3_f22_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta5_f)^((1-cond_reg)*actsf22_reg))
D1P3_f23_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta5_f)^((1-cond_reg)*actsf23_reg))
D1P3_f32_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta5_f)^((1-cond_reg)*actsf32_reg))
D1P3_f33_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta5_f)^((1-cond_reg)*actsf33_reg))
D1P3_m22_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta5_m)^((1-cond_reg)*actsm22_reg))
D1P3_m23_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta5_m)^((1-cond_reg)*actsm23_reg))
D1P3_m32_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta5_m)^((1-cond_reg)*actsm32_reg))
D1P3_m33_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta5_m)^((1-cond_reg)*actsm33_reg))
D1P3_f33_casual = 1 - ((1-beta5_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta5_f)^((1-cond_casual)*actsf33_casual))
D1P3_m33_casual = 1 - ((1-beta5_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta5_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in preAIDS HV stage, not on HCT/ART

D1X1_f11_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf11_main) * (1-beta2_f)^((1-cond_main)*actsf11_main))
D1X1_f12_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf12_main) * (1-beta2_f)^((1-cond_main)*actsf12_main))
D1X1_f13_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf13_main) * (1-beta2_f)^((1-cond_main)*actsf13_main))
D1X1_f21_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf21_main) * (1-beta2_f)^((1-cond_main)*actsf21_main))
D1X1_f22_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf22_main) * (1-beta2_f)^((1-cond_main)*actsf22_main))
D1X1_f23_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf23_main) * (1-beta2_f)^((1-cond_main)*actsf23_main))
D1X1_f31_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf31_main) * (1-beta2_f)^((1-cond_main)*actsf31_main))
D1X1_f32_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf32_main) * (1-beta2_f)^((1-cond_main)*actsf32_main))
D1X1_f33_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf33_main) * (1-beta2_f)^((1-cond_main)*actsf33_main))
D1X1_m11_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm11_main) * (1-beta2_m)^((1-cond_main)*actsm11_main))
D1X1_m12_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm12_main) * (1-beta2_m)^((1-cond_main)*actsm12_main))
D1X1_m13_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm13_main) * (1-beta2_m)^((1-cond_main)*actsm13_main))
D1X1_m21_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm21_main) * (1-beta2_m)^((1-cond_main)*actsm21_main))
D1X1_m22_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm22_main) * (1-beta2_m)^((1-cond_main)*actsm22_main))
D1X1_m23_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm23_main) * (1-beta2_m)^((1-cond_main)*actsm23_main))
D1X1_m31_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm31_main) * (1-beta2_m)^((1-cond_main)*actsm31_main))
D1X1_m32_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm32_main) * (1-beta2_m)^((1-cond_main)*actsm32_main))
D1X1_m33_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm33_main) * (1-beta2_m)^((1-cond_main)*actsm33_main))
D1X1_f22_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta2_f)^((1-cond_reg)*actsf22_reg))
D1X1_f23_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta2_f)^((1-cond_reg)*actsf23_reg))
D1X1_f32_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta2_f)^((1-cond_reg)*actsf32_reg))
D1X1_f33_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta2_f)^((1-cond_reg)*actsf33_reg))
D1X1_m22_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta2_m)^((1-cond_reg)*actsm22_reg))
D1X1_m23_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta2_m)^((1-cond_reg)*actsm23_reg))
D1X1_m32_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta2_m)^((1-cond_reg)*actsm32_reg))
D1X1_m33_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta2_m)^((1-cond_reg)*actsm33_reg))
D1X1_f33_casual = 1 - ((1-beta2_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta2_f)^((1-cond_casual)*actsf33_casual))
D1X1_m33_casual = 1 - ((1-beta2_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta2_m)^((1-cond_casual)*actsm33_casual))

## Partner in preAIDS HV stage, on HCT

D1X2_f11_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf11_main))
D1X2_f12_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf12_main))
D1X2_f13_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf13_main))
D1X2_f21_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf21_main))
D1X2_f22_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf22_main))
D1X2_f23_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf23_main))
D1X2_f31_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf31_main))
D1X2_f32_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf32_main))
D1X2_f33_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf33_main))
D1X2_m11_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
D1X2_m12_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
D1X2_m13_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
D1X2_m21_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
D1X2_m22_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
D1X2_m23_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
D1X2_m31_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
D1X2_m32_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
D1X2_m33_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
D1X2_f22_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf22_reg))
D1X2_f23_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf23_reg))
D1X2_f32_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf32_reg))
D1X2_f33_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf33_reg))
D1X2_m22_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
D1X2_m23_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
D1X2_m32_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
D1X2_m33_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
D1X2_f33_casual = 1 - ((1-beta2_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta2_f)^((1-cond_casual*HCT_cond)*HCT_red*actsf33_casual))
D1X2_m33_casual = 1 - ((1-beta2_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta2_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))

## Partner in preAIDS HV stage, on ART

D1X3_f11_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf11_main) * (1-beta6_f)^((1-cond_main)*actsf11_main))
D1X3_f12_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf12_main) * (1-beta6_f)^((1-cond_main)*actsf12_main))
D1X3_f13_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf13_main) * (1-beta6_f)^((1-cond_main)*actsf13_main))
D1X3_f21_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf21_main) * (1-beta6_f)^((1-cond_main)*actsf21_main))
D1X3_f22_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf22_main) * (1-beta6_f)^((1-cond_main)*actsf22_main))
D1X3_f23_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf23_main) * (1-beta6_f)^((1-cond_main)*actsf23_main))
D1X3_f31_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf31_main) * (1-beta6_f)^((1-cond_main)*actsf31_main))
D1X3_f32_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf32_main) * (1-beta6_f)^((1-cond_main)*actsf32_main))
D1X3_f33_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf33_main) * (1-beta6_f)^((1-cond_main)*actsf33_main))
D1X3_m11_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm11_main) * (1-beta6_m)^((1-cond_main)*actsm11_main))
D1X3_m12_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm12_main) * (1-beta6_m)^((1-cond_main)*actsm12_main))
D1X3_m13_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm13_main) * (1-beta6_m)^((1-cond_main)*actsm13_main))
D1X3_m21_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm21_main) * (1-beta6_m)^((1-cond_main)*actsm21_main))
D1X3_m22_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm22_main) * (1-beta6_m)^((1-cond_main)*actsm22_main))
D1X3_m23_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm23_main) * (1-beta6_m)^((1-cond_main)*actsm23_main))
D1X3_m31_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm31_main) * (1-beta6_m)^((1-cond_main)*actsm31_main))
D1X3_m32_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm32_main) * (1-beta6_m)^((1-cond_main)*actsm32_main))
D1X3_m33_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm33_main) * (1-beta6_m)^((1-cond_main)*actsm33_main))
D1X3_f22_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta6_f)^((1-cond_reg)*actsf22_reg))
D1X3_f23_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta6_f)^((1-cond_reg)*actsf23_reg))
D1X3_f32_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta6_f)^((1-cond_reg)*actsf32_reg))
D1X3_f33_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta6_f)^((1-cond_reg)*actsf33_reg))
D1X3_m22_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta6_m)^((1-cond_reg)*actsm22_reg))
D1X3_m23_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta6_m)^((1-cond_reg)*actsm23_reg))
D1X3_m32_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta6_m)^((1-cond_reg)*actsm32_reg))
D1X3_m33_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta6_m)^((1-cond_reg)*actsm33_reg))
D1X3_f33_casual = 1 - ((1-beta6_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta6_f)^((1-cond_casual)*actsf33_casual))
D1X3_m33_casual = 1 - ((1-beta6_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta6_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in AIDS stage, not on HCT/ART

D1K1_f11_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf11_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf11_main))
D1K1_f12_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf12_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf12_main))
D1K1_f13_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf13_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf13_main))
D1K1_f21_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf21_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf21_main))
D1K1_f22_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf22_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf22_main))
D1K1_f23_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf23_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf23_main))
D1K1_f31_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf31_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf31_main))
D1K1_f32_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf32_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf32_main))
D1K1_f33_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf33_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf33_main))
D1K1_m11_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm11_main))
D1K1_m12_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm12_main))
D1K1_m13_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm13_main))
D1K1_m21_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm21_main))
D1K1_m22_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm22_main))
D1K1_m23_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm23_main))
D1K1_m31_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm31_main))
D1K1_m32_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm32_main))
D1K1_m33_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm33_main))
D1K1_f22_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf22_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf22_reg))
D1K1_f23_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf23_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf23_reg))
D1K1_f32_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf32_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf32_reg))
D1K1_f33_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf33_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf33_reg))
D1K1_m22_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
D1K1_m23_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
D1K1_m32_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
D1K1_m33_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
D1K1_f33_casual = 1 - ((1-beta3_f*cond_eff)^(cond_casual*AIDS_red*actsf33_casual) * (1-beta3_f)^((1-cond_casual)*AIDS_red*actsf33_casual))
D1K1_m33_casual = 1 - ((1-beta3_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta3_m)^((1-cond_casual)*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on HCT

D1K2_f11_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf11_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf11_main))
D1K2_f12_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf12_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf12_main))
D1K2_f13_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf13_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf13_main))
D1K2_f21_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf21_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf21_main))
D1K2_f22_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf22_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf22_main))
D1K2_f23_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf23_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf23_main))
D1K2_f31_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf31_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf31_main))
D1K2_f32_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf32_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf32_main))
D1K2_f33_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf33_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf33_main))
D1K2_m11_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm11_main))
D1K2_m12_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm12_main))
D1K2_m13_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm13_main))
D1K2_m21_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm21_main))
D1K2_m22_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm22_main))
D1K2_m23_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm23_main))
D1K2_m31_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm31_main))
D1K2_m32_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm32_main))
D1K2_m33_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm33_main))
D1K2_f22_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf22_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf22_reg))
D1K2_f23_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf23_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf23_reg))
D1K2_f32_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf32_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf32_reg))
D1K2_f33_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf33_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf33_reg))
D1K2_m22_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm22_reg))
D1K2_m23_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm23_reg))
D1K2_m32_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm32_reg))
D1K2_m33_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm33_reg))
D1K2_f33_casual = 1 - ((1-beta3_f*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsf33_casual) * (1-beta3_f)^((1-cond_casual*HCT_cond)*AIDS_red*actsf33_casual))
D1K2_m33_casual = 1 - ((1-beta3_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta3_m)^((1-cond_casual*HCT_cond)*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on ART

D1K3_f11_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf11_main) * (1-beta7_f)^((1-cond_main)*actsf11_main))
D1K3_f12_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf12_main) * (1-beta7_f)^((1-cond_main)*actsf12_main))
D1K3_f13_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf13_main) * (1-beta7_f)^((1-cond_main)*actsf13_main))
D1K3_f21_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf21_main) * (1-beta7_f)^((1-cond_main)*actsf21_main))
D1K3_f22_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf22_main) * (1-beta7_f)^((1-cond_main)*actsf22_main))
D1K3_f23_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf23_main) * (1-beta7_f)^((1-cond_main)*actsf23_main))
D1K3_f31_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf31_main) * (1-beta7_f)^((1-cond_main)*actsf31_main))
D1K3_f32_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf32_main) * (1-beta7_f)^((1-cond_main)*actsf32_main))
D1K3_f33_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf33_main) * (1-beta7_f)^((1-cond_main)*actsf33_main))
D1K3_m11_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm11_main) * (1-beta7_m)^((1-cond_main)*actsm11_main))
D1K3_m12_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm12_main) * (1-beta7_m)^((1-cond_main)*actsm12_main))
D1K3_m13_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm13_main) * (1-beta7_m)^((1-cond_main)*actsm13_main))
D1K3_m21_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm21_main) * (1-beta7_m)^((1-cond_main)*actsm21_main))
D1K3_m22_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm22_main) * (1-beta7_m)^((1-cond_main)*actsm22_main))
D1K3_m23_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm23_main) * (1-beta7_m)^((1-cond_main)*actsm23_main))
D1K3_m31_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm31_main) * (1-beta7_m)^((1-cond_main)*actsm31_main))
D1K3_m32_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm32_main) * (1-beta7_m)^((1-cond_main)*actsm32_main))
D1K3_m33_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm33_main) * (1-beta7_m)^((1-cond_main)*actsm33_main))
D1K3_f22_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta7_f)^((1-cond_reg)*actsf22_reg))
D1K3_f23_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta7_f)^((1-cond_reg)*actsf23_reg))
D1K3_f32_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta7_f)^((1-cond_reg)*actsf32_reg))
D1K3_f33_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta7_f)^((1-cond_reg)*actsf33_reg))
D1K3_m22_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta7_m)^((1-cond_reg)*actsm22_reg))
D1K3_m23_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta7_m)^((1-cond_reg)*actsm23_reg))
D1K3_m32_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta7_m)^((1-cond_reg)*actsm32_reg))
D1K3_m33_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta7_m)^((1-cond_reg)*actsm33_reg))
D1K3_f33_casual = 1 - ((1-beta7_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta7_f)^((1-cond_casual)*actsf33_casual))
D1K3_m33_casual = 1 - ((1-beta7_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta7_m)^((1-cond_casual)*actsm33_casual))

### Susceptibles on PrEP/not circumcised ###

## Partner in LV stage, not on HCT/ART

D2G1_f11_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf11_main) * (1-beta8_f)^((1-cond_main)*actsf11_main))
D2G1_f12_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf12_main) * (1-beta8_f)^((1-cond_main)*actsf12_main))
D2G1_f13_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf13_main) * (1-beta8_f)^((1-cond_main)*actsf13_main))
D2G1_f21_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf21_main) * (1-beta8_f)^((1-cond_main)*actsf21_main))
D2G1_f22_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf22_main) * (1-beta8_f)^((1-cond_main)*actsf22_main))
D2G1_f23_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf23_main) * (1-beta8_f)^((1-cond_main)*actsf23_main))
D2G1_f31_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf31_main) * (1-beta8_f)^((1-cond_main)*actsf31_main))
D2G1_f32_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf32_main) * (1-beta8_f)^((1-cond_main)*actsf32_main))
D2G1_f33_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf33_main) * (1-beta8_f)^((1-cond_main)*actsf33_main))
D2G1_m11_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm11_main) * (1-beta8_m)^((1-cond_main)*actsm11_main))
D2G1_m12_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm12_main) * (1-beta8_m)^((1-cond_main)*actsm12_main))
D2G1_m13_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm13_main) * (1-beta8_m)^((1-cond_main)*actsm13_main))
D2G1_m21_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm21_main) * (1-beta8_m)^((1-cond_main)*actsm21_main))
D2G1_m22_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm22_main) * (1-beta8_m)^((1-cond_main)*actsm22_main))
D2G1_m23_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm23_main) * (1-beta8_m)^((1-cond_main)*actsm23_main))
D2G1_m31_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm31_main) * (1-beta8_m)^((1-cond_main)*actsm31_main))
D2G1_m32_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm32_main) * (1-beta8_m)^((1-cond_main)*actsm32_main))
D2G1_m33_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm33_main) * (1-beta8_m)^((1-cond_main)*actsm33_main))
D2G1_f22_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta8_f)^((1-cond_reg)*actsf22_reg))
D2G1_f23_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta8_f)^((1-cond_reg)*actsf23_reg))
D2G1_f32_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta8_f)^((1-cond_reg)*actsf32_reg))
D2G1_f33_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta8_f)^((1-cond_reg)*actsf33_reg))
D2G1_m22_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta8_m)^((1-cond_reg)*actsm22_reg))
D2G1_m23_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta8_m)^((1-cond_reg)*actsm23_reg))
D2G1_m32_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta8_m)^((1-cond_reg)*actsm32_reg))
D2G1_m33_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta8_m)^((1-cond_reg)*actsm33_reg))
D2G1_f33_casual = 1 - ((1-beta8_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta8_f)^((1-cond_casual)*actsf33_casual))
D2G1_m33_casual = 1 - ((1-beta8_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta8_m)^((1-cond_casual)*actsm33_casual))

## Partner in LV stage, on HCT

D2G2_f11_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf11_main))
D2G2_f12_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf12_main))
D2G2_f13_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf13_main))
D2G2_f21_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf21_main))
D2G2_f22_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf22_main))
D2G2_f23_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf23_main))
D2G2_f31_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf31_main))
D2G2_f32_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf32_main))
D2G2_f33_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf33_main))
D2G2_m11_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
D2G2_m12_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
D2G2_m13_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
D2G2_m21_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
D2G2_m22_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
D2G2_m23_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
D2G2_m31_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
D2G2_m32_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
D2G2_m33_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
D2G2_f22_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf22_reg))
D2G2_f23_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf23_reg))
D2G2_f32_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf32_reg))
D2G2_f33_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf33_reg))
D2G2_m22_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
D2G2_m23_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
D2G2_m32_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
D2G2_m33_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
D2G2_f33_casual = 1 - ((1-beta8_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta8_f)^((1-cond_casual)*HCT_cond*HCT_red*actsf33_casual))
D2G2_m33_casual = 1 - ((1-beta8_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta8_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))

## Partner in LV stage, on ART

D2G3_f11_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf11_main) * (1-beta16_f)^((1-cond_main)*actsf11_main))
D2G3_f12_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf12_main) * (1-beta16_f)^((1-cond_main)*actsf12_main))
D2G3_f13_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf13_main) * (1-beta16_f)^((1-cond_main)*actsf13_main))
D2G3_f21_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf21_main) * (1-beta16_f)^((1-cond_main)*actsf21_main))
D2G3_f22_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf22_main) * (1-beta16_f)^((1-cond_main)*actsf22_main))
D2G3_f23_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf23_main) * (1-beta16_f)^((1-cond_main)*actsf23_main))
D2G3_f31_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf31_main) * (1-beta16_f)^((1-cond_main)*actsf31_main))
D2G3_f32_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf32_main) * (1-beta16_f)^((1-cond_main)*actsf32_main))
D2G3_f33_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf33_main) * (1-beta16_f)^((1-cond_main)*actsf33_main))
D2G3_m11_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm11_main) * (1-beta16_m)^((1-cond_main)*actsm11_main))
D2G3_m12_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm12_main) * (1-beta16_m)^((1-cond_main)*actsm12_main))
D2G3_m13_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm13_main) * (1-beta16_m)^((1-cond_main)*actsm13_main))
D2G3_m21_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm21_main) * (1-beta16_m)^((1-cond_main)*actsm21_main))
D2G3_m22_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm22_main) * (1-beta16_m)^((1-cond_main)*actsm22_main))
D2G3_m23_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm23_main) * (1-beta16_m)^((1-cond_main)*actsm23_main))
D2G3_m31_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm31_main) * (1-beta16_m)^((1-cond_main)*actsm31_main))
D2G3_m32_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm32_main) * (1-beta16_m)^((1-cond_main)*actsm32_main))
D2G3_m33_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm33_main) * (1-beta16_m)^((1-cond_main)*actsm33_main))
D2G3_f22_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta16_f)^((1-cond_reg)*actsf22_reg))
D2G3_f23_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta16_f)^((1-cond_reg)*actsf23_reg))
D2G3_f32_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta16_f)^((1-cond_reg)*actsf32_reg))
D2G3_f33_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta16_f)^((1-cond_reg)*actsf33_reg))
D2G3_m22_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta16_m)^((1-cond_reg)*actsm22_reg))
D2G3_m23_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta16_m)^((1-cond_reg)*actsm23_reg))
D2G3_m32_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta16_m)^((1-cond_reg)*actsm32_reg))
D2G3_m33_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta16_m)^((1-cond_reg)*actsm33_reg))
D2G3_f33_casual = 1 - ((1-beta16_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta16_f)^((1-cond_casual)*actsf33_casual))
D2G3_m33_casual = 1 - ((1-beta16_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta16_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in initial HV stage, not on HCT/ART

D2P1_f11_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf11_main) * (1-beta9_f)^((1-cond_main)*actsf11_main))
D2P1_f12_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf12_main) * (1-beta9_f)^((1-cond_main)*actsf12_main))
D2P1_f13_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf13_main) * (1-beta9_f)^((1-cond_main)*actsf13_main))
D2P1_f21_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf21_main) * (1-beta9_f)^((1-cond_main)*actsf21_main))
D2P1_f22_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf22_main) * (1-beta9_f)^((1-cond_main)*actsf22_main))
D2P1_f23_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf23_main) * (1-beta9_f)^((1-cond_main)*actsf23_main))
D2P1_f31_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf31_main) * (1-beta9_f)^((1-cond_main)*actsf31_main))
D2P1_f32_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf32_main) * (1-beta9_f)^((1-cond_main)*actsf32_main))
D2P1_f33_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf33_main) * (1-beta9_f)^((1-cond_main)*actsf33_main))
D2P1_m11_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm11_main) * (1-beta9_m)^((1-cond_main)*actsm11_main))
D2P1_m12_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm12_main) * (1-beta9_m)^((1-cond_main)*actsm12_main))
D2P1_m13_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm13_main) * (1-beta9_m)^((1-cond_main)*actsm13_main))
D2P1_m21_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm21_main) * (1-beta9_m)^((1-cond_main)*actsm21_main))
D2P1_m22_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm22_main) * (1-beta9_m)^((1-cond_main)*actsm22_main))
D2P1_m23_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm23_main) * (1-beta9_m)^((1-cond_main)*actsm23_main))
D2P1_m31_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm31_main) * (1-beta9_m)^((1-cond_main)*actsm31_main))
D2P1_m32_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm32_main) * (1-beta9_m)^((1-cond_main)*actsm32_main))
D2P1_m33_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm33_main) * (1-beta9_m)^((1-cond_main)*actsm33_main))
D2P1_f22_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta9_f)^((1-cond_reg)*actsf22_reg))
D2P1_f23_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta9_f)^((1-cond_reg)*actsf23_reg))
D2P1_f32_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta9_f)^((1-cond_reg)*actsf32_reg))
D2P1_f33_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta9_f)^((1-cond_reg)*actsf33_reg))
D2P1_m22_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta9_m)^((1-cond_reg)*actsm22_reg))
D2P1_m23_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta9_m)^((1-cond_reg)*actsm23_reg))
D2P1_m32_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta9_m)^((1-cond_reg)*actsm32_reg))
D2P1_m33_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta9_m)^((1-cond_reg)*actsm33_reg))
D2P1_f33_casual = 1 - ((1-beta9_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta9_f)^((1-cond_casual)*actsf33_casual))
D2P1_m33_casual = 1 - ((1-beta9_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta9_m)^((1-cond_casual)*actsm33_casual))

## Partner in initial HV stage, on HCT

D2P2_f11_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf11_main))
D2P2_f12_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf12_main))
D2P2_f13_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf13_main))
D2P2_f21_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf21_main))
D2P2_f22_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf22_main))
D2P2_f23_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf23_main))
D2P2_f31_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf31_main))
D2P2_f32_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf32_main))
D2P2_f33_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf33_main))
D2P2_m11_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
D2P2_m12_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
D2P2_m13_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
D2P2_m21_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
D2P2_m22_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
D2P2_m23_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
D2P2_m31_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
D2P2_m32_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
D2P2_m33_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
D2P2_f22_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf22_reg))
D2P2_f23_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf23_reg))
D2P2_f32_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf32_reg))
D2P2_f33_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf33_reg))
D2P2_m22_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
D2P2_m23_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
D2P2_m32_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
D2P2_m33_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
D2P2_f33_casual = 1 - ((1-beta9_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta9_f)^((1-cond_casual)*HCT_cond*HCT_red*actsf33_casual))
D2P2_m33_casual = 1 - ((1-beta9_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta9_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))

## Partner in initial HV stage, on ART

D2P3_f11_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf11_main) * (1-beta17_f)^((1-cond_main)*actsf11_main))
D2P3_f12_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf12_main) * (1-beta17_f)^((1-cond_main)*actsf12_main))
D2P3_f13_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf13_main) * (1-beta17_f)^((1-cond_main)*actsf13_main))
D2P3_f21_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf21_main) * (1-beta17_f)^((1-cond_main)*actsf21_main))
D2P3_f22_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf22_main) * (1-beta17_f)^((1-cond_main)*actsf22_main))
D2P3_f23_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf23_main) * (1-beta17_f)^((1-cond_main)*actsf23_main))
D2P3_f31_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf31_main) * (1-beta17_f)^((1-cond_main)*actsf31_main))
D2P3_f32_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf32_main) * (1-beta17_f)^((1-cond_main)*actsf32_main))
D2P3_f33_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf33_main) * (1-beta17_f)^((1-cond_main)*actsf33_main))
D2P3_m11_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm11_main) * (1-beta17_m)^((1-cond_main)*actsm11_main))
D2P3_m12_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm12_main) * (1-beta17_m)^((1-cond_main)*actsm12_main))
D2P3_m13_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm13_main) * (1-beta17_m)^((1-cond_main)*actsm13_main))
D2P3_m21_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm21_main) * (1-beta17_m)^((1-cond_main)*actsm21_main))
D2P3_m22_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm22_main) * (1-beta17_m)^((1-cond_main)*actsm22_main))
D2P3_m23_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm23_main) * (1-beta17_m)^((1-cond_main)*actsm23_main))
D2P3_m31_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm31_main) * (1-beta17_m)^((1-cond_main)*actsm31_main))
D2P3_m32_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm32_main) * (1-beta17_m)^((1-cond_main)*actsm32_main))
D2P3_m33_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm33_main) * (1-beta17_m)^((1-cond_main)*actsm33_main))
D2P3_f22_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta17_f)^((1-cond_reg)*actsf22_reg))
D2P3_f23_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta17_f)^((1-cond_reg)*actsf23_reg))
D2P3_f32_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta17_f)^((1-cond_reg)*actsf32_reg))
D2P3_f33_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta17_f)^((1-cond_reg)*actsf33_reg))
D2P3_m22_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta17_m)^((1-cond_reg)*actsm22_reg))
D2P3_m23_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta17_m)^((1-cond_reg)*actsm23_reg))
D2P3_m32_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta17_m)^((1-cond_reg)*actsm32_reg))
D2P3_m33_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta17_m)^((1-cond_reg)*actsm33_reg))
D2P3_f33_casual = 1 - ((1-beta17_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta17_f)^((1-cond_casual)*actsf33_casual))
D2P3_m33_casual = 1 - ((1-beta17_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta17_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in preAIDS HV stage, not on HCT/ART

D2X1_f11_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf11_main) * (1-beta10_f)^((1-cond_main)*actsf11_main))
D2X1_f12_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf12_main) * (1-beta10_f)^((1-cond_main)*actsf12_main))
D2X1_f13_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf13_main) * (1-beta10_f)^((1-cond_main)*actsf13_main))
D2X1_f21_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf21_main) * (1-beta10_f)^((1-cond_main)*actsf21_main))
D2X1_f22_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf22_main) * (1-beta10_f)^((1-cond_main)*actsf22_main))
D2X1_f23_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf23_main) * (1-beta10_f)^((1-cond_main)*actsf23_main))
D2X1_f31_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf31_main) * (1-beta10_f)^((1-cond_main)*actsf31_main))
D2X1_f32_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf32_main) * (1-beta10_f)^((1-cond_main)*actsf32_main))
D2X1_f33_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf33_main) * (1-beta10_f)^((1-cond_main)*actsf33_main))
D2X1_m11_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm11_main) * (1-beta10_m)^((1-cond_main)*actsm11_main))
D2X1_m12_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm12_main) * (1-beta10_m)^((1-cond_main)*actsm12_main))
D2X1_m13_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm13_main) * (1-beta10_m)^((1-cond_main)*actsm13_main))
D2X1_m21_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm21_main) * (1-beta10_m)^((1-cond_main)*actsm21_main))
D2X1_m22_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm22_main) * (1-beta10_m)^((1-cond_main)*actsm22_main))
D2X1_m23_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm23_main) * (1-beta10_m)^((1-cond_main)*actsm23_main))
D2X1_m31_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm31_main) * (1-beta10_m)^((1-cond_main)*actsm31_main))
D2X1_m32_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm32_main) * (1-beta10_m)^((1-cond_main)*actsm32_main))
D2X1_m33_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm33_main) * (1-beta10_m)^((1-cond_main)*actsm33_main))
D2X1_f22_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta10_f)^((1-cond_reg)*actsf22_reg))
D2X1_f23_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta10_f)^((1-cond_reg)*actsf23_reg))
D2X1_f32_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta10_f)^((1-cond_reg)*actsf32_reg))
D2X1_f33_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta10_f)^((1-cond_reg)*actsf33_reg))
D2X1_m22_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta10_m)^((1-cond_reg)*actsm22_reg))
D2X1_m23_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta10_m)^((1-cond_reg)*actsm23_reg))
D2X1_m32_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta10_m)^((1-cond_reg)*actsm32_reg))
D2X1_m33_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta10_m)^((1-cond_reg)*actsm33_reg))
D2X1_f33_casual = 1 - ((1-beta10_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta10_f)^((1-cond_casual)*actsf33_casual))
D2X1_m33_casual = 1 - ((1-beta10_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta10_m)^((1-cond_casual)*actsm33_casual))

## Partner in preAIDS HV stage, on HCT

D2X2_f11_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf11_main))
D2X2_f12_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf12_main))
D2X2_f13_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf13_main))
D2X2_f21_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf21_main))
D2X2_f22_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf22_main))
D2X2_f23_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf23_main))
D2X2_f31_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf31_main))
D2X2_f32_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf32_main))
D2X2_f33_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf33_main))
D2X2_m11_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
D2X2_m12_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
D2X2_m13_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
D2X2_m21_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
D2X2_m22_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
D2X2_m23_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
D2X2_m31_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
D2X2_m32_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
D2X2_m33_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
D2X2_f22_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf22_reg))
D2X2_f23_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf23_reg))
D2X2_f32_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf32_reg))
D2X2_f33_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf33_reg))
D2X2_m22_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
D2X2_m23_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
D2X2_m32_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
D2X2_m33_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
D2X2_f33_casual = 1 - ((1-beta10_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta10_f)^((1-cond_casual)*HCT_cond*HCT_red*actsf33_casual))
D2X2_m33_casual = 1 - ((1-beta10_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta10_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))

## Partner in preAIDS HV stage, on ART

D2X3_f11_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf11_main) * (1-beta18_f)^((1-cond_main)*actsf11_main))
D2X3_f12_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf12_main) * (1-beta18_f)^((1-cond_main)*actsf12_main))
D2X3_f13_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf13_main) * (1-beta18_f)^((1-cond_main)*actsf13_main))
D2X3_f21_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf21_main) * (1-beta18_f)^((1-cond_main)*actsf21_main))
D2X3_f22_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf22_main) * (1-beta18_f)^((1-cond_main)*actsf22_main))
D2X3_f23_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf23_main) * (1-beta18_f)^((1-cond_main)*actsf23_main))
D2X3_f31_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf31_main) * (1-beta18_f)^((1-cond_main)*actsf31_main))
D2X3_f32_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf32_main) * (1-beta18_f)^((1-cond_main)*actsf32_main))
D2X3_f33_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf33_main) * (1-beta18_f)^((1-cond_main)*actsf33_main))
D2X3_m11_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm11_main) * (1-beta18_m)^((1-cond_main)*actsm11_main))
D2X3_m12_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm12_main) * (1-beta18_m)^((1-cond_main)*actsm12_main))
D2X3_m13_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm13_main) * (1-beta18_m)^((1-cond_main)*actsm13_main))
D2X3_m21_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm21_main) * (1-beta18_m)^((1-cond_main)*actsm21_main))
D2X3_m22_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm22_main) * (1-beta18_m)^((1-cond_main)*actsm22_main))
D2X3_m23_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm23_main) * (1-beta18_m)^((1-cond_main)*actsm23_main))
D2X3_m31_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm31_main) * (1-beta18_m)^((1-cond_main)*actsm31_main))
D2X3_m32_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm32_main) * (1-beta18_m)^((1-cond_main)*actsm32_main))
D2X3_m33_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm33_main) * (1-beta18_m)^((1-cond_main)*actsm33_main))
D2X3_f22_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta18_f)^((1-cond_reg)*actsf22_reg))
D2X3_f23_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta18_f)^((1-cond_reg)*actsf23_reg))
D2X3_f32_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta18_f)^((1-cond_reg)*actsf32_reg))
D2X3_f33_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta18_f)^((1-cond_reg)*actsf33_reg))
D2X3_m22_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta18_m)^((1-cond_reg)*actsm22_reg))
D2X3_m23_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta18_m)^((1-cond_reg)*actsm23_reg))
D2X3_m32_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta18_m)^((1-cond_reg)*actsm32_reg))
D2X3_m33_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta18_m)^((1-cond_reg)*actsm33_reg))
D2X3_f33_casual = 1 - ((1-beta18_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta18_f)^((1-cond_casual)*actsf33_casual))
D2X3_m33_casual = 1 - ((1-beta18_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta18_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in AIDS stage, not on HCT/ART

D2K1_f11_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf11_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf11_main))
D2K1_f12_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf12_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf12_main))
D2K1_f13_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf13_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf13_main))
D2K1_f21_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf21_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf21_main))
D2K1_f22_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf22_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf22_main))
D2K1_f23_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf23_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf23_main))
D2K1_f31_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf31_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf31_main))
D2K1_f32_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf32_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf32_main))
D2K1_f33_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf33_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf33_main))
D2K1_m11_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm11_main))
D2K1_m12_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm12_main))
D2K1_m13_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm13_main))
D2K1_m21_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm21_main))
D2K1_m22_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm22_main))
D2K1_m23_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm23_main))
D2K1_m31_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm31_main))
D2K1_m32_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm32_main))
D2K1_m33_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm33_main))
D2K1_f22_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf22_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf22_reg))
D2K1_f23_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf23_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf23_reg))
D2K1_f32_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf32_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf32_reg))
D2K1_f33_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf33_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf33_reg))
D2K1_m22_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
D2K1_m23_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
D2K1_m32_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
D2K1_m33_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
D2K1_f33_casual = 1 - ((1-beta11_f*cond_eff)^(cond_casual*AIDS_red*actsf33_casual) * (1-beta11_f)^((1-cond_casual)*AIDS_red*actsf33_casual))
D2K1_m33_casual = 1 - ((1-beta11_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta11_m)^((1-cond_casual)*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on HCT

D2K2_f11_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf11_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf11_main))
D2K2_f12_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf12_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf12_main))
D2K2_f13_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf13_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf13_main))
D2K2_f21_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf21_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf21_main))
D2K2_f22_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf22_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf22_main))
D2K2_f23_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf23_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf23_main))
D2K2_f31_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf31_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf31_main))
D2K2_f32_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf32_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf32_main))
D2K2_f33_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf33_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf33_main))
D2K2_m11_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm11_main))
D2K2_m12_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm12_main))
D2K2_m13_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm13_main))
D2K2_m21_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm21_main))
D2K2_m22_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm22_main))
D2K2_m23_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm23_main))
D2K2_m31_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm31_main))
D2K2_m32_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm32_main))
D2K2_m33_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm33_main))
D2K2_f22_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf22_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf22_reg))
D2K2_f23_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf23_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf23_reg))
D2K2_f32_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf32_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf32_reg))
D2K2_f33_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf33_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf33_reg))
D2K2_m22_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm22_reg))
D2K2_m23_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm23_reg))
D2K2_m32_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm32_reg))
D2K2_m33_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm33_reg))
D2K2_f33_casual = 1 - ((1-beta11_f*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsf33_casual) * (1-beta11_f)^((1-cond_casual)*HCT_cond*AIDS_red*actsf33_casual))
D2K2_m33_casual = 1 - ((1-beta11_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta11_m)^((1-cond_casual)*HCT_cond*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on ART

D2K3_f11_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf11_main) * (1-beta19_f)^((1-cond_main)*actsf11_main))
D2K3_f12_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf12_main) * (1-beta19_f)^((1-cond_main)*actsf12_main))
D2K3_f13_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf13_main) * (1-beta19_f)^((1-cond_main)*actsf13_main))
D2K3_f21_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf21_main) * (1-beta19_f)^((1-cond_main)*actsf21_main))
D2K3_f22_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf22_main) * (1-beta19_f)^((1-cond_main)*actsf22_main))
D2K3_f23_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf23_main) * (1-beta19_f)^((1-cond_main)*actsf23_main))
D2K3_f31_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf31_main) * (1-beta19_f)^((1-cond_main)*actsf31_main))
D2K3_f32_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf32_main) * (1-beta19_f)^((1-cond_main)*actsf32_main))
D2K3_f33_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf33_main) * (1-beta19_f)^((1-cond_main)*actsf33_main))
D2K3_m11_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm11_main) * (1-beta19_m)^((1-cond_main)*actsm11_main))
D2K3_m12_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm12_main) * (1-beta19_m)^((1-cond_main)*actsm12_main))
D2K3_m13_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm13_main) * (1-beta19_m)^((1-cond_main)*actsm13_main))
D2K3_m21_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm21_main) * (1-beta19_m)^((1-cond_main)*actsm21_main))
D2K3_m22_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm22_main) * (1-beta19_m)^((1-cond_main)*actsm22_main))
D2K3_m23_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm23_main) * (1-beta19_m)^((1-cond_main)*actsm23_main))
D2K3_m31_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm31_main) * (1-beta19_m)^((1-cond_main)*actsm31_main))
D2K3_m32_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm32_main) * (1-beta19_m)^((1-cond_main)*actsm32_main))
D2K3_m33_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm33_main) * (1-beta19_m)^((1-cond_main)*actsm33_main))
D2K3_f22_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta19_f)^((1-cond_reg)*actsf22_reg))
D2K3_f23_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta19_f)^((1-cond_reg)*actsf23_reg))
D2K3_f32_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta19_f)^((1-cond_reg)*actsf32_reg))
D2K3_f33_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta19_f)^((1-cond_reg)*actsf33_reg))
D2K3_m22_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta19_m)^((1-cond_reg)*actsm22_reg))
D2K3_m23_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta19_m)^((1-cond_reg)*actsm23_reg))
D2K3_m32_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta19_m)^((1-cond_reg)*actsm32_reg))
D2K3_m33_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta19_m)^((1-cond_reg)*actsm33_reg))
D2K3_f33_casual = 1 - ((1-beta19_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta19_f)^((1-cond_casual)*actsf33_casual))
D2K3_m33_casual = 1 - ((1-beta19_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta19_m)^((1-cond_casual)*actsm33_casual))

### Circumcised males not on PrEP ###

## Partner in LV stage, not on HCT/ART

D3G1_m11_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm11_main) * (1-beta12_m)^((1-cond_main)*actsm11_main))
D3G1_m12_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm12_main) * (1-beta12_m)^((1-cond_main)*actsm12_main))
D3G1_m13_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm13_main) * (1-beta12_m)^((1-cond_main)*actsm13_main))
D3G1_m21_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm21_main) * (1-beta12_m)^((1-cond_main)*actsm21_main))
D3G1_m22_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm22_main) * (1-beta12_m)^((1-cond_main)*actsm22_main))
D3G1_m23_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm23_main) * (1-beta12_m)^((1-cond_main)*actsm23_main))
D3G1_m31_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm31_main) * (1-beta12_m)^((1-cond_main)*actsm31_main))
D3G1_m32_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm32_main) * (1-beta12_m)^((1-cond_main)*actsm32_main))
D3G1_m33_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm33_main) * (1-beta12_m)^((1-cond_main)*actsm33_main))
D3G1_m22_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta12_m)^((1-cond_reg)*actsm22_reg))
D3G1_m23_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta12_m)^((1-cond_reg)*actsm23_reg))
D3G1_m32_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta12_m)^((1-cond_reg)*actsm32_reg))
D3G1_m33_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta12_m)^((1-cond_reg)*actsm33_reg))
D3G1_m33_casual = 1 - ((1-beta12_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta12_m)^((1-cond_casual)*actsm33_casual))

## Partner in LV stage, on HCT

D3G2_m11_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
D3G2_m12_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
D3G2_m13_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
D3G2_m21_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
D3G2_m22_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
D3G2_m23_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
D3G2_m31_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
D3G2_m32_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
D3G2_m33_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
D3G2_m22_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
D3G2_m23_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
D3G2_m32_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
D3G2_m33_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
D3G2_m33_casual = 1 - ((1-beta12_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta12_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))

## Partner in LV stage, on ART

D3G3_m11_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm11_main) * (1-beta20_m)^((1-cond_main)*actsm11_main))
D3G3_m12_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm12_main) * (1-beta20_m)^((1-cond_main)*actsm12_main))
D3G3_m13_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm13_main) * (1-beta20_m)^((1-cond_main)*actsm13_main))
D3G3_m21_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm21_main) * (1-beta20_m)^((1-cond_main)*actsm21_main))
D3G3_m22_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm22_main) * (1-beta20_m)^((1-cond_main)*actsm22_main))
D3G3_m23_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm23_main) * (1-beta20_m)^((1-cond_main)*actsm23_main))
D3G3_m31_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm31_main) * (1-beta20_m)^((1-cond_main)*actsm31_main))
D3G3_m32_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm32_main) * (1-beta20_m)^((1-cond_main)*actsm32_main))
D3G3_m33_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm33_main) * (1-beta20_m)^((1-cond_main)*actsm33_main))
D3G3_m22_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta20_m)^((1-cond_reg)*actsm22_reg))
D3G3_m23_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta20_m)^((1-cond_reg)*actsm23_reg))
D3G3_m32_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta20_m)^((1-cond_reg)*actsm32_reg))
D3G3_m33_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta20_m)^((1-cond_reg)*actsm33_reg))
D3G3_m33_casual = 1 - ((1-beta20_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta20_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in initial HV stage, not on HCT/ART

D3P1_m11_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm11_main) * (1-beta13_m)^((1-cond_main)*actsm11_main))
D3P1_m12_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm12_main) * (1-beta13_m)^((1-cond_main)*actsm12_main))
D3P1_m13_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm13_main) * (1-beta13_m)^((1-cond_main)*actsm13_main))
D3P1_m21_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm21_main) * (1-beta13_m)^((1-cond_main)*actsm21_main))
D3P1_m22_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm22_main) * (1-beta13_m)^((1-cond_main)*actsm22_main))
D3P1_m23_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm23_main) * (1-beta13_m)^((1-cond_main)*actsm23_main))
D3P1_m31_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm31_main) * (1-beta13_m)^((1-cond_main)*actsm31_main))
D3P1_m32_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm32_main) * (1-beta13_m)^((1-cond_main)*actsm32_main))
D3P1_m33_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm33_main) * (1-beta13_m)^((1-cond_main)*actsm33_main))
D3P1_m22_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta13_m)^((1-cond_reg)*actsm22_reg))
D3P1_m23_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta13_m)^((1-cond_reg)*actsm23_reg))
D3P1_m32_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta13_m)^((1-cond_reg)*actsm32_reg))
D3P1_m33_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta13_m)^((1-cond_reg)*actsm33_reg))
D3P1_m33_casual = 1 - ((1-beta13_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta13_m)^((1-cond_casual)*actsm33_casual))

## Partner in initial HV stage, on HCT

D3P2_m11_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
D3P2_m12_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
D3P2_m13_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
D3P2_m21_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
D3P2_m22_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
D3P2_m23_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
D3P2_m31_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
D3P2_m32_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
D3P2_m33_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
D3P2_m22_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
D3P2_m23_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
D3P2_m32_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
D3P2_m33_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
D3P2_m33_casual = 1 - ((1-beta13_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta13_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))

## Partner in initial HV stage, on ART

D3P3_m11_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm11_main) * (1-beta21_m)^((1-cond_main)*actsm11_main))
D3P3_m12_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm12_main) * (1-beta21_m)^((1-cond_main)*actsm12_main))
D3P3_m13_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm13_main) * (1-beta21_m)^((1-cond_main)*actsm13_main))
D3P3_m21_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm21_main) * (1-beta21_m)^((1-cond_main)*actsm21_main))
D3P3_m22_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm22_main) * (1-beta21_m)^((1-cond_main)*actsm22_main))
D3P3_m23_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm23_main) * (1-beta21_m)^((1-cond_main)*actsm23_main))
D3P3_m31_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm31_main) * (1-beta21_m)^((1-cond_main)*actsm31_main))
D3P3_m32_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm32_main) * (1-beta21_m)^((1-cond_main)*actsm32_main))
D3P3_m33_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm33_main) * (1-beta21_m)^((1-cond_main)*actsm33_main))
D3P3_m22_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta21_m)^((1-cond_reg)*actsm22_reg))
D3P3_m23_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta21_m)^((1-cond_reg)*actsm23_reg))
D3P3_m32_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta21_m)^((1-cond_reg)*actsm32_reg))
D3P3_m33_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta21_m)^((1-cond_reg)*actsm33_reg))
D3P3_m33_casual = 1 - ((1-beta21_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta21_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in preAIDS HV stage, not on HCT/ART

D3X1_m11_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm11_main) * (1-beta14_m)^((1-cond_main)*actsm11_main))
D3X1_m12_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm12_main) * (1-beta14_m)^((1-cond_main)*actsm12_main))
D3X1_m13_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm13_main) * (1-beta14_m)^((1-cond_main)*actsm13_main))
D3X1_m21_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm21_main) * (1-beta14_m)^((1-cond_main)*actsm21_main))
D3X1_m22_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm22_main) * (1-beta14_m)^((1-cond_main)*actsm22_main))
D3X1_m23_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm23_main) * (1-beta14_m)^((1-cond_main)*actsm23_main))
D3X1_m31_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm31_main) * (1-beta14_m)^((1-cond_main)*actsm31_main))
D3X1_m32_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm32_main) * (1-beta14_m)^((1-cond_main)*actsm32_main))
D3X1_m33_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm33_main) * (1-beta14_m)^((1-cond_main)*actsm33_main))
D3X1_m22_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta14_m)^((1-cond_reg)*actsm22_reg))
D3X1_m23_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta14_m)^((1-cond_reg)*actsm23_reg))
D3X1_m32_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta14_m)^((1-cond_reg)*actsm32_reg))
D3X1_m33_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta14_m)^((1-cond_reg)*actsm33_reg))
D3X1_m33_casual = 1 - ((1-beta14_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta14_m)^((1-cond_casual)*actsm33_casual))

## Partner in preAIDS HV stage, on HCT

D3X2_m11_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
D3X2_m12_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
D3X2_m13_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
D3X2_m21_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
D3X2_m22_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
D3X2_m23_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
D3X2_m31_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
D3X2_m32_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
D3X2_m33_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
D3X2_m22_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
D3X2_m23_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
D3X2_m32_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
D3X2_m33_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
D3X2_m33_casual = 1 - ((1-beta14_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta14_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))

## Partner in preAIDS HV stage, on ART

D3X3_m11_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm11_main) * (1-beta22_m)^((1-cond_main)*actsm11_main))
D3X3_m12_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm12_main) * (1-beta22_m)^((1-cond_main)*actsm12_main))
D3X3_m13_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm13_main) * (1-beta22_m)^((1-cond_main)*actsm13_main))
D3X3_m21_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm21_main) * (1-beta22_m)^((1-cond_main)*actsm21_main))
D3X3_m22_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm22_main) * (1-beta22_m)^((1-cond_main)*actsm22_main))
D3X3_m23_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm23_main) * (1-beta22_m)^((1-cond_main)*actsm23_main))
D3X3_m31_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm31_main) * (1-beta22_m)^((1-cond_main)*actsm31_main))
D3X3_m32_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm32_main) * (1-beta22_m)^((1-cond_main)*actsm32_main))
D3X3_m33_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm33_main) * (1-beta22_m)^((1-cond_main)*actsm33_main))
D3X3_m22_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta22_m)^((1-cond_reg)*actsm22_reg))
D3X3_m23_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta22_m)^((1-cond_reg)*actsm23_reg))
D3X3_m32_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta22_m)^((1-cond_reg)*actsm32_reg))
D3X3_m33_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta22_m)^((1-cond_reg)*actsm33_reg))
D3X3_m33_casual = 1 - ((1-beta22_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta22_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in AIDS stage, not on HCT/ART

D3K1_m11_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm11_main))
D3K1_m12_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm12_main))
D3K1_m13_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm13_main))
D3K1_m21_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm21_main))
D3K1_m22_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm22_main))
D3K1_m23_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm23_main))
D3K1_m31_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm31_main))
D3K1_m32_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm32_main))
D3K1_m33_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm33_main))
D3K1_m22_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
D3K1_m23_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
D3K1_m32_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
D3K1_m33_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
D3K1_m33_casual = 1 - ((1-beta15_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta15_m)^((1-cond_casual)*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on HCT

D3K2_m11_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm11_main))
D3K2_m12_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm12_main))
D3K2_m13_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm13_main))
D3K2_m21_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm21_main))
D3K2_m22_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm22_main))
D3K2_m23_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm23_main))
D3K2_m31_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm31_main))
D3K2_m32_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm32_main))
D3K2_m33_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm33_main))
D3K2_m22_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm22_reg))
D3K2_m23_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm23_reg))
D3K2_m32_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm32_reg))
D3K2_m33_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm33_reg))
D3K2_m33_casual = 1 - ((1-beta15_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta15_m)^((1-cond_casual*HCT_cond)*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on ART

D3K3_m11_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm11_main) * (1-beta23_m)^((1-cond_main)*actsm11_main))
D3K3_m12_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm12_main) * (1-beta23_m)^((1-cond_main)*actsm12_main))
D3K3_m13_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm13_main) * (1-beta23_m)^((1-cond_main)*actsm13_main))
D3K3_m21_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm21_main) * (1-beta23_m)^((1-cond_main)*actsm21_main))
D3K3_m22_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm22_main) * (1-beta23_m)^((1-cond_main)*actsm22_main))
D3K3_m23_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm23_main) * (1-beta23_m)^((1-cond_main)*actsm23_main))
D3K3_m31_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm31_main) * (1-beta23_m)^((1-cond_main)*actsm31_main))
D3K3_m32_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm32_main) * (1-beta23_m)^((1-cond_main)*actsm32_main))
D3K3_m33_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm33_main) * (1-beta23_m)^((1-cond_main)*actsm33_main))
D3K3_m22_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta23_m)^((1-cond_reg)*actsm22_reg))
D3K3_m23_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta23_m)^((1-cond_reg)*actsm23_reg))
D3K3_m32_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta23_m)^((1-cond_reg)*actsm32_reg))
D3K3_m33_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta23_m)^((1-cond_reg)*actsm33_reg))
D3K3_m33_casual = 1 - ((1-beta23_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta23_m)^((1-cond_casual)*actsm33_casual))

### Circumcised males on PrEP ###

## Partner in LV stage, not on HCT/ART

D4G1_m11_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm11_main) * (1-beta28_m)^((1-cond_main)*actsm11_main))
D4G1_m12_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm12_main) * (1-beta28_m)^((1-cond_main)*actsm12_main))
D4G1_m13_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm13_main) * (1-beta28_m)^((1-cond_main)*actsm13_main))
D4G1_m21_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm21_main) * (1-beta28_m)^((1-cond_main)*actsm21_main))
D4G1_m22_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm22_main) * (1-beta28_m)^((1-cond_main)*actsm22_main))
D4G1_m23_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm23_main) * (1-beta28_m)^((1-cond_main)*actsm23_main))
D4G1_m31_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm31_main) * (1-beta28_m)^((1-cond_main)*actsm31_main))
D4G1_m32_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm32_main) * (1-beta28_m)^((1-cond_main)*actsm32_main))
D4G1_m33_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm33_main) * (1-beta28_m)^((1-cond_main)*actsm33_main))
D4G1_m22_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta28_m)^((1-cond_reg)*actsm22_reg))
D4G1_m23_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta28_m)^((1-cond_reg)*actsm23_reg))
D4G1_m32_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta28_m)^((1-cond_reg)*actsm32_reg))
D4G1_m33_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta28_m)^((1-cond_reg)*actsm33_reg))
D4G1_m33_casual = 1 - ((1-beta28_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta28_m)^((1-cond_casual)*actsm33_casual))

## Partner in LV stage, on HCT

D4G2_m11_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
D4G2_m12_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
D4G2_m13_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
D4G2_m21_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
D4G2_m22_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
D4G2_m23_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
D4G2_m31_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
D4G2_m32_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
D4G2_m33_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
D4G2_m22_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
D4G2_m23_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
D4G2_m32_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
D4G2_m33_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
D4G2_m33_casual = 1 - ((1-beta28_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta28_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))

## Partner in LV stage, on ART

D4G3_m11_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm11_main) * (1-beta24_m)^((1-cond_main)*actsm11_main))
D4G3_m12_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm12_main) * (1-beta24_m)^((1-cond_main)*actsm12_main))
D4G3_m13_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm13_main) * (1-beta24_m)^((1-cond_main)*actsm13_main))
D4G3_m21_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm21_main) * (1-beta24_m)^((1-cond_main)*actsm21_main))
D4G3_m22_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm22_main) * (1-beta24_m)^((1-cond_main)*actsm22_main))
D4G3_m23_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm23_main) * (1-beta24_m)^((1-cond_main)*actsm23_main))
D4G3_m31_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm31_main) * (1-beta24_m)^((1-cond_main)*actsm31_main))
D4G3_m32_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm32_main) * (1-beta24_m)^((1-cond_main)*actsm32_main))
D4G3_m33_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm33_main) * (1-beta24_m)^((1-cond_main)*actsm33_main))
D4G3_m22_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta24_m)^((1-cond_reg)*actsm22_reg))
D4G3_m23_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta24_m)^((1-cond_reg)*actsm23_reg))
D4G3_m32_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta24_m)^((1-cond_reg)*actsm32_reg))
D4G3_m33_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta24_m)^((1-cond_reg)*actsm33_reg))
D4G3_m33_casual = 1 - ((1-beta24_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta24_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in initial HV stage, not on HCT/ART

D4P1_m11_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm11_main) * (1-beta29_m)^((1-cond_main)*actsm11_main))
D4P1_m12_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm12_main) * (1-beta29_m)^((1-cond_main)*actsm12_main))
D4P1_m13_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm13_main) * (1-beta29_m)^((1-cond_main)*actsm13_main))
D4P1_m21_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm21_main) * (1-beta29_m)^((1-cond_main)*actsm21_main))
D4P1_m22_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm22_main) * (1-beta29_m)^((1-cond_main)*actsm22_main))
D4P1_m23_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm23_main) * (1-beta29_m)^((1-cond_main)*actsm23_main))
D4P1_m31_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm31_main) * (1-beta29_m)^((1-cond_main)*actsm31_main))
D4P1_m32_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm32_main) * (1-beta29_m)^((1-cond_main)*actsm32_main))
D4P1_m33_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm33_main) * (1-beta29_m)^((1-cond_main)*actsm33_main))
D4P1_m22_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta29_m)^((1-cond_reg)*actsm22_reg))
D4P1_m23_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta29_m)^((1-cond_reg)*actsm23_reg))
D4P1_m32_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta29_m)^((1-cond_reg)*actsm32_reg))
D4P1_m33_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta29_m)^((1-cond_reg)*actsm33_reg))
D4P1_m33_casual = 1 - ((1-beta29_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta29_m)^((1-cond_casual)*actsm33_casual))

## Partner in initial HV stage, on HCT

D4P2_m11_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
D4P2_m12_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
D4P2_m13_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
D4P2_m21_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
D4P2_m22_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
D4P2_m23_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
D4P2_m31_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
D4P2_m32_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
D4P2_m33_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
D4P2_m22_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
D4P2_m23_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
D4P2_m32_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
D4P2_m33_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
D4P2_m33_casual = 1 - ((1-beta29_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta29_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))

## Partner in initial HV stage, on ART

D4P3_m11_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm11_main) * (1-beta25_m)^((1-cond_main)*actsm11_main))
D4P3_m12_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm12_main) * (1-beta25_m)^((1-cond_main)*actsm12_main))
D4P3_m13_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm13_main) * (1-beta25_m)^((1-cond_main)*actsm13_main))
D4P3_m21_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm21_main) * (1-beta25_m)^((1-cond_main)*actsm21_main))
D4P3_m22_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm22_main) * (1-beta25_m)^((1-cond_main)*actsm22_main))
D4P3_m23_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm23_main) * (1-beta25_m)^((1-cond_main)*actsm23_main))
D4P3_m31_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm31_main) * (1-beta25_m)^((1-cond_main)*actsm31_main))
D4P3_m32_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm32_main) * (1-beta25_m)^((1-cond_main)*actsm32_main))
D4P3_m33_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm33_main) * (1-beta25_m)^((1-cond_main)*actsm33_main))
D4P3_m22_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta25_m)^((1-cond_reg)*actsm22_reg))
D4P3_m23_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta25_m)^((1-cond_reg)*actsm23_reg))
D4P3_m32_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta25_m)^((1-cond_reg)*actsm32_reg))
D4P3_m33_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta25_m)^((1-cond_reg)*actsm33_reg))
D4P3_m33_casual = 1 - ((1-beta25_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta25_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in preAIDS HV stage, not on HCT/ART

D4X1_m11_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm11_main) * (1-beta30_m)^((1-cond_main)*actsm11_main))
D4X1_m12_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm12_main) * (1-beta30_m)^((1-cond_main)*actsm12_main))
D4X1_m13_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm13_main) * (1-beta30_m)^((1-cond_main)*actsm13_main))
D4X1_m21_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm21_main) * (1-beta30_m)^((1-cond_main)*actsm21_main))
D4X1_m22_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm22_main) * (1-beta30_m)^((1-cond_main)*actsm22_main))
D4X1_m23_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm23_main) * (1-beta30_m)^((1-cond_main)*actsm23_main))
D4X1_m31_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm31_main) * (1-beta30_m)^((1-cond_main)*actsm31_main))
D4X1_m32_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm32_main) * (1-beta30_m)^((1-cond_main)*actsm32_main))
D4X1_m33_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm33_main) * (1-beta30_m)^((1-cond_main)*actsm33_main))
D4X1_m22_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta30_m)^((1-cond_reg)*actsm22_reg))
D4X1_m23_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta30_m)^((1-cond_reg)*actsm23_reg))
D4X1_m32_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta30_m)^((1-cond_reg)*actsm32_reg))
D4X1_m33_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta30_m)^((1-cond_reg)*actsm33_reg))
D4X1_m33_casual = 1 - ((1-beta30_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta30_m)^((1-cond_casual)*actsm33_casual))

## Partner in preAIDS HV stage, on HCT

D4X2_m11_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
D4X2_m12_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
D4X2_m13_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
D4X2_m21_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
D4X2_m22_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
D4X2_m23_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
D4X2_m31_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
D4X2_m32_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
D4X2_m33_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
D4X2_m22_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
D4X2_m23_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
D4X2_m32_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
D4X2_m33_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
D4X2_m33_casual = 1 - ((1-beta30_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta30_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))

## Partner in preAIDS HV stage, on ART

D4X3_m11_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm11_main) * (1-beta26_m)^((1-cond_main)*actsm11_main))
D4X3_m12_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm12_main) * (1-beta26_m)^((1-cond_main)*actsm12_main))
D4X3_m13_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm13_main) * (1-beta26_m)^((1-cond_main)*actsm13_main))
D4X3_m21_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm21_main) * (1-beta26_m)^((1-cond_main)*actsm21_main))
D4X3_m22_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm22_main) * (1-beta26_m)^((1-cond_main)*actsm22_main))
D4X3_m23_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm23_main) * (1-beta26_m)^((1-cond_main)*actsm23_main))
D4X3_m31_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm31_main) * (1-beta26_m)^((1-cond_main)*actsm31_main))
D4X3_m32_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm32_main) * (1-beta26_m)^((1-cond_main)*actsm32_main))
D4X3_m33_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm33_main) * (1-beta26_m)^((1-cond_main)*actsm33_main))
D4X3_m22_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta26_m)^((1-cond_reg)*actsm22_reg))
D4X3_m23_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta26_m)^((1-cond_reg)*actsm23_reg))
D4X3_m32_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta26_m)^((1-cond_reg)*actsm32_reg))
D4X3_m33_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta26_m)^((1-cond_reg)*actsm33_reg))
D4X3_m33_casual = 1 - ((1-beta26_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta26_m)^((1-cond_casual)*actsm33_casual))

#####
## Partner in AIDS stage, not on HCT/ART

D4K1_m11_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm11_main))
D4K1_m12_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm12_main))
D4K1_m13_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm13_main))
D4K1_m21_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm21_main))
D4K1_m22_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm22_main))
D4K1_m23_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm23_main))
D4K1_m31_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm31_main))
D4K1_m32_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm32_main))
D4K1_m33_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm33_main))
D4K1_m22_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
D4K1_m23_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
D4K1_m32_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
D4K1_m33_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
D4K1_m33_casual = 1 - ((1-beta31_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta31_m)^((1-cond_casual)*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on HCT

D4K2_m11_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm11_main))
D4K2_m12_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm12_main))
D4K2_m13_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm13_main))
D4K2_m21_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm21_main))
D4K2_m22_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm22_main))
D4K2_m23_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm23_main))
D4K2_m31_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm31_main))
D4K2_m32_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm32_main))
D4K2_m33_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm33_main))
D4K2_m22_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm22_reg))
D4K2_m23_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm23_reg))
D4K2_m32_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm32_reg))
D4K2_m33_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm33_reg))
D4K2_m33_casual = 1 - ((1-beta31_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta31_m)^((1-cond_casual)*HCT_cond*AIDS_red*actsm33_casual))

## Partner in AIDS stage, on ART

D4K3_m11_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm11_main) * (1-beta27_m)^((1-cond_main)*actsm11_main))
D4K3_m12_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm12_main) * (1-beta27_m)^((1-cond_main)*actsm12_main))
D4K3_m13_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm13_main) * (1-beta27_m)^((1-cond_main)*actsm13_main))
D4K3_m21_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm21_main) * (1-beta27_m)^((1-cond_main)*actsm21_main))
D4K3_m22_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm22_main) * (1-beta27_m)^((1-cond_main)*actsm22_main))
D4K3_m23_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm23_main) * (1-beta27_m)^((1-cond_main)*actsm23_main))
D4K3_m31_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm31_main) * (1-beta27_m)^((1-cond_main)*actsm31_main))
D4K3_m32_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm32_main) * (1-beta27_m)^((1-cond_main)*actsm32_main))
D4K3_m33_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm33_main) * (1-beta27_m)^((1-cond_main)*actsm33_main))
D4K3_m22_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta27_m)^((1-cond_reg)*actsm22_reg))
D4K3_m23_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta27_m)^((1-cond_reg)*actsm23_reg))
D4K3_m32_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta27_m)^((1-cond_reg)*actsm32_reg))
D4K3_m33_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta27_m)^((1-cond_reg)*actsm33_reg))
D4K3_m33_casual = 1 - ((1-beta27_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta27_m)^((1-cond_casual)*actsm33_casual))

##########################
### FORCE OF INFECTION ###

# Females not on PrEP
pi1_f1 = (
  # low-low
  (df11_main*rho_f1_1_main)*(
    # infected partner not on HCT/ART
      D1P1_f11_main*(H1m1/Nm1) + D1G1_f11_main*(Y1m1/Nm1) + D1X1_f11_main*(Z1m1/Nm1) + D1K1_f11_main*(A1m1/Nm1) # = (D1P1_f11_main*H1m1 + D1G1_f11_main*Y1m1 + D1X1_f11_main*Z1m1 + D1K1_f11_main*A1m1)/Nm1
    # infected partner on HCT
    + D1P2_f11_main*(H2m1/Nm1) + D1G2_f11_main*(Y2m1/Nm1) + D1X2_f11_main*(Z2m1/Nm1) + D1K2_f11_main*(A2m1/Nm1)
    # infected partner on ART
    + D1P3_f11_main*(H3m1/Nm1) + D1G3_f11_main*(Y3m1/Nm1) + D1X3_f11_main*(Z3m1/Nm1) + D1K3_f11_main*(A3m1/Nm1)
  )
  # low-medium
  +(df12_main*rho_f1_2_main)*(
    # infected partner not on HCT/ART
      D1P1_f12_main*(H1m2/Nm2) + D1G1_f12_main*(Y1m2/Nm2) + D1X1_f12_main*(Z1m2/Nm2) + D1K1_f12_main*(A1m2/Nm2)
    # infected partner on HCT
    + D1P2_f12_main*(H2m2/Nm2) + D1G2_f12_main*(Y2m2/Nm2) + D1X2_f12_main*(Z2m2/Nm2) + D1K2_f12_main*(A2m2/Nm2)
    # infected partner on ART
    + D1P3_f12_main*(H3m2/Nm2) + D1G3_f12_main*(Y3m2/Nm2) + D1X3_f12_main*(Z3m2/Nm2) + D1K3_f12_main*(A3m2/Nm2)
  )
  # low-high
  +(df13_main*rho_f1_3_main)*(
    # infected partner not on HCT/ART
      D1P1_f13_main*(H1m3/Nm3) + D1G1_f13_main*(Y1m3/Nm3) + D1X1_f13_main*(Z1m3/Nm3) + D1K1_f13_main*(A1m3/Nm3)
    # infected partner on HCT
    + D1P2_f13_main*(H2m3/Nm3) + D1G2_f13_main*(Y2m3/Nm3) + D1X2_f13_main*(Z2m3/Nm3) + D1K2_f13_main*(A2m3/Nm3)
    # infected partner on ART
    + D1P3_f13_main*(H3m3/Nm3) + D1G3_f13_main*(Y3m3/Nm3) + D1X3_f13_main*(Z3m3/Nm3) + D1K3_f13_main*(A3m3/Nm3)
  )
)


pi1_f2 = (
  # main
  (df21_main*rho_f2_1_main)*(
    # infected partner not on HCT/ART
      D1P1_f21_main*(H1m1/Nm1) + D1G1_f21_main*(Y1m1/Nm1) + D1X1_f21_main*(Z1m1/Nm1) + D1K1_f21_main*(A1m1/Nm1)
    # infected partner on HCT
    + D1P2_f21_main*(H2m1/Nm1) + D1G2_f21_main*(Y2m1/Nm1) + D1X2_f21_main*(Z2m1/Nm1) + D1K2_f21_main*(A2m1/Nm1)
    # infected partner on ART
    + D1P3_f21_main*(H3m1/Nm1) + D1G3_f21_main*(Y3m1/Nm1) + D1X3_f21_main*(Z3m1/Nm1) + D1K3_f21_main*(A3m1/Nm1)
  )
  +(df22_main*rho_f2_2_main)*(
    # infected partner not on HCT/ART
      D1P1_f22_main*(H1m2/Nm2) + D1G1_f22_main*(Y1m2/Nm2) + D1X1_f22_main*(Z1m2/Nm2) + D1K1_f22_main*(A1m2/Nm2)
    # infected partner on HCT
    + D1P2_f22_main*(H2m2/Nm2) + D1G2_f22_main*(Y2m2/Nm2) + D1X2_f22_main*(Z2m2/Nm2) + D1K2_f22_main*(A2m2/Nm2)
    # infected partner on ART
    + D1P3_f22_main*(H3m2/Nm2) + D1G3_f22_main*(Y3m2/Nm2) + D1X3_f22_main*(Z3m2/Nm2) + D1K3_f22_main*(A3m2/Nm2)
  )
  +(df23_main*rho_f2_3_main)*(
    # infected partner not on HCT/ART
      D1P1_f23_main*(H1m3/Nm3) + D1G1_f23_main*(Y1m3/Nm3) + D1X1_f23_main*(Z1m3/Nm3) + D1K1_f23_main*(A1m3/Nm3)
    # infected partner on HCT
    + D1P2_f23_main*(H2m3/Nm3) + D1G2_f23_main*(Y2m3/Nm3) + D1X2_f23_main*(Z2m3/Nm3) + D1K2_f23_main*(A2m3/Nm3)
    # infected partner on ART
    + D1P3_f23_main*(H3m3/Nm3) + D1G3_f23_main*(Y3m3/Nm3) + D1X3_f23_main*(Z3m3/Nm3) + D1K3_f23_main*(A3m3/Nm3)
  )
  # regular
  +(df22_reg*rho_f2_2_reg)*(
    # infected partner not on HCT/ART
      D1P1_f22_reg*(H1m2/Nm2) + D1G1_f22_reg*(Y1m2/Nm2) + D1X1_f22_reg*(Z1m2/Nm2) + D1K1_f22_reg*(A1m2/Nm2)
    # infected partner on HCT
    + D1P2_f22_reg*(H2m2/Nm2) + D1G2_f22_reg*(Y2m2/Nm2) + D1X2_f22_reg*(Z2m2/Nm2) + D1K2_f22_reg*(A2m2/Nm2)
    # infected partner on ART
    + D1P3_f22_reg*(H3m2/Nm2) + D1G3_f22_reg*(Y3m2/Nm2) + D1X3_f22_reg*(Z3m2/Nm2) + D1K3_f22_reg*(A3m2/Nm2)
  )
  +(df23_reg*rho_f2_3_reg)*(
    # infected partner not on HCT/ART
      D1P1_f23_reg*(H1m3/Nm3) + D1G1_f23_reg*(Y1m3/Nm3) + D1X1_f23_reg*(Z1m3/Nm3) + D1K1_f23_reg*(A1m3/Nm3)
    # infected partner on HCT
    + D1P2_f23_reg*(H2m3/Nm3) + D1G2_f23_reg*(Y2m3/Nm3) + D1X2_f23_reg*(Z2m3/Nm3) + D1K2_f23_reg*(A2m3/Nm3)
    # infected partner on ART
    + D1P3_f23_reg*(H3m3/Nm3) + D1G3_f23_reg*(Y3m3/Nm3) + D1X3_f23_reg*(Z3m3/Nm3) + D1K3_f23_reg*(A3m3/Nm3)
  )
)

# pitest = (
#   (1 - (((D1P1_f21_main*(H1m1/Nm1))^(df21_main*rho_f2_1_main)) * ((D1G1_f21_main*(Y1m1/Nm1))^(df21_main*rho_f2_1_main)) * ((D1X1_f21_main*(Z1m1/Nm1))^(df21_main*rho_f2_1_main)) * ((D1K1_f11_main*(A1m1/Nm1))^(df21_main*rho_f2_1_main))))
# )

pi1_f3 = (
  # main
  (df31_main*rho_f3_1_main)*(
    # infected partner not on HCT/ART
      D1P1_f31_main*(H1m1/Nm1) + D1G1_f31_main*(Y1m1/Nm1) + D1X1_f31_main*(Z1m1/Nm1) + D1K1_f31_main*(A1m1/Nm1)
    # infected partner on HCT
    + D1P2_f31_main*(H2m1/Nm1) + D1G2_f31_main*(Y2m1/Nm1) + D1X2_f31_main*(Z2m1/Nm1) + D1K2_f31_main*(A2m1/Nm1)
    # infected partner on ART
    + D1P3_f31_main*(H3m1/Nm1) + D1G3_f31_main*(Y3m1/Nm1) + D1X3_f31_main*(Z3m1/Nm1) + D1K3_f31_main*(A3m1/Nm1)
  )
  +(df32_main*rho_f3_2_main)*(
    # infected partner not on HCT/ART
      D1P1_f32_main*(H1m2/Nm2) + D1G1_f32_main*(Y1m2/Nm2) + D1X1_f32_main*(Z1m2/Nm2) + D1K1_f32_main*(A1m2/Nm2)
    # infected partner on HCT
    + D1P2_f32_main*(H2m2/Nm2) + D1G2_f32_main*(Y2m2/Nm2) + D1X2_f32_main*(Z2m2/Nm2) + D1K2_f32_main*(A2m2/Nm2)
    # infected partner on ART
    + D1P3_f32_main*(H3m2/Nm2) + D1G3_f32_main*(Y3m2/Nm2) + D1X3_f32_main*(Z3m2/Nm2) + D1K3_f32_main*(A3m2/Nm2)
  )
  +(df33_main*rho_f3_3_main)*(
    # infected partner not on HCT/ART
      D1P1_f33_main*(H1m3/Nm3) + D1G1_f33_main*(Y1m3/Nm3) + D1X1_f33_main*(Z1m3/Nm3) + D1K1_f33_main*(A1m3/Nm3)
    # infected partner on HCT
    + D1P2_f33_main*(H2m3/Nm3) + D1G2_f33_main*(Y2m3/Nm3) + D1X2_f33_main*(Z2m3/Nm3) + D1K2_f33_main*(A2m3/Nm3)
    # infected partner on ART
    + D1P3_f33_main*(H3m3/Nm3) + D1G3_f33_main*(Y3m3/Nm3) + D1X3_f33_main*(Z3m3/Nm3) + D1K3_f33_main*(A3m3/Nm3)
  )
  # regular
  +(df32_reg*rho_f3_2_reg)*(
    # infected partner not on HCT/ART
      D1P1_f32_reg*(H1m2/Nm2) + D1G1_f32_reg*(Y1m2/Nm2) + D1X1_f32_reg*(Z1m2/Nm2) + D1K1_f32_reg*(A1m2/Nm2)
    # infected partner on HCT
    + D1P2_f32_reg*(H2m2/Nm2) + D1G2_f32_reg*(Y2m2/Nm2) + D1X2_f32_reg*(Z2m2/Nm2) + D1K2_f32_reg*(A2m2/Nm2)
    # infected partner on ART
    + D1P3_f32_reg*(H3m2/Nm2) + D1G3_f32_reg*(Y3m2/Nm2) + D1X3_f32_reg*(Z3m2/Nm2) + D1K3_f32_reg*(A3m2/Nm2)
  )
  +(df33_reg*rho_f3_3_reg)*(
    # infected partner not on HCT/ART
      D1P1_f33_reg*(H1m3/Nm3) + D1G1_f33_reg*(Y1m3/Nm3) + D1X1_f33_reg*(Z1m3/Nm3) + D1K1_f33_reg*(A1m3/Nm3)
    # infected partner on HCT
    + D1P2_f33_reg*(H2m3/Nm3) + D1G2_f33_reg*(Y2m3/Nm3) + D1X2_f33_reg*(Z2m3/Nm3) + D1K2_f33_reg*(A2m3/Nm3)
    # infected partner on ART
    + D1P3_f33_reg*(H3m3/Nm3) + D1G3_f33_reg*(Y3m3/Nm3) + D1X3_f33_reg*(Z3m3/Nm3) + D1K3_f33_reg*(A3m3/Nm3)
  )
  # casual
  +(df33_casual*rho_f3_3_casual)*(
    # infected partner not on HCT/ART
      D1P1_f33_casual*(H1m3/Nm3) + D1G1_f33_casual*(Y1m3/Nm3) + D1X1_f33_casual*(Z1m3/Nm3) + D1K1_f33_casual*(A1m3/Nm3)
    # infected partner on HCT
    + D1P2_f33_casual*(H2m3/Nm3) + D1G2_f33_casual*(Y2m3/Nm3) + D1X2_f33_casual*(Z2m3/Nm3) + D1K2_f33_casual*(A2m3/Nm3)
    # infected partner on ART
    + D1P3_f33_casual*(H3m3/Nm3) + D1G3_f33_casual*(Y3m3/Nm3) + D1X3_f33_casual*(Z3m3/Nm3) + D1K3_f33_casual*(A3m3/Nm3)
  )
)



# Uncircumcised not on PrEP
pi1_m1 = (
  # low-low
  (dm11_main*rho_m1_1_main)*(
    # infected partner not on HCT/ART
      D1P1_m11_main*(H1f1/Nf1) + D1G1_m11_main*(Y1f1/Nf1) + D1X1_m11_main*(Z1f1/Nf1) + D1K1_m11_main*(A1f1/Nf1)
    # infected partner on HCT
    + D1P2_m11_main*(H2f1/Nf1) + D1G2_m11_main*(Y2f1/Nf1) + D1X2_m11_main*(Z2f1/Nf1) + D1K2_m11_main*(A2f1/Nf1)
    # infected partner on ART
    + D1P3_m11_main*(H3f1/Nf1) + D1G3_m11_main*(Y3f1/Nf1) + D1X3_m11_main*(Z3f1/Nf1) + D1K3_m11_main*(A3f1/Nf1)
  )
  # low-medium
  +(dm12_main*rho_m1_2_main)*(
    # infected partner not on HCT/ART
      D1P1_m12_main*(H1f2/Nf2) + D1G1_m12_main*(Y1f2/Nf2) + D1X1_m12_main*(Z1f2/Nf2) + D1K1_m12_main*(A1f2/Nf2)
    # infected partner on HCT
    + D1P2_m12_main*(H2f2/Nf2) + D1G2_m12_main*(Y2f2/Nf2) + D1X2_m12_main*(Z2f2/Nf2) + D1K2_m12_main*(A2f2/Nf2)
    # infected partner on ART
    + D1P3_m12_main*(H3f2/Nf2) + D1G3_m12_main*(Y3f2/Nf2) + D1X3_m12_main*(Z3f2/Nf2) + D1K3_m12_main*(A3f2/Nf2)
  )
  # low-high
  +(dm13_main*rho_m1_3_main)*(
    # infected partner not on HCT/ART
      D1P1_m13_main*(H1f3/Nf3) + D1G1_m13_main*(Y1f3/Nf3) + D1X1_m13_main*(Z1f3/Nf3) + D1K1_m13_main*(A1f3/Nf3)
    # infected partner on HCT
    + D1P2_m13_main*(H2f3/Nf3) + D1G2_m13_main*(Y2f3/Nf3) + D1X2_m13_main*(Z2f3/Nf3) + D1K2_m13_main*(A2f3/Nf3)
    # infected partner on ART
    + D1P3_m13_main*(H3f3/Nf3) + D1G3_m13_main*(Y3f3/Nf3) + D1X3_m13_main*(Z3f3/Nf3) + D1K3_m13_main*(A3f3/Nf3)
  )
)

pi1_m2 = (
  # main
  (dm21_main*rho_m2_1_main)*(
    # infected partner not on HCT/ART
      D1P1_m21_main*(H1f1/Nf1) + D1G1_m21_main*(Y1f1/Nf1) + D1X1_m21_main*(Z1f1/Nf1) + D1K1_m21_main*(A1f1/Nf1)
    # infected partner on HCT
    + D1P2_m21_main*(H2f1/Nf1) + D1G2_m21_main*(Y2f1/Nf1) + D1X2_m21_main*(Z2f1/Nf1) + D1K2_m21_main*(A2f1/Nf1)
    # infected partner on ART
    + D1P3_m21_main*(H3f1/Nf1) + D1G3_m21_main*(Y3f1/Nf1) + D1X3_m21_main*(Z3f1/Nf1) + D1K3_m21_main*(A3f1/Nf1)
  )
  +(dm22_main*rho_m2_2_main)*(
    # infected partner not on HCT/ART
      D1P1_m22_main*(H1f2/Nf2) + D1G1_m22_main*(Y1f2/Nf2) + D1X1_m22_main*(Z1f2/Nf2) + D1K1_m22_main*(A1f2/Nf2)
    # infected partner on HCT
    + D1P2_m22_main*(H2f2/Nf2) + D1G2_m22_main*(Y2f2/Nf2) + D1X2_m22_main*(Z2f2/Nf2) + D1K2_m22_main*(A2f2/Nf2)
    # infected partner on ART
    + D1P3_m22_main*(H3f2/Nf2) + D1G3_m22_main*(Y3f2/Nf2) + D1X3_m22_main*(Z3f2/Nf2) + D1K3_m22_main*(A3f2/Nf2)
  )
  +(dm23_main*rho_m2_3_main)*(
    # infected partner not on HCT/ART
      D1P1_m23_main*(H1f3/Nf3) + D1G1_m23_main*(Y1f3/Nf3) + D1X1_m23_main*(Z1f3/Nf3) + D1K1_m23_main*(A1f3/Nf3)
    # infected partner on HCT
    + D1P2_m23_main*(H2f3/Nf3) + D1G2_m23_main*(Y2f3/Nf3) + D1X2_m23_main*(Z2f3/Nf3) + D1K2_m23_main*(A2f3/Nf3)
    # infected partner on ART
    + D1P3_m23_main*(H3f3/Nf3) + D1G3_m23_main*(Y3f3/Nf3) + D1X3_m23_main*(Z3f3/Nf3) + D1K3_m23_main*(A3f3/Nf3)
  )
  # regular
  +(dm22_reg*rho_m2_2_reg)*(
    # infected partner not on HCT/ART
      D1P1_m22_reg*(H1f2/Nf2) + D1G1_m22_reg*(Y1f2/Nf2) + D1X1_m22_reg*(Z1f2/Nf2) + D1K1_m22_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D1P2_m22_reg*(H2f2/Nf2) + D1G2_m22_reg*(Y2f2/Nf2) + D1X2_m22_reg*(Z2f2/Nf2) + D1K2_m22_reg*(A2f2/Nf2)
    # infected partner on ART
    + D1P3_m22_reg*(H3f2/Nf2) + D1G3_m22_reg*(Y3f2/Nf2) + D1X3_m22_reg*(Z3f2/Nf2) + D1K3_m22_reg*(A3f2/Nf2)
  )
  +(dm23_reg*rho_m2_3_reg)*(
    # infected partner not on HCT/ART
      D1P1_m23_reg*(H1f3/Nf3) + D1G1_m23_reg*(Y1f3/Nf3) + D1X1_m23_reg*(Z1f3/Nf3) + D1K1_m23_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D1P2_m23_reg*(H2f3/Nf3) + D1G2_m23_reg*(Y2f3/Nf3) + D1X2_m23_reg*(Z2f3/Nf3) + D1K2_m23_reg*(A2f3/Nf3)
    # infected partner on ART
    + D1P3_m23_reg*(H3f3/Nf3) + D1G3_m23_reg*(Y3f3/Nf3) + D1X3_m23_reg*(Z3f3/Nf3) + D1K3_m23_reg*(A3f3/Nf3)
  )
)

pi1_m3 = (
  # main
  (dm31_main*rho_m3_1_main)*(
    # infected partner not on HCT/ART
      D1P1_m31_main*(H1f1/Nf1) + D1G1_m31_main*(Y1f1/Nf1) + D1X1_m31_main*(Z1f1/Nf1) + D1K1_m31_main*(A1f1/Nf1)
    # infected partner on HCT
    + D1P2_m31_main*(H2f1/Nf1) + D1G2_m31_main*(Y2f1/Nf1) + D1X2_m31_main*(Z2f1/Nf1) + D1K2_m31_main*(A2f1/Nf1)
    # infected partner on ART
    + D1P3_m31_main*(H3f1/Nf1) + D1G3_m31_main*(Y3f1/Nf1) + D1X3_m31_main*(Z3f1/Nf1) + D1K3_m31_main*(A3f1/Nf1)
  )
  +(dm32_main*rho_m3_2_main)*(
    # infected partner not on HCT/ART
      D1P1_m32_main*(H1f2/Nf2) + D1G1_m32_main*(Y1f2/Nf2) + D1X1_m32_main*(Z1f2/Nf2) + D1K1_m32_main*(A1f2/Nf2)
    # infected partner on HCT
    + D1P2_m32_main*(H2f2/Nf2) + D1G2_m32_main*(Y2f2/Nf2) + D1X2_m32_main*(Z2f2/Nf2) + D1K2_m32_main*(A2f2/Nf2)
    # infected partner on ART
    + D1P3_m32_main*(H3f2/Nf2) + D1G3_m32_main*(Y3f2/Nf2) + D1X3_m32_main*(Z3f2/Nf2) + D1K3_m32_main*(A3f2/Nf2)
  )
  +(dm33_main*rho_m3_3_main)*(
    # infected partner not on HCT/ART
      D1P1_m33_main*(H1f3/Nf3) + D1G1_m33_main*(Y1f3/Nf3) + D1X1_m33_main*(Z1f3/Nf3) + D1K1_m33_main*(A1f3/Nf3)
    # infected partner on HCT
    + D1P2_m33_main*(H2f3/Nf3) + D1G2_m33_main*(Y2f3/Nf3) + D1X2_m33_main*(Z2f3/Nf3) + D1K2_m33_main*(A2f3/Nf3)
    # infected partner on ART
    + D1P3_m33_main*(H3f3/Nf3) + D1G3_m33_main*(Y3f3/Nf3) + D1X3_m33_main*(Z3f3/Nf3) + D1K3_m33_main*(A3f3/Nf3)
  )
  # regular
  +(dm32_reg*rho_m3_2_reg)*(
    # infected partner not on HCT/ART
      D1P1_m32_reg*(H1f2/Nf2) + D1G1_m32_reg*(Y1f2/Nf2) + D1X1_m32_reg*(Z1f2/Nf2) + D1K1_m32_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D1P2_m32_reg*(H2f2/Nf2) + D1G2_m32_reg*(Y2f2/Nf2) + D1X2_m32_reg*(Z2f2/Nf2) + D1K2_m32_reg*(A2f2/Nf2)
    # infected partner on ART
    + D1P3_m32_reg*(H3f2/Nf2) + D1G3_m32_reg*(Y3f2/Nf2) + D1X3_m32_reg*(Z3f2/Nf2) + D1K3_m32_reg*(A3f2/Nf2)
  )
  +(dm33_reg*rho_m3_3_reg)*(
    # infected partner not on HCT/ART
      D1P1_m33_reg*(H1f3/Nf3) + D1G1_m33_reg*(Y1f3/Nf3) + D1X1_m33_reg*(Z1f3/Nf3) + D1K1_m33_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D1P2_m33_reg*(H2f3/Nf3) + D1G2_m33_reg*(Y2f3/Nf3) + D1X2_m33_reg*(Z2f3/Nf3) + D1K2_m33_reg*(A2f3/Nf3)
    # infected partner on ART
    + D1P3_m33_reg*(H3f3/Nf3) + D1G3_m33_reg*(Y3f3/Nf3) + D1X3_m33_reg*(Z3f3/Nf3) + D1K3_m33_reg*(A3f3/Nf3)
  )
  # casual
  +(dm33_casual*rho_m3_3_casual)*(
    # infected partner not on HCT/ART
      D1P1_m33_casual*(H1f3/Nf3) + D1G1_m33_casual*(Y1f3/Nf3) + D1X1_m33_casual*(Z1f3/Nf3) + D1K1_m33_casual*(A1f3/Nf3)
    # infected partner on HCT
    + D1P2_m33_casual*(H2f3/Nf3) + D1G2_m33_casual*(Y2f3/Nf3) + D1X2_m33_casual*(Z2f3/Nf3) + D1K2_m33_casual*(A2f3/Nf3)
    # infected partner on ART
    + D1P3_m33_casual*(H3f3/Nf3) + D1G3_m33_casual*(Y3f3/Nf3) + D1X3_m33_casual*(Z3f3/Nf3) + D1K3_m33_casual*(A3f3/Nf3)
  )
)


# Females on PrEP
pi2_f1 = (
  # low-low
  (df11_main*rho_f1_1_main)*(
    # infected partner not on HCT/ART
      D2P1_f11_main*(H1m1/Nm1) + D2G1_f11_main*(Y1m1/Nm1) + D2X1_f11_main*(Z1m1/Nm1) + D2K1_f11_main*(A1m1/Nm1)
    # infected partner on HCT
    + D2P2_f11_main*(H2m1/Nm1) + D2G2_f11_main*(Y2m1/Nm1) + D2X2_f11_main*(Z2m1/Nm1) + D2K2_f11_main*(A2m1/Nm1)
    # infected partner on ART
    + D2P3_f11_main*(H3m1/Nm1) + D2G3_f11_main*(Y3m1/Nm1) + D2X3_f11_main*(Z3m1/Nm1) + D2K3_f11_main*(A3m1/Nm1)
  )
  # low-medium
  +(df12_main*rho_f1_2_main)*(
    # infected partner not on HCT/ART
      D2P1_f12_main*(H1m2/Nm2) + D2G1_f12_main*(Y1m2/Nm2) + D2X1_f12_main*(Z1m2/Nm2) + D2K1_f12_main*(A1m2/Nm2)
    # infected partner on HCT
    + D2P2_f12_main*(H2m2/Nm2) + D2G2_f12_main*(Y2m2/Nm2) + D2X2_f12_main*(Z2m2/Nm2) + D2K2_f12_main*(A2m2/Nm2)
    # infected partner on ART
    + D2P3_f12_main*(H3m2/Nm2) + D2G3_f12_main*(Y3m2/Nm2) + D2X3_f12_main*(Z3m2/Nm2) + D2K3_f12_main*(A3m2/Nm2)
  )
  # low-high
  +(df13_main*rho_f1_3_main)*(
    # infected partner not on HCT/ART
      D2P1_f13_main*(H1m3/Nm3) + D2G1_f13_main*(Y1m3/Nm3) + D2X1_f13_main*(Z1m3/Nm3) + D2K1_f13_main*(A1m3/Nm3)
    # infected partner on HCT
    + D2P2_f13_main*(H2m3/Nm3) + D2G2_f13_main*(Y2m3/Nm3) + D2X2_f13_main*(Z2m3/Nm3) + D2K2_f13_main*(A2m3/Nm3)
    # infected partner on ART
    + D2P3_f13_main*(H3m3/Nm3) + D2G3_f13_main*(Y3m3/Nm3) + D2X3_f13_main*(Z3m3/Nm3) + D2K3_f13_main*(A3m3/Nm3)
  )
)

pi2_f2 = (
  # main
  (df21_main*rho_f2_1_main)*(
    # infected partner not on HCT/ART
      D2P1_f21_main*(H1m1/Nm1) + D2G1_f21_main*(Y1m1/Nm1) + D2X1_f21_main*(Z1m1/Nm1) + D2K1_f21_main*(A1m1/Nm1)
    # infected partner on HCT
    + D2P2_f21_main*(H2m1/Nm1) + D2G2_f21_main*(Y2m1/Nm1) + D2X2_f21_main*(Z2m1/Nm1) + D2K2_f21_main*(A2m1/Nm1)
    # infected partner on ART
    + D2P3_f21_main*(H3m1/Nm1) + D2G3_f21_main*(Y3m1/Nm1) + D2X3_f21_main*(Z3m1/Nm1) + D2K3_f21_main*(A3m1/Nm1)
  )
  +(df22_main*rho_f2_2_main)*(
    # infected partner not on HCT/ART
      D2P1_f22_main*(H1m2/Nm2) + D2G1_f22_main*(Y1m2/Nm2) + D2X1_f22_main*(Z1m2/Nm2) + D2K1_f22_main*(A1m2/Nm2)
    # infected partner on HCT
    + D2P2_f22_main*(H2m2/Nm2) + D2G2_f22_main*(Y2m2/Nm2) + D2X2_f22_main*(Z2m2/Nm2) + D2K2_f22_main*(A2m2/Nm2)
    # infected partner on ART
    + D2P3_f22_main*(H3m2/Nm2) + D2G3_f22_main*(Y3m2/Nm2) + D2X3_f22_main*(Z3m2/Nm2) + D2K3_f22_main*(A3m2/Nm2)
  )
  +(df23_main*rho_f2_3_main)*(
    # infected partner not on HCT/ART
      D2P1_f23_main*(H1m3/Nm3) + D2G1_f23_main*(Y1m3/Nm3) + D2X1_f23_main*(Z1m3/Nm3) + D2K1_f23_main*(A1m3/Nm3)
    # infected partner on HCT
    + D2P2_f23_main*(H2m3/Nm3) + D2G2_f23_main*(Y2m3/Nm3) + D2X2_f23_main*(Z2m3/Nm3) + D2K2_f23_main*(A2m3/Nm3)
    # infected partner on ART
    + D2P3_f23_main*(H3m3/Nm3) + D2G3_f23_main*(Y3m3/Nm3) + D2X3_f23_main*(Z3m3/Nm3) + D2K3_f23_main*(A3m3/Nm3)
  )
  # regular
  +(df22_reg*rho_f2_2_reg)*(
    # infected partner not on HCT/ART
      D2P1_f22_reg*(H1m2/Nm2) + D2G1_f22_reg*(Y1m2/Nm2) + D2X1_f22_reg*(Z1m2/Nm2) + D2K1_f22_reg*(A1m2/Nm2)
    # infected partner on HCT
    + D2P2_f22_reg*(H2m2/Nm2) + D2G2_f22_reg*(Y2m2/Nm2) + D2X2_f22_reg*(Z2m2/Nm2) + D2K2_f22_reg*(A2m2/Nm2)
    # infected partner on ART
    + D2P3_f22_reg*(H3m2/Nm2) + D2G3_f22_reg*(Y3m2/Nm2) + D2X3_f22_reg*(Z3m2/Nm2) + D2K3_f22_reg*(A3m2/Nm2)
  )
  +(df23_reg*rho_f2_3_reg)*(
    # infected partner not on HCT/ART
      D2P1_f23_reg*(H1m3/Nm3) + D2G1_f23_reg*(Y1m3/Nm3) + D2X1_f23_reg*(Z1m3/Nm3) + D2K1_f23_reg*(A1m3/Nm3)
    # infected partner on HCT
    + D2P2_f23_reg*(H2m3/Nm3) + D2G2_f23_reg*(Y2m3/Nm3) + D2X2_f23_reg*(Z2m3/Nm3) + D2K2_f23_reg*(A2m3/Nm3)
    # infected partner on ART
    + D2P3_f23_reg*(H3m3/Nm3) + D2G3_f23_reg*(Y3m3/Nm3) + D2X3_f23_reg*(Z3m3/Nm3) + D2K3_f23_reg*(A3m3/Nm3)
  )
)

pi2_f3 = (
  # main
  (df31_main*rho_f3_1_main)*(
    # infected partner not on HCT/ART
      D2P1_f31_main*(H1m1/Nm1) + D2G1_f31_main*(Y1m1/Nm1) + D2X1_f31_main*(Z1m1/Nm1) + D2K1_f31_main*(A1m1/Nm1)
    # infected partner on HCT
    + D2P2_f31_main*(H2m1/Nm1) + D2G2_f31_main*(Y2m1/Nm1) + D2X2_f31_main*(Z2m1/Nm1) + D2K2_f31_main*(A2m1/Nm1)
    # infected partner on ART
    + D2P3_f31_main*(H3m1/Nm1) + D2G3_f31_main*(Y3m1/Nm1) + D2X3_f31_main*(Z3m1/Nm1) + D2K3_f31_main*(A3m1/Nm1)
  )
  +(df32_main*rho_f3_2_main)*(
    # infected partner not on HCT/ART
      D2P1_f32_main*(H1m2/Nm2) + D2G1_f32_main*(Y1m2/Nm2) + D2X1_f32_main*(Z1m2/Nm2) + D2K1_f32_main*(A1m2/Nm2)
    # infected partner on HCT
    + D2P2_f32_main*(H2m2/Nm2) + D2G2_f32_main*(Y2m2/Nm2) + D2X2_f32_main*(Z2m2/Nm2) + D2K2_f32_main*(A2m2/Nm2)
    # infected partner on ART
    + D2P3_f32_main*(H3m2/Nm2) + D2G3_f32_main*(Y3m2/Nm2) + D2X3_f32_main*(Z3m2/Nm2) + D2K3_f32_main*(A3m2/Nm2)
  )
  +(df33_main*rho_f3_3_main)*(
    # infected partner not on HCT/ART
      D2P1_f33_main*(H1m3/Nm3) + D2G1_f33_main*(Y1m3/Nm3) + D2X1_f33_main*(Z1m3/Nm3) + D2K1_f33_main*(A1m3/Nm3)
    # infected partner on HCT
    + D2P2_f33_main*(H2m3/Nm3) + D2G2_f33_main*(Y2m3/Nm3) + D2X2_f33_main*(Z2m3/Nm3) + D2K2_f33_main*(A2m3/Nm3)
    # infected partner on ART
    + D2P3_f33_main*(H3m3/Nm3) + D2G3_f33_main*(Y3m3/Nm3) + D2X3_f33_main*(Z3m3/Nm3) + D2K3_f33_main*(A3m3/Nm3)
  )
  # regular
  +(df32_reg*rho_f3_2_reg)*(
    # infected partner not on HCT/ART
      D2P1_f32_reg*(H1m2/Nm2) + D2G1_f32_reg*(Y1m2/Nm2) + D2X1_f32_reg*(Z1m2/Nm2) + D2K1_f32_reg*(A1m2/Nm2)
    # infected partner on HCT
    + D2P2_f32_reg*(H2m2/Nm2) + D2G2_f32_reg*(Y2m2/Nm2) + D2X2_f32_reg*(Z2m2/Nm2) + D2K2_f32_reg*(A2m2/Nm2)
    # infected partner on ART
    + D2P3_f32_reg*(H3m2/Nm2) + D2G3_f32_reg*(Y3m2/Nm2) + D2X3_f32_reg*(Z3m2/Nm2) + D2K3_f32_reg*(A3m2/Nm2)
  )
  +(df33_reg*rho_f3_3_reg)*(
    # infected partner not on HCT/ART
      D2P1_f33_reg*(H1m3/Nm3) + D2G1_f33_reg*(Y1m3/Nm3) + D2X1_f33_reg*(Z1m3/Nm3) + D2K1_f33_reg*(A1m3/Nm3)
    # infected partner on HCT
    + D2P2_f33_reg*(H2m3/Nm3) + D2G2_f33_reg*(Y2m3/Nm3) + D2X2_f33_reg*(Z2m3/Nm3) + D2K2_f33_reg*(A2m3/Nm3)
    # infected partner on ART
    + D2P3_f33_reg*(H3m3/Nm3) + D2G3_f33_reg*(Y3m3/Nm3) + D2X3_f33_reg*(Z3m3/Nm3) + D2K3_f33_reg*(A3m3/Nm3)
  )
  # casual
  +(df33_casual*rho_f3_3_casual)*(
    # infected partner not on HCT/ART
      D2P1_f33_casual*(H1m3/Nm3) + D2G1_f33_casual*(Y1m3/Nm3) + D2X1_f33_casual*(Z1m3/Nm3) + D2K1_f33_casual*(A1m3/Nm3)
    # infected partner on HCT
    + D2P2_f33_casual*(H2m3/Nm3) + D2G2_f33_casual*(Y2m3/Nm3) + D2X2_f33_casual*(Z2m3/Nm3) + D2K2_f33_casual*(A2m3/Nm3)
    # infected partner on ART
    + D2P3_f33_casual*(H3m3/Nm3) + D2G3_f33_casual*(Y3m3/Nm3) + D2X3_f33_casual*(Z3m3/Nm3) + D2K3_f33_casual*(A3m3/Nm3)
  )
)

###
# Uncircumcised males on PrEP
pi2_m1 = (
  # low-low
  (dm11_main*rho_m1_1_main)*(
    # infected partner not on HCT/ART
      D2P1_m11_main*(H1f1/Nf1) + D2G1_m11_main*(Y1f1/Nf1) + D2X1_m11_main*(Z1f1/Nf1) + D2K1_m11_main*(A1f1/Nf1)
    # infected partner on HCT
    + D2P2_m11_main*(H2f1/Nf1) + D2G2_m11_main*(Y2f1/Nf1) + D2X2_m11_main*(Z2f1/Nf1) + D2K2_m11_main*(A2f1/Nf1)
    # infected partner on ART
    + D2P3_m11_main*(H3f1/Nf1) + D2G3_m11_main*(Y3f1/Nf1) + D2X3_m11_main*(Z3f1/Nf1) + D2K3_m11_main*(A3f1/Nf1)
  )
  # low-medium
  +(dm12_main*rho_m1_2_main)*(
    # infected partner not on HCT/ART
      D2P1_m12_main*(H1f2/Nf2) + D2G1_m12_main*(Y1f2/Nf2) + D2X1_m12_main*(Z1f2/Nf2) + D2K1_m12_main*(A1f2/Nf2)
    # infected partner on HCT
    + D2P2_m12_main*(H2f2/Nf2) + D2G2_m12_main*(Y2f2/Nf2) + D2X2_m12_main*(Z2f2/Nf2) + D2K2_m12_main*(A2f2/Nf2)
    # infected partner on ART
    + D2P3_m12_main*(H3f2/Nf2) + D2G3_m12_main*(Y3f2/Nf2) + D2X3_m12_main*(Z3f2/Nf2) + D2K3_m12_main*(A3f2/Nf2)
  )
  # low-high
  +(dm13_main*rho_m1_3_main)*(
    # infected partner not on HCT/ART
      D2P1_m13_main*(H1f3/Nf3) + D2G1_m13_main*(Y1f3/Nf3) + D2X1_m13_main*(Z1f3/Nf3) + D2K1_m13_main*(A1f3/Nf3)
    # infected partner on HCT
    + D2P2_m13_main*(H2f3/Nf3) + D2G2_m13_main*(Y2f3/Nf3) + D2X2_m13_main*(Z2f3/Nf3) + D2K2_m13_main*(A2f3/Nf3)
    # infected partner on ART
    + D2P3_m13_main*(H3f3/Nf3) + D2G3_m13_main*(Y3f3/Nf3) + D2X3_m13_main*(Z3f3/Nf3) + D2K3_m13_main*(A3f3/Nf3)
  )
)

pi2_m2 = (
  # main
  (dm21_main*rho_m2_1_main)*(
    # infected partner not on HCT/ART
      D2P1_m21_main*(H1f1/Nf1) + D2G1_m21_main*(Y1f1/Nf1) + D2X1_m21_main*(Z1f1/Nf1) + D2K1_m21_main*(A1f1/Nf1)
    # infected partner on HCT
    + D2P2_m21_main*(H2f1/Nf1) + D2G2_m21_main*(Y2f1/Nf1) + D2X2_m21_main*(Z2f1/Nf1) + D2K2_m21_main*(A2f1/Nf1)
    # infected partner on ART
    + D2P3_m21_main*(H3f1/Nf1) + D2G3_m21_main*(Y3f1/Nf1) + D2X3_m21_main*(Z3f1/Nf1) + D2K3_m21_main*(A3f1/Nf1)
  )
  +(dm22_main*rho_m2_2_main)*(
    # infected partner not on HCT/ART
      D2P1_m22_main*(H1f2/Nf2) + D2G1_m22_main*(Y1f2/Nf2) + D2X1_m22_main*(Z1f2/Nf2) + D2K1_m22_main*(A1f2/Nf2)
    # infected partner on HCT
    + D2P2_m22_main*(H2f2/Nf2) + D2G2_m22_main*(Y2f2/Nf2) + D2X2_m22_main*(Z2f2/Nf2) + D2K2_m22_main*(A2f2/Nf2)
    # infected partner on ART
    + D2P3_m22_main*(H3f2/Nf2) + D2G3_m22_main*(Y3f2/Nf2) + D2X3_m22_main*(Z3f2/Nf2) + D2K3_m22_main*(A3f2/Nf2)
  )
  +(dm23_main*rho_m2_3_main)*(
    # infected partner not on HCT/ART
      D2P1_m23_main*(H1f3/Nf3) + D2G1_m23_main*(Y1f3/Nf3) + D2X1_m23_main*(Z1f3/Nf3) + D2K1_m23_main*(A1f3/Nf3)
    # infected partner on HCT
    + D2P2_m23_main*(H2f3/Nf3) + D2G2_m23_main*(Y2f3/Nf3) + D2X2_m23_main*(Z2f3/Nf3) + D2K2_m23_main*(A2f3/Nf3)
    # infected partner on ART
    + D2P3_m23_main*(H3f3/Nf3) + D2G3_m23_main*(Y3f3/Nf3) + D2X3_m23_main*(Z3f3/Nf3) + D2K3_m23_main*(A3f3/Nf3)
  )
  # regular
  +(dm22_reg*rho_m2_2_reg)*(
    # infected partner not on HCT/ART
      D2P1_m22_reg*(H1f2/Nf2) + D2G1_m22_reg*(Y1f2/Nf2) + D2X1_m22_reg*(Z1f2/Nf2) + D2K1_m22_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D2P2_m22_reg*(H2f2/Nf2) + D2G2_m22_reg*(Y2f2/Nf2) + D2X2_m22_reg*(Z2f2/Nf2) + D2K2_m22_reg*(A2f2/Nf2)
    # infected partner on ART
    + D2P3_m22_reg*(H3f2/Nf2) + D2G3_m22_reg*(Y3f2/Nf2) + D2X3_m22_reg*(Z3f2/Nf2) + D2K3_m22_reg*(A3f2/Nf2)
  )
  +(dm23_reg*rho_m2_3_reg)*(
    # infected partner not on HCT/ART
      D2P1_m23_reg*(H1f3/Nf3) + D2G1_m23_reg*(Y1f3/Nf3) + D2X1_m23_reg*(Z1f3/Nf3) + D2K1_m23_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D2P2_m23_reg*(H2f3/Nf3) + D2G2_m23_reg*(Y2f3/Nf3) + D2X2_m23_reg*(Z2f3/Nf3) + D2K2_m23_reg*(A2f3/Nf3)
    # infected partner on ART
    + D2P3_m23_reg*(H3f3/Nf3) + D2G3_m23_reg*(Y3f3/Nf3) + D2X3_m23_reg*(Z3f3/Nf3) + D2K3_m23_reg*(A3f3/Nf3)
  )
)

pi2_m3 = (
  # main
  (dm31_main*rho_m3_1_main)*(
    # infected partner not on HCT/ART
      D2P1_m31_main*(H1f1/Nf1) + D2G1_m31_main*(Y1f1/Nf1) + D2X1_m31_main*(Z1f1/Nf1) + D2K1_m31_main*(A1f1/Nf1)
    # infected partner on HCT
    + D2P2_m31_main*(H2f1/Nf1) + D2G2_m31_main*(Y2f1/Nf1) + D2X2_m31_main*(Z2f1/Nf1) + D2K2_m31_main*(A2f1/Nf1)
    # infected partner on ART
    + D2P3_m31_main*(H3f1/Nf1) + D2G3_m31_main*(Y3f1/Nf1) + D2X3_m31_main*(Z3f1/Nf1) + D2K3_m31_main*(A3f1/Nf1)
  )
  +(dm32_main*rho_m3_2_main)*(
    # infected partner not on HCT/ART
      D2P1_m32_main*(H1f2/Nf2) + D2G1_m32_main*(Y1f2/Nf2) + D2X1_m32_main*(Z1f2/Nf2) + D2K1_m32_main*(A1f2/Nf2)
    # infected partner on HCT
    + D2P2_m32_main*(H2f2/Nf2) + D2G2_m32_main*(Y2f2/Nf2) + D2X2_m32_main*(Z2f2/Nf2) + D2K2_m32_main*(A2f2/Nf2)
    # infected partner on ART
    + D2P3_m32_main*(H3f2/Nf2) + D2G3_m32_main*(Y3f2/Nf2) + D2X3_m32_main*(Z3f2/Nf2) + D2K3_m32_main*(A3f2/Nf2)
  )
  +(dm33_main*rho_m3_3_main)*(
    # infected partner not on HCT/ART
      D2P1_m33_main*(H1f3/Nf3) + D2G1_m33_main*(Y1f3/Nf3) + D2X1_m33_main*(Z1f3/Nf3) + D2K1_m33_main*(A1f3/Nf3)
    # infected partner on HCT
    + D2P2_m33_main*(H2f3/Nf3) + D2G2_m33_main*(Y2f3/Nf3) + D2X2_m33_main*(Z2f3/Nf3) + D2K2_m33_main*(A2f3/Nf3)
    # infected partner on ART
    + D2P3_m33_main*(H3f3/Nf3) + D2G3_m33_main*(Y3f3/Nf3) + D2X3_m33_main*(Z3f3/Nf3) + D2K3_m33_main*(A3f3/Nf3)
  )
  # regular
  +(dm32_reg*rho_m3_2_reg)*(
    # infected partner not on HCT/ART
      D2P1_m32_reg*(H1f2/Nf2) + D2G1_m32_reg*(Y1f2/Nf2) + D2X1_m32_reg*(Z1f2/Nf2) + D2K1_m32_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D2P2_m32_reg*(H2f2/Nf2) + D2G2_m32_reg*(Y2f2/Nf2) + D2X2_m32_reg*(Z2f2/Nf2) + D2K2_m32_reg*(A2f2/Nf2)
    # infected partner on ART
    + D2P3_m32_reg*(H3f2/Nf2) + D2G3_m32_reg*(Y3f2/Nf2) + D2X3_m32_reg*(Z3f2/Nf2) + D2K3_m32_reg*(A3f2/Nf2)
  )
  +(dm33_reg*rho_m3_3_reg)*(
    # infected partner not on HCT/ART
      D2P1_m33_reg*(H1f3/Nf3) + D2G1_m33_reg*(Y1f3/Nf3) + D2X1_m33_reg*(Z1f3/Nf3) + D2K1_m33_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D2P2_m33_reg*(H2f3/Nf3) + D2G2_m33_reg*(Y2f3/Nf3) + D2X2_m33_reg*(Z2f3/Nf3) + D2K2_m33_reg*(A2f3/Nf3)
    # infected partner on ART
    + D2P3_m33_reg*(H3f3/Nf3) + D2G3_m33_reg*(Y3f3/Nf3) + D2X3_m33_reg*(Z3f3/Nf3) + D2K3_m33_reg*(A3f3/Nf3)
  )
  # casual
  +(dm33_casual*rho_m3_3_casual)*(
    # infected partner not on HCT/ART
      D2P1_m33_casual*(H1f3/Nf3) + D2G1_m33_casual*(Y1f3/Nf3) + D2X1_m33_casual*(Z1f3/Nf3) + D2K1_m33_casual*(A1f3/Nf3)
    # infected partner on HCT
    + D2P2_m33_casual*(H2f3/Nf3) + D2G2_m33_casual*(Y2f3/Nf3) + D2X2_m33_casual*(Z2f3/Nf3) + D2K2_m33_casual*(A2f3/Nf3)
    # infected partner on ART
    + D2P3_m33_casual*(H3f3/Nf3) + D2G3_m33_casual*(Y3f3/Nf3) + D2X3_m33_casual*(Z3f3/Nf3) + D2K3_m33_casual*(A3f3/Nf3)
  )
)

###
# Circumcised males not on PrEP
pi3_m1 = (
  # low-low
  (dm11_main*rho_m1_1_main)*(
    # infected partner not on HCT/ART
      D3P1_m11_main*(H1f1/Nf1) + D3G1_m11_main*(Y1f1/Nf1) + D3X1_m11_main*(Z1f1/Nf1) + D3K1_m11_main*(A1f1/Nf1)
    # infected partner on HCT
    + D3P2_m11_main*(H2f1/Nf1) + D3G2_m11_main*(Y2f1/Nf1) + D3X2_m11_main*(Z2f1/Nf1) + D3K2_m11_main*(A2f1/Nf1)
    # infected partner on ART
    + D3P3_m11_main*(H3f1/Nf1) + D3G3_m11_main*(Y3f1/Nf1) + D3X3_m11_main*(Z3f1/Nf1) + D3K3_m11_main*(A3f1/Nf1)
  )
  # low-medium
  +(dm12_main*rho_m1_2_main)*(
    # infected partner not on HCT/ART
      D3P1_m12_main*(H1f2/Nf2) + D3G1_m12_main*(Y1f2/Nf2) + D3X1_m12_main*(Z1f2/Nf2) + D3K1_m12_main*(A1f2/Nf2)
    # infected partner on HCT
    + D3P2_m12_main*(H2f2/Nf2) + D3G2_m12_main*(Y2f2/Nf2) + D3X2_m12_main*(Z2f2/Nf2) + D3K2_m12_main*(A2f2/Nf2)
    # infected partner on ART
    + D3P3_m12_main*(H3f2/Nf2) + D3G3_m12_main*(Y3f2/Nf2) + D3X3_m12_main*(Z3f2/Nf2) + D3K3_m12_main*(A3f2/Nf2)
  )
  # low-high
  +(dm13_main*rho_m1_3_main)*(
    # infected partner not on HCT/ART
      D3P1_m13_main*(H1f3/Nf3) + D3G1_m13_main*(Y1f3/Nf3) + D3X1_m13_main*(Z1f3/Nf3) + D3K1_m13_main*(A1f3/Nf3)
    # infected partner on HCT
    + D3P2_m13_main*(H2f3/Nf3) + D3G2_m13_main*(Y2f3/Nf3) + D3X2_m13_main*(Z2f3/Nf3) + D3K2_m13_main*(A2f3/Nf3)
    # infected partner on ART
    + D3P3_m13_main*(H3f3/Nf3) + D3G3_m13_main*(Y3f3/Nf3) + D3X3_m13_main*(Z3f3/Nf3) + D3K3_m13_main*(A3f3/Nf3)
  )
)

pi3_m2 = (
  # main
  (dm21_main*rho_m2_1_main)*(
    # infected partner not on HCT/ART
      D3P1_m21_main*(H1f1/Nf1) + D3G1_m21_main*(Y1f1/Nf1) + D3X1_m21_main*(Z1f1/Nf1) + D3K1_m21_main*(A1f1/Nf1)
    # infected partner on HCT
    + D3P2_m21_main*(H2f1/Nf1) + D3G2_m21_main*(Y2f1/Nf1) + D3X2_m21_main*(Z2f1/Nf1) + D3K2_m21_main*(A2f1/Nf1)
    # infected partner on ART
    + D3P3_m21_main*(H3f1/Nf1) + D3G3_m21_main*(Y3f1/Nf1) + D3X3_m21_main*(Z3f1/Nf1) + D3K3_m21_main*(A3f1/Nf1)
  )
  +(dm22_main*rho_m2_2_main)*(
    # infected partner not on HCT/ART
      D3P1_m22_main*(H1f2/Nf2) + D3G1_m22_main*(Y1f2/Nf2) + D3X1_m22_main*(Z1f2/Nf2) + D3K1_m22_main*(A1f2/Nf2)
    # infected partner on HCT
    + D3P2_m22_main*(H2f2/Nf2) + D3G2_m22_main*(Y2f2/Nf2) + D3X2_m22_main*(Z2f2/Nf2) + D3K2_m22_main*(A2f2/Nf2)
    # infected partner on ART
    + D3P3_m22_main*(H3f2/Nf2) + D3G3_m22_main*(Y3f2/Nf2) + D3X3_m22_main*(Z3f2/Nf2) + D3K3_m22_main*(A3f2/Nf2)
  )
  +(dm23_main*rho_m2_3_main)*(
    # infected partner not on HCT/ART
      D3P1_m23_main*(H1f3/Nf3) + D3G1_m23_main*(Y1f3/Nf3) + D3X1_m23_main*(Z1f3/Nf3) + D3K1_m23_main*(A1f3/Nf3)
    # infected partner on HCT
    + D3P2_m23_main*(H2f3/Nf3) + D3G2_m23_main*(Y2f3/Nf3) + D3X2_m23_main*(Z2f3/Nf3) + D3K2_m23_main*(A2f3/Nf3)
    # infected partner on ART
    + D3P3_m23_main*(H3f3/Nf3) + D3G3_m23_main*(Y3f3/Nf3) + D3X3_m23_main*(Z3f3/Nf3) + D3K3_m23_main*(A3f3/Nf3)
  )
  # regular
  +(dm22_reg*rho_m2_2_reg)*(
    # infected partner not on HCT/ART
      D3P1_m22_reg*(H1f2/Nf2) + D3G1_m22_reg*(Y1f2/Nf2) + D3X1_m22_reg*(Z1f2/Nf2) + D3K1_m22_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D3P2_m22_reg*(H2f2/Nf2) + D3G2_m22_reg*(Y2f2/Nf2) + D3X2_m22_reg*(Z2f2/Nf2) + D3K2_m22_reg*(A2f2/Nf2)
    # infected partner on ART
    + D3P3_m22_reg*(H3f2/Nf2) + D3G3_m22_reg*(Y3f2/Nf2) + D3X3_m22_reg*(Z3f2/Nf2) + D3K3_m22_reg*(A3f2/Nf2)
  )
  +(dm23_reg*rho_m2_3_reg)*(
    # infected partner not on HCT/ART
      D3P1_m23_reg*(H1f3/Nf3) + D3G1_m23_reg*(Y1f3/Nf3) + D3X1_m23_reg*(Z1f3/Nf3) + D3K1_m23_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D3P2_m23_reg*(H2f3/Nf3) + D3G2_m23_reg*(Y2f3/Nf3) + D3X2_m23_reg*(Z2f3/Nf3) + D3K2_m23_reg*(A2f3/Nf3)
    # infected partner on ART
    + D3P3_m23_reg*(H3f3/Nf3) + D3G3_m23_reg*(Y3f3/Nf3) + D3X3_m23_reg*(Z3f3/Nf3) + D3K3_m23_reg*(A3f3/Nf3)
  )
)

pi3_m3 = (
  # main
  (dm31_main*rho_m3_1_main)*(
    # infected partner not on HCT/ART
      D3P1_m31_main*(H1f1/Nf1) + D3G1_m31_main*(Y1f1/Nf1) + D3X1_m31_main*(Z1f1/Nf1) + D3K1_m31_main*(A1f1/Nf1)
    # infected partner on HCT
    + D3P2_m31_main*(H2f1/Nf1) + D3G2_m31_main*(Y2f1/Nf1) + D3X2_m31_main*(Z2f1/Nf1) + D3K2_m31_main*(A2f1/Nf1)
    # infected partner on ART
    + D3P3_m31_main*(H3f1/Nf1) + D3G3_m31_main*(Y3f1/Nf1) + D3X3_m31_main*(Z3f1/Nf1) + D3K3_m31_main*(A3f1/Nf1)
  )
  +(dm32_main*rho_m3_2_main)*(
    # infected partner not on HCT/ART
      D3P1_m32_main*(H1f2/Nf2) + D3G1_m32_main*(Y1f2/Nf2) + D3X1_m32_main*(Z1f2/Nf2) + D3K1_m32_main*(A1f2/Nf2)
    # infected partner on HCT
    + D3P2_m32_main*(H2f2/Nf2) + D3G2_m32_main*(Y2f2/Nf2) + D3X2_m32_main*(Z2f2/Nf2) + D3K2_m32_main*(A2f2/Nf2)
    # infected partner on ART
    + D3P3_m32_main*(H3f2/Nf2) + D3G3_m32_main*(Y3f2/Nf2) + D3X3_m32_main*(Z3f2/Nf2) + D3K3_m32_main*(A3f2/Nf2)
  )
  +(dm33_main*rho_m3_3_main)*(
    # infected partner not on HCT/ART
      D3P1_m33_main*(H1f3/Nf3) + D3G1_m33_main*(Y1f3/Nf3) + D3X1_m33_main*(Z1f3/Nf3) + D3K1_m33_main*(A1f3/Nf3)
    # infected partner on HCT
    + D3P2_m33_main*(H2f3/Nf3) + D3G2_m33_main*(Y2f3/Nf3) + D3X2_m33_main*(Z2f3/Nf3) + D3K2_m33_main*(A2f3/Nf3)
    # infected partner on ART
    + D3P3_m33_main*(H3f3/Nf3) + D3G3_m33_main*(Y3f3/Nf3) + D3X3_m33_main*(Z3f3/Nf3) + D3K3_m33_main*(A3f3/Nf3)
  )
  # regular
  +(dm32_reg*rho_m3_2_reg)*(
    # infected partner not on HCT/ART
      D3P1_m32_reg*(H1f2/Nf2) + D3G1_m32_reg*(Y1f2/Nf2) + D3X1_m32_reg*(Z1f2/Nf2) + D3K1_m32_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D3P2_m32_reg*(H2f2/Nf2) + D3G2_m32_reg*(Y2f2/Nf2) + D3X2_m32_reg*(Z2f2/Nf2) + D3K2_m32_reg*(A2f2/Nf2)
    # infected partner on ART
    + D3P3_m32_reg*(H3f2/Nf2) + D3G3_m32_reg*(Y3f2/Nf2) + D3X3_m32_reg*(Z3f2/Nf2) + D3K3_m32_reg*(A3f2/Nf2)
  )
  +(dm33_reg*rho_m3_3_reg)*(
    # infected partner not on HCT/ART
      D3P1_m33_reg*(H1f3/Nf3) + D3G1_m33_reg*(Y1f3/Nf3) + D3X1_m33_reg*(Z1f3/Nf3) + D3K1_m33_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D3P2_m33_reg*(H2f3/Nf3) + D3G2_m33_reg*(Y2f3/Nf3) + D3X2_m33_reg*(Z2f3/Nf3) + D3K2_m33_reg*(A2f3/Nf3)
    # infected partner on ART
    + D3P3_m33_reg*(H3f3/Nf3) + D3G3_m33_reg*(Y3f3/Nf3) + D3X3_m33_reg*(Z3f3/Nf3) + D3K3_m33_reg*(A3f3/Nf3)
  )
  # casual
  +(dm33_casual*rho_m3_3_casual)*(
    # infected partner not on HCT/ART
      D3P1_m33_casual*(H1f3/Nf3) + D3G1_m33_casual*(Y1f3/Nf3) + D3X1_m33_casual*(Z1f3/Nf3) + D3K1_m33_casual*(A1f3/Nf3)
    # infected partner on HCT
    + D3P2_m33_casual*(H2f3/Nf3) + D3G2_m33_casual*(Y2f3/Nf3) + D3X2_m33_casual*(Z2f3/Nf3) + D3K2_m33_casual*(A2f3/Nf3)
    # infected partner on ART
    + D3P3_m33_casual*(H3f3/Nf3) + D3G3_m33_casual*(Y3f3/Nf3) + D3X3_m33_casual*(Z3f3/Nf3) + D3K3_m33_casual*(A3f3/Nf3)
  )
)


###
# Circumcised males on PrEP
pi4_m1 = (
  # low-low
  (dm11_main*rho_m1_1_main)*(
    # infected partner not on HCT/ART
      D4P1_m11_main*(H1f1/Nf1) + D4G1_m11_main*(Y1f1/Nf1) + D4X1_m11_main*(Z1f1/Nf1) + D4K1_m11_main*(A1f1/Nf1)
    # infected partner on HCT
    + D4P2_m11_main*(H2f1/Nf1) + D4G2_m11_main*(Y2f1/Nf1) + D4X2_m11_main*(Z2f1/Nf1) + D4K2_m11_main*(A2f1/Nf1)
    # infected partner on ART
    + D4P3_m11_main*(H3f1/Nf1) + D4G3_m11_main*(Y3f1/Nf1) + D4X3_m11_main*(Z3f1/Nf1) + D4K3_m11_main*(A3f1/Nf1)
  )
  # low-medium
  +(dm12_main*rho_m1_2_main)*(
    # infected partner not on HCT/ART
      D4P1_m12_main*(H1f2/Nf2) + D4G1_m12_main*(Y1f2/Nf2) + D4X1_m12_main*(Z1f2/Nf2) + D4K1_m12_main*(A1f2/Nf2)
    # infected partner on HCT
    + D4P2_m12_main*(H2f2/Nf2) + D4G2_m12_main*(Y2f2/Nf2) + D4X2_m12_main*(Z2f2/Nf2) + D4K2_m12_main*(A2f2/Nf2)
    # infected partner on ART
    + D4P3_m12_main*(H3f2/Nf2) + D4G3_m12_main*(Y3f2/Nf2) + D4X3_m12_main*(Z3f2/Nf2) + D4K3_m12_main*(A3f2/Nf2)
  )
  # low-high
  +(dm13_main*rho_m1_3_main)*(
    # infected partner not on HCT/ART
      D4P1_m13_main*(H1f3/Nf3) + D4G1_m13_main*(Y1f3/Nf3) + D4X1_m13_main*(Z1f3/Nf3) + D4K1_m13_main*(A1f3/Nf3)
    # infected partner on HCT
    + D4P2_m13_main*(H2f3/Nf3) + D4G2_m13_main*(Y2f3/Nf3) + D4X2_m13_main*(Z2f3/Nf3) + D4K2_m13_main*(A2f3/Nf3)
    # infected partner on ART
    + D4P3_m13_main*(H3f3/Nf3) + D4G3_m13_main*(Y3f3/Nf3) + D4X3_m13_main*(Z3f3/Nf3) + D4K3_m13_main*(A3f3/Nf3)
  )
)

pi4_m2 = (
  # main
  (dm21_main*rho_m2_1_main)*(
    # infected partner not on HCT/ART
      D4P1_m21_main*(H1f1/Nf1) + D4G1_m21_main*(Y1f1/Nf1) + D4X1_m21_main*(Z1f1/Nf1) + D4K1_m21_main*(A1f1/Nf1)
    # infected partner on HCT
    + D4P2_m21_main*(H2f1/Nf1) + D4G2_m21_main*(Y2f1/Nf1) + D4X2_m21_main*(Z2f1/Nf1) + D4K2_m21_main*(A2f1/Nf1)
    # infected partner on ART
    + D4P3_m21_main*(H3f1/Nf1) + D4G3_m21_main*(Y3f1/Nf1) + D4X3_m21_main*(Z3f1/Nf1) + D4K3_m21_main*(A3f1/Nf1)
  )
  +(dm22_main*rho_m2_2_main)*(
    # infected partner not on HCT/ART
      D4P1_m22_main*(H1f2/Nf2) + D4G1_m22_main*(Y1f2/Nf2) + D4X1_m22_main*(Z1f2/Nf2) + D4K1_m22_main*(A1f2/Nf2)
    # infected partner on HCT
    + D4P2_m22_main*(H2f2/Nf2) + D4G2_m22_main*(Y2f2/Nf2) + D4X2_m22_main*(Z2f2/Nf2) + D4K2_m22_main*(A2f2/Nf2)
    # infected partner on ART
    + D4P3_m22_main*(H3f2/Nf2) + D4G3_m22_main*(Y3f2/Nf2) + D4X3_m22_main*(Z3f2/Nf2) + D4K3_m22_main*(A3f2/Nf2)
  )
  +(dm23_main*rho_m2_3_main)*(
    # infected partner not on HCT/ART
      D4P1_m23_main*(H1f3/Nf3) + D4G1_m23_main*(Y1f3/Nf3) + D4X1_m23_main*(Z1f3/Nf3) + D4K1_m23_main*(A1f3/Nf3)
    # infected partner on HCT
    + D4P2_m23_main*(H2f3/Nf3) + D4G2_m23_main*(Y2f3/Nf3) + D4X2_m23_main*(Z2f3/Nf3) + D4K2_m23_main*(A2f3/Nf3)
    # infected partner on ART
    + D4P3_m23_main*(H3f3/Nf3) + D4G3_m23_main*(Y3f3/Nf3) + D4X3_m23_main*(Z3f3/Nf3) + D4K3_m23_main*(A3f3/Nf3)
  )
  # regular
  +(dm22_reg*rho_m2_2_reg)*(
    # infected partner not on HCT/ART
      D4P1_m22_reg*(H1f2/Nf2) + D4G1_m22_reg*(Y1f2/Nf2) + D4X1_m22_reg*(Z1f2/Nf2) + D4K1_m22_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D4P2_m22_reg*(H2f2/Nf2) + D4G2_m22_reg*(Y2f2/Nf2) + D4X2_m22_reg*(Z2f2/Nf2) + D4K2_m22_reg*(A2f2/Nf2)
    # infected partner on ART
    + D4P3_m22_reg*(H3f2/Nf2) + D4G3_m22_reg*(Y3f2/Nf2) + D4X3_m22_reg*(Z3f2/Nf2) + D4K3_m22_reg*(A3f2/Nf2)
  )
  +(dm23_reg*rho_m2_3_reg)*(
    # infected partner not on HCT/ART
      D4P1_m23_reg*(H1f3/Nf3) + D4G1_m23_reg*(Y1f3/Nf3) + D4X1_m23_reg*(Z1f3/Nf3) + D4K1_m23_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D4P2_m23_reg*(H2f3/Nf3) + D4G2_m23_reg*(Y2f3/Nf3) + D4X2_m23_reg*(Z2f3/Nf3) + D4K2_m23_reg*(A2f3/Nf3)
    # infected partner on ART
    + D4P3_m23_reg*(H3f3/Nf3) + D4G3_m23_reg*(Y3f3/Nf3) + D4X3_m23_reg*(Z3f3/Nf3) + D4K3_m23_reg*(A3f3/Nf3)
  )
)

pi4_m3 = (
  # main
  (dm31_main*rho_m3_1_main)*(
    # infected partner not on HCT/ART
      D4P1_m31_main*(H1f1/Nf1) + D4G1_m31_main*(Y1f1/Nf1) + D4X1_m31_main*(Z1f1/Nf1) + D4K1_m31_main*(A1f1/Nf1)
    # infected partner on HCT
    + D4P2_m31_main*(H2f1/Nf1) + D4G2_m31_main*(Y2f1/Nf1) + D4X2_m31_main*(Z2f1/Nf1) + D4K2_m31_main*(A2f1/Nf1)
    # infected partner on ART
    + D4P3_m31_main*(H3f1/Nf1) + D4G3_m31_main*(Y3f1/Nf1) + D4X3_m31_main*(Z3f1/Nf1) + D4K3_m31_main*(A3f1/Nf1)
  )
  +(dm32_main*rho_m3_2_main)*(
    # infected partner not on HCT/ART
      D4P1_m32_main*(H1f2/Nf2) + D4G1_m32_main*(Y1f2/Nf2) + D4X1_m32_main*(Z1f2/Nf2) + D4K1_m32_main*(A1f2/Nf2)
    # infected partner on HCT
    + D4P2_m32_main*(H2f2/Nf2) + D4G2_m32_main*(Y2f2/Nf2) + D4X2_m32_main*(Z2f2/Nf2) + D4K2_m32_main*(A2f2/Nf2)
    # infected partner on ART
    + D4P3_m32_main*(H3f2/Nf2) + D4G3_m32_main*(Y3f2/Nf2) + D4X3_m32_main*(Z3f2/Nf2) + D4K3_m32_main*(A3f2/Nf2)
  )
  +(dm33_main*rho_m3_3_main)*(
    # infected partner not on HCT/ART
      D4P1_m33_main*(H1f3/Nf3) + D4G1_m33_main*(Y1f3/Nf3) + D4X1_m33_main*(Z1f3/Nf3) + D4K1_m33_main*(A1f3/Nf3)
    # infected partner on HCT
    + D4P2_m33_main*(H2f3/Nf3) + D4G2_m33_main*(Y2f3/Nf3) + D4X2_m33_main*(Z2f3/Nf3) + D4K2_m33_main*(A2f3/Nf3)
    # infected partner on ART
    + D4P3_m33_main*(H3f3/Nf3) + D4G3_m33_main*(Y3f3/Nf3) + D4X3_m33_main*(Z3f3/Nf3) + D4K3_m33_main*(A3f3/Nf3)
  )
  # regular
  +(dm32_reg*rho_m3_2_reg)*(
    # infected partner not on HCT/ART
      D4P1_m32_reg*(H1f2/Nf2) + D4G1_m32_reg*(Y1f2/Nf2) + D4X1_m32_reg*(Z1f2/Nf2) + D4K1_m32_reg*(A1f2/Nf2)
    # infected partner on HCT
    + D4P2_m32_reg*(H2f2/Nf2) + D4G2_m32_reg*(Y2f2/Nf2) + D4X2_m32_reg*(Z2f2/Nf2) + D4K2_m32_reg*(A2f2/Nf2)
    # infected partner on ART
    + D4P3_m32_reg*(H3f2/Nf2) + D4G3_m32_reg*(Y3f2/Nf2) + D4X3_m32_reg*(Z3f2/Nf2) + D4K3_m32_reg*(A3f2/Nf2)
  )
  +(dm33_reg*rho_m3_3_reg)*(
    # infected partner not on HCT/ART
      D4P1_m33_reg*(H1f3/Nf3) + D4G1_m33_reg*(Y1f3/Nf3) + D4X1_m33_reg*(Z1f3/Nf3) + D4K1_m33_reg*(A1f3/Nf3)
    # infected partner on HCT
    + D4P2_m33_reg*(H2f3/Nf3) + D4G2_m33_reg*(Y2f3/Nf3) + D4X2_m33_reg*(Z2f3/Nf3) + D4K2_m33_reg*(A2f3/Nf3)
    # infected partner on ART
    + D4P3_m33_reg*(H3f3/Nf3) + D4G3_m33_reg*(Y3f3/Nf3) + D4X3_m33_reg*(Z3f3/Nf3) + D4K3_m33_reg*(A3f3/Nf3)
  )
  # casual
  +(dm33_casual*rho_m3_3_casual)*(
    # infected partner not on HCT/ART
      D4P1_m33_casual*(H1f3/Nf3) + D4G1_m33_casual*(Y1f3/Nf3) + D4X1_m33_casual*(Z1f3/Nf3) + D4K1_m33_casual*(A1f3/Nf3)
    # infected partner on HCT
    + D4P2_m33_casual*(H2f3/Nf3) + D4G2_m33_casual*(Y2f3/Nf3) + D4X2_m33_casual*(Z2f3/Nf3) + D4K2_m33_casual*(A2f3/Nf3)
    # infected partner on ART
    + D4P3_m33_casual*(H3f3/Nf3) + D4G3_m33_casual*(Y3f3/Nf3) + D4X3_m33_casual*(Z3f3/Nf3) + D4K3_m33_casual*(A3f3/Nf3)
  )
)


## Incidence per 100 PY
IRf1 = ((pi1_f1*S1f1 + pi2_f1*S2f1)/(S1f1+S2f1)*100)*12
IRf2 = ((pi1_f2*S1f2 + pi2_f2*S2f2)/(S1f2+S2f2)*100)*12
IRf3 = ((pi1_f3*S1f3 + pi2_f3*S2f3)/(S1f3+S2f3)*100)*12
IRm1 = ((pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1)/(S1m1+S2m1+S3m1+S4m1)*100)*12
IRm2 = ((pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2)/(S1m2+S2m2+S3m2+S4m2)*100)*12
IRm3 = ((pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)/(S1m3+S2m3+S3m3+S4m3)*100)*12

person.time = S1f1+S2f1+S1m1+S2m1+S3m1+S4m1
             +S1f2+S2f2+S1m2+S2m2+S3m2+S4m2
             +S1f3+S2f3+S1m3+S2m3+S3m3+S4m3
         
IRtot = ((pi1_f1*S1f1 + pi2_f1*S2f1 + pi1_f2*S1f2 + pi2_f2*S2f2 + pi1_f3*S1f3 + pi2_f3*S2f3
        + pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1 + pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2
        + pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)/person.time*100)*12

## New infections
NItot = (pi1_f1*S1f1 + pi2_f1*S2f1 + pi1_f2*S1f2 + pi2_f2*S2f2 + pi1_f3*S1f3 + pi2_f3*S2f3
        + pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1 + pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2
        + pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)

NIf1 = (pi1_f1*S1f1 + pi2_f1*S2f1)
NIf2 = (pi1_f2*S1f2 + pi2_f2*S2f2)
NIf3 = (pi1_f3*S1f3 + pi2_f3*S2f3)
NIm1 = (pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1)         
NIm2 = (pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2)         
NIm3 = (pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)  

# AIDS-related mortality per 100 PY
AMtot = ((delta1*A1f1 + delta1*A2f1 + delta2*A3f1 + delta1*A1f2 + delta1*A2f2 + delta2*A3f2 + delta1*A1f3 + delta1*A2f3 + delta2*A3f3
        + delta1*A1m1 + delta1*A2m1 + delta2*A3m1 + delta1*A1m2 + delta1*A2m2 + delta2*A3m2 + delta1*A1m3 + delta1*A2m3 + delta2*A3m3)/N*100)*12

AMf1 = ((delta1*A1f1 + delta1*A2f1 + delta2*A3f1)/Nf1*100)*12
AMf2 = ((delta1*A1f2 + delta1*A2f2 + delta2*A3f2)/Nf2*100)*12
AMf3 = ((delta1*A1f3 + delta1*A2f3 + delta2*A3f3)/Nf3*100)*12
AMm1 = ((delta1*A1m1 + delta1*A2m1 + delta2*A3m1)/Nm1*100)*12
AMm2 = ((delta1*A1m2 + delta1*A2m2 + delta2*A3m2)/Nm2*100)*12
AMm3 = ((delta1*A1m3 + delta1*A2m3 + delta2*A3m3)/Nm3*100)*12

# ## Treatment uptake current time
# # start ART in 2004
# if(t<=419) {tau3 = 0; tau4 = 0; gamma = 0}
# if(t<504) {tau1 = 0; tau2 = 0}
# # if(t>552) {tau1=1; tau2=1; tau3=1; tau4=1} # test-and-treat 2015
# 
# 
# tauC = HCT_m*fun.tauT(t,rHCT,pHCT,propaware)*fun.circ(t,pCIRC)
# tauA1 = fun.tauA12(t,tau1)
# tauA2 = fun.tauA12(t,tau2)
# tauA3 = fun.tauA34(t,tau3)
# tauA4 = fun.tauA34(t,tau4)
# tauT = fun.tauT(t,rHCT,pHCT,propaware)

###########
### ODE ###
###########
         
         ### Females       
         dS1f1 = mu*Nf1 + om*S2f1 - (pi1_f1 + tauTS*theta + psi)*S1f1
         dS2f1 = tauTS*theta*S1f1 - (pi2_f1 + om + psi)*S2f1
         dH1f1 = pi1_f1*S1f1 + pi2_f1*S2f1 - (nu1 + tauT + psi)*H1f1
         dH2f1 = tauT*H1f1 + gamma*H3f1 - (nu1 + tauA1 + psi)*H2f1
         dH3f1 = tauA1*H2f1 - (nu2 + gamma + psi)*H3f1
         dY1f1 = nu1*H1f1 - (omega1 + tauT + psi)*Y1f1
         dY2f1 = nu1*H2f1 + tauT*Y1f1 + gamma*Y3f1 - (omega1 + tauA2 + psi)*Y2f1
         dY3f1 = nu2*H3f1 + tauA2*Y2f1 - (omega2 + gamma + psi)*Y3f1
         dZ1f1 = omega1*Y1f1 - (epi1 + tauT + psi)*Z1f1
         dZ2f1 = omega1*Y2f1 + tauT*Z1f1 + gamma*Z3f1 - (epi1 + tauA3 + psi)*Z2f1
         dZ3f1 = omega2*Y3f1 + tauA3*Z2f1 - (epi2 + gamma + psi)*Z3f1
         dA1f1 = epi1*Z1f1 - (delta1 + tauT + psi)*A1f1
         dA2f1 = epi1*Z2f1 + tauT*A1f1 + gamma*A3f1 - (delta1 + tauA4 + psi)*A2f1
         dA3f1 = epi2*Z3f1 + tauA4*A2f1 - (delta2 + gamma + psi)*A3f1
         
         dS1f2 = mu*Nf2 + om*S2f2 - (pi1_f2 + tauTS*theta + psi)*S1f2
         dS2f2 = tauTS*theta*S1f2 - (pi2_f2 + om + psi)*S2f2
         dH1f2 = pi1_f2*S1f2 + pi2_f2*S2f2 - (nu1 + tauT + psi)*H1f2
         dH2f2 = tauT*H1f2 + gamma*H3f2 - (nu1 + tauA1 + psi)*H2f2
         dH3f2 = tauA1*H2f2 - (nu2 + gamma + psi)*H3f2
         dY1f2 = nu1*H1f2 - (omega1 + tauT + psi)*Y1f2
         dY2f2 = nu1*H2f2 + tauT*Y1f2 + gamma*Y3f2 - (omega1 + tauA2 + psi)*Y2f2
         dY3f2 = nu2*H3f2 + tauA2*Y2f2 - (omega2 + gamma + psi)*Y3f2
         dZ1f2 = omega1*Y1f2 - (epi1 + tauT + psi)*Z1f2
         dZ2f2 = omega1*Y2f2 + tauT*Z1f2 + gamma*Z3f2 - (epi1 + tauA3 + psi)*Z2f2
         dZ3f2 = omega2*Y3f2 + tauA3*Z2f2 - (epi2 + gamma + psi)*Z3f2
         dA1f2 = epi1*Z1f2 - (delta1 + tauT + psi)*A1f2
         dA2f2 = epi1*Z2f2 + tauT*A1f2 + gamma*A3f2 - (delta1 + tauA4 + psi)*A2f2
         dA3f2 = epi2*Z3f2 + tauA4*A2f2 - (delta2 + gamma + psi)*A3f2
         
         dS1f3 = mu*Nf3 + om*S2f3 - (pi1_f3 + tauTS*theta + psi)*S1f3
         dS2f3 = tauTS*theta*S1f3 - (pi2_f3 + om + psi)*S2f3
         dH1f3 = pi1_f3*S1f3 + pi2_f3*S2f3 - (nu1 + tauT + psi)*H1f3
         dH2f3 = tauT*H1f3 + gamma*H3f3 - (nu1 + tauA1 + psi)*H2f3
         dH3f3 = tauA1*H2f3 - (nu2 + gamma + psi)*H3f3
         dY1f3 = nu1*H1f3 - (omega1 + tauT + psi)*Y1f3
         dY2f3 = nu1*H2f3 + tauT*Y1f3 + gamma*Y3f3 - (omega1 + tauA2 + psi)*Y2f3
         dY3f3 = nu2*H3f3 + tauA2*Y2f3 - (omega2 + gamma + psi)*Y3f3
         dZ1f3 = omega1*Y1f3 - (epi1 + tauT + psi)*Z1f3
         dZ2f3 = omega1*Y2f3 + tauT*Z1f3 + gamma*Z3f3 - (epi1 + tauA3 + psi)*Z2f3
         dZ3f3 = omega2*Y3f3 + tauA3*Z2f3 - (epi2 + gamma + psi)*Z3f3
         dA1f3 = epi1*Z1f3 - (delta1 + tauT + psi)*A1f3
         dA2f3 = epi1*Z2f3 + tauT*A1f3 + gamma*A3f3 - (delta1 + tauA4 + psi)*A2f3
         dA3f3 = epi2*Z3f3 + tauA4*A2f3 - (delta2 + gamma + psi)*A3f3
         
         ### Males
         
         dS1m1 = (1-prop_circ)*mu*Nm1 + om*S2m1 - (pi1_m1 + tauC + HCT_m*tauTS*theta + psi)*S1m1
         dS2m1 = HCT_m*tauTS*theta*S1m1 - (pi2_m1 + tauC + om + psi)*S2m1
         dS3m1 = prop_circ*mu*Nm1 + tauC*S1m1 + om*S4m1 - (pi3_m1 + HCT_m*tauTS*theta + psi)*S3m1
         dS4m1 = HCT_m*tauTS*theta*S3m1 + tauC*S2m1 - (pi4_m1 + om + psi)*S4m1
         dH1m1 = pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1 - (nu1 + HCT_m*tauT + psi)*H1m1
         dH2m1 = HCT_m*tauT*H1m1 + gamma*H3m1 - (nu1 + ART_m*tauA1 + psi)*H2m1
         dH3m1 = ART_m*tauA1*H2m1 - (nu2 + gamma + psi)*H3m1
         dY1m1 = nu1*H1m1 - (omega1 + HCT_m*tauT + psi)*Y1m1
         dY2m1 = nu1*H2m1 + HCT_m*tauT*Y1m1 + gamma*Y3m1 - (omega1 + ART_m*tauA2 + psi)*Y2m1
         dY3m1 = nu2*H3m1 + ART_m*tauA2*Y2m1 - (omega2 + gamma + psi)*Y3m1
         dZ1m1 = omega1*Y1m1 - (epi1 + HCT_m*tauT + psi)*Z1m1
         dZ2m1 = omega1*Y2m1 + HCT_m*tauT*Z1m1 + gamma*Z3m1 - (epi1 + ART_m*tauA3 + psi)*Z2m1
         dZ3m1 = omega2*Y3m1 + ART_m*tauA3*Z2m1 - (epi2 + gamma + psi)*Z3m1
         dA1m1 = epi1*Z1m1 - (delta1 + HCT_m*tauT + psi)*A1m1
         dA2m1 = epi1*Z2m1 + HCT_m*tauT*A1m1 + gamma*A3m1 - (delta1 + ART_m*tauA4 + psi)*A2m1
         dA3m1 = epi2*Z3m1 + ART_m*tauA4*A2m1 - (delta2 + gamma + psi)*A3m1
   
         
         dS1m2 = (1-prop_circ)*mu*Nm2 + om*S2m2 - (pi1_m2 + tauC + HCT_m*tauTS*theta + psi)*S1m2
         dS2m2 = HCT_m*tauTS*theta*S1m2 - (pi2_m2 + tauC + om + psi)*S2m2
         dS3m2 = prop_circ*mu*Nm2 + tauC*S1m2 + om*S4m2 - (pi3_m2 + HCT_m*tauTS*theta + psi)*S3m2
         dS4m2 = HCT_m*tauTS*theta*S3m2 + tauC*S2m2 - (pi4_m2 + om + psi)*S4m2
         dH1m2 = pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2 - (nu1 + HCT_m*tauT + psi)*H1m2
         dH2m2 = HCT_m*tauT*H1m2 + gamma*H3m2 - (nu1 + ART_m*tauA1 + psi)*H2m2
         dH3m2 = ART_m*tauA1*H2m2 - (nu2 + gamma + psi)*H3m2
         dY1m2 = nu1*H1m2 - (omega1 + HCT_m*tauT + psi)*Y1m2
         dY2m2 = nu1*H2m2 + HCT_m*tauT*Y1m2 + gamma*Y3m2 - (omega1 + ART_m*tauA2 + psi)*Y2m2
         dY3m2 = nu2*H3m2 + ART_m*tauA2*Y2m2 - (omega2 + gamma + psi)*Y3m2
         dZ1m2 = omega1*Y1m2 - (epi1 + HCT_m*tauT + psi)*Z1m2
         dZ2m2 = omega1*Y2m2 + HCT_m*tauT*Z1m2 + gamma*Z3m2 - (epi1 + ART_m*tauA3 + psi)*Z2m2
         dZ3m2 = omega2*Y3m2 + ART_m*tauA3*Z2m2 - (epi2 + gamma + psi)*Z3m2
         dA1m2 = epi1*Z1m2 - (delta1 + HCT_m*tauT + psi)*A1m2
         dA2m2 = epi1*Z2m2 + HCT_m*tauT*A1m2 + gamma*A3m2 - (delta1 + ART_m*tauA4 + psi)*A2m2
         dA3m2 = epi2*Z3m2 + ART_m*tauA4*A2m2 - (delta2 + gamma + psi)*A3m2
         
         dS1m3 = (1-prop_circ)*mu*Nm3 + om*S2m3 - (pi1_m3 + tauC + HCT_m*tauTS*theta + psi)*S1m3
         dS2m3 = HCT_m*tauTS*theta*S1m3 - (pi2_m3 + tauC + om + psi)*S2m3
         dS3m3 = prop_circ*mu*Nm3 + tauC*S1m3 + om*S4m3 - (pi3_m3 + HCT_m*tauTS*theta + psi)*S3m3
         dS4m3 = HCT_m*tauTS*theta*S3m3 + tauC*S2m3 - (pi4_m3 + om + psi)*S4m3
         dH1m3 = pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3 - (nu1 + HCT_m*tauT + psi)*H1m3
         dH2m3 = HCT_m*tauT*H1m3 + gamma*H3m3 - (nu1 + ART_m*tauA1 + psi)*H2m3
         dH3m3 = ART_m*tauA1*H2m3 - (nu2 + gamma + psi)*H3m3
         dY1m3 = nu1*H1m3 - (omega1 + HCT_m*tauT + psi)*Y1m3
         dY2m3 = nu1*H2m3 + HCT_m*tauT*Y1m3 + gamma*Y3m3 - (omega1 + ART_m*tauA2 + psi)*Y2m3
         dY3m3 = nu2*H3m3 + ART_m*tauA2*Y2m3 - (omega2 + gamma + psi)*Y3m3
         dZ1m3 = omega1*Y1m3 - (epi1 + HCT_m*tauT + psi)*Z1m3
         dZ2m3 = omega1*Y2m3 + HCT_m*tauT*Z1m3 + gamma*Z3m3 - (epi1 + ART_m*tauA3 + psi)*Z2m3
         dZ3m3 = omega2*Y3m3 + ART_m*tauA3*Z2m3 - (epi2 + gamma + psi)*Z3m3
         dA1m3 = epi1*Z1m3 - (delta1 + HCT_m*tauT + psi)*A1m3
         dA2m3 = epi1*Z2m3 + HCT_m*tauT*A1m3 + gamma*A3m3 - (delta1 + ART_m*tauA4 + psi)*A2m3
         dA3m3 = epi2*Z3m3 + ART_m*tauA4*A2m3 - (delta2 + gamma + psi)*A3m3
         
## Prevalence
PRf1 = (H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1)/Nf1
PRf2 = (H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2)/Nf2
PRf3 = (H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3)/Nf3
PRm1 = (H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1)/Nm1
PRm2 = (H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2)/Nm2
PRm3 = (H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/Nm3

PRtot = (( H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1
          +H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2
          +H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3
          +H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1
          +H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2
          +H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/N)
         

# Proportion aware of HIV+ status
prop.aware = (H2f1+H3f1+H2f2+H3f2+H2f3+H3f3+H2m1+H3m1+H2m2+H3m2+H2m3+H3m3+
              Y2f1+Y3f1+Y2f2+Y3f2+Y2f3+Y3f3+Y2m1+Y3m1+Y2m2+Y3m2+Y2m3+Y3m3+
              Z2f1+Z3f1+Z2f2+Z3f2+Z2f3+Z3f3+Z2m1+Z3m1+Z2m2+Z3m2+Z2m3+Z3m3+
              A2f1+A3f1+A2f2+A3f2+A2f3+A3f3+A2m1+A3m1+A2m2+A3m2+A2m3+A3m3)/(Htot+Ytot+Ztot+Atot)



         ## output
         list(c(dS1f1,dH1f1,dY1f1,dZ1f1,dA1f1,
                dS1f2,dH1f2,dY1f2,dZ1f2,dA1f2,
                dS1f3,dH1f3,dY1f3,dZ1f3,dA1f3,
                dS1m1,dH1m1,dY1m1,dZ1m1,dA1m1,
                dS1m2,dH1m2,dY1m2,dZ1m2,dA1m2,
                dS1m3,dH1m3,dY1m3,dZ1m3,dA1m3,
                dS2f1,dH2f1,dY2f1,dZ2f1,dA2f1,
                dS2f2,dH2f2,dY2f2,dZ2f2,dA2f2,
                dS2f3,dH2f3,dY2f3,dZ2f3,dA2f3,
                dS2m1,dH2m1,dY2m1,dZ2m1,dA2m1,
                dS2m2,dH2m2,dY2m2,dZ2m2,dA2m2,
                dS2m3,dH2m3,dY2m3,dZ2m3,dA2m3,
                dS3m1,dS3m2,dS3m3,dS4m1,dS4m2,dS4m3,
                dH3f1,dY3f1,dZ3f1,dA3f1,
                dH3f2,dY3f2,dZ3f2,dA3f2,
                dH3f3,dY3f3,dZ3f3,dA3f3,
                dH3m1,dY3m1,dZ3m1,dA3m1,
                dH3m2,dY3m2,dZ3m2,dA3m2,
                dH3m3,dY3m3,dZ3m3,dA3m3),
              IRf1=IRf1,IRf2=IRf2,IRf3=IRf3,IRm1=IRm1,IRm2=IRm2,IRm3=IRm3,
              PRf1=PRf1,PRf2=PRf2,PRf3=PRf3,PRm1=PRm1,PRm2=PRm2,PRm3=PRm3,
              IRtot=IRtot,PRtot=PRtot,AMtot=AMtot,
              AMf1=AMf1, AMf2=AMf2, AMf3=AMf3, AMm1=AMm1, AMm2=AMm2, AMm3=AMm3,
              Stot=Stot,Htot=Htot,Ytot=Ytot,Ztot=Ztot,Atot=Atot,
              NItot=NItot,NIf1=NIf1, NIf2=NIf2, NIf3=NIf3, NIm1=NIm1, NIm2=NIm2, NIm3=NIm3)
       }
  )
}

########################
### SOLVE ODE SYSTEM ###
########################

state=c(S1f1=S1f1,H1f1=H1f1,Y1f1=Y1f1,Z1f1=Z1f1,A1f1=A1f1,
        S1f2=S1f2,H1f2=H1f2,Y1f2=Y1f2,Z1f2=Z1f2,A1f2=A1f2,
        S1f3=S1f3,H1f3=H1f3,Y1f3=Y1f3,Z1f3=Z1f3,A1f3=A1f3,
        S1m1=S1m1,H1m1=H1m1,Y1m1=Y1m1,Z1m1=Z1m1,A1m1=A1m1,
        S1m2=S1m2,H1m2=H1m2,Y1m2=Y1m2,Z1m2=Z1m2,A1m2=A1m2,
        S1m3=S1m3,H1m3=H1m3,Y1m3=Y1m3,Z1m3=Z1m3,A1m3=A1m3,
        S2f1=S2f1,H2f1=H2f1,Y2f1=Y2f1,Z2f1=Z2f1,A2f1=A2f1,
        S2f2=S2f2,H2f2=H2f2,Y2f2=Y2f2,Z2f2=Z2f2,A2f2=A2f2,
        S2f3=S2f3,H2f3=H2f3,Y2f3=Y2f3,Z2f3=Z2f3,A2f3=A2f3,
        S2m1=S2m1,H2m1=H2m1,Y2m1=Y2m1,Z2m1=Z2m1,A2m1=A2m1,
        S2m2=S2m2,H2m2=H2m2,Y2m2=Y2m2,Z2m2=Z2m2,A2m2=A2m2,
        S2m3=S2m3,H2m3=H2m3,Y2m3=Y2m3,Z2m3=Z2m3,A2m3=A2m3,
        S3m1=S3m1,S3m2=S3m2,S3m3=S3m3,S4m1=S4m1,S4m2=S4m2,S4m3=S4m3,
        H3f1=H3f1,Y3f1=Y3f1,Z3f1=Z3f1,A3f1=A3f1,
        H3f2=H3f2,Y3f2=Y3f2,Z3f2=Z3f2,A3f2=A3f2,
        H3f3=H3f3,Y3f3=Y3f3,Z3f3=Z3f3,A3f3=A3f3,
        H3m1=H3m1,Y3m1=Y3m1,Z3m1=Z3m1,A3m1=A3m1,
        H3m2=H3m2,Y3m2=Y3m2,Z3m2=Z3m2,A3m2=A3m2,
        H3m3=H3m3,Y3m3=Y3m3,Z3m3=Z3m3,A3m3=A3m3)
#state=ifelse(state<1e-6,0,state)

#UNIT=months
times=seq(0,720,by=1) #run for 45 years (= 1 generation of sexually active persons)


parameters=c(mu=mu, psi=psi,
             nu1=1/avg.nu1, omega1=1/avg.omega1, epi1=1/avg.epi1, delta1=1/avg.delta1,
             nu2=1/avg.nu2, omega2=1/avg.omega2, epi2=1/avg.epi2, delta2=1/avg.delta2)

require(deSolve) 

trial1_FM_BASE=as.data.frame(ode(y=state,times=times,func=SHYZA_FM,parms=parameters))
tail(trial1_FM_BASE)


######################################################################################################
### INTERVENTION SCALE-UP

## HCT uptake (monthly rate) / pHCT = desired proportion aware
fun.tauT = function(t,rHCT,pHCT,prop, rHCT2, pHCT2){
  tau_t=0
  if(t>217 & t<=(217+1/rHCT)){ # testing since 1987
    tau_t = ((t-217)*rHCT*pHCT) - prop
  } 
  else if(t>(217+1/rHCT)){
    tau_t = tau_t + ((t-540)*rHCT2*(pHCT2-pHCT)) - (prop-pHCT)
  }
  return(tau_t)
}

## HCT for susceptibles (after 2014)
fun.hct = function(t,pmax){
  tau_t = 0
  if(t>540){
    h = 5
    tau_t = pmax*((t-540)^h/((t-540)^h+100^h))
  }
  return(tau_t)
}

## ART uptake (increasing probability of accepting treatment after positive test)
fun.tauA34 = function(t,pART,scaleART){ # pART = max prob of accepting treatment
  tau_t = 0
  if(t>420 & t<=540){
    h = 2
    tau_t = pART*((t-420)^h/((t-420)^h+10^h))
  }
  else if(t>540){
    h = 6
    tau_t = pART + scaleART*((t-540)^h/((t-540)^h+95^h))
  }
  return(tau_t)
}

fun.tauA12 = function(t,pART,scaleART){ # pART = max prob of accepting treatment
  tau_t = 0
  if(t>504 & t<=540){
    h = 2
    tau_t = pART*((t-504)^h/((t-504)^h+10^h))
  }
  else if(t>540){
    h = 6
    tau_t = pART + scaleART*((t-540)^h/((t-540)^h+95^h))
  }
  return(tau_t)
}

## PrEP from 2015 onward
fun.prep = function(t, theta){
  tau_p = 0
  if(t>552){
    h = 5
    tau_p = theta*((t-552)^h/((t-552)^h+85^h))
  }
  return(tau_p)
}

## Circumcision (no 'extra' circumcision before HCT in 1987)
fun.circ = function(t,pCIRC,scaleCIRC){
  tau_t = 0
  if(t>217 & t<=540){
    h = 4
    tau_t = pCIRC*((t-217)^h/((t-217)^h+100^h))
  }
  else if(t>540){
    h = 6
    tau_t = pCIRC + scaleCIRC*((t-540)^h/((t-540)^h+85^h))
  }
  return(tau_t)
}

## Condom use
fun.cond = function(t,cond_max,scaleCOND){
  h=4
  cond = cond_max*((t-0)^h/((t-0)^h+200^h))
  if(t>540){
    h = 6
    cond = cond + scaleCOND*((t-540)^h/((t-540)^h+85^h))
  }
  return(cond)
}


#####################################
library(deSolve)

### INITIAL STATES

init=0.01

N = 2000

Nf1=0.43*N; Nm1=0.33*N
Nf2=0.02*N; Nm2=0.07*N
Nf3=0.05*N ; Nm3=0.1*N 

H2m1=0; H2m2=0; H2m3=0; H2f1=0; H2f2=0; H2f3=0
H3m1=0; H3m2=0; H3m3=0; H3f1=0; H3f2=0; H3f3=0

Y2m1=0; Y2m2=0; Y2m3=0; Y2f1=0; Y2f2=0; Y2f3=0
Y3m1=0; Y3m2=0; Y3m3=0; Y3f1=0; Y3f2=0; Y3f3=0

Z2m1=0; Z2m2=0; Z2m3=0; Z2f1=0; Z2f2=0; Z2f3=0
Z3m1=0; Z3m2=0; Z3m3=0; Z3f1=0; Z3f2=0; Z3f3=0

A2m1=0; A2m2=0; A2m3=0; A2f1=0; A2f2=0; A2f3=0
A3m1=0; A3m2=0; A3m3=0; A3f1=0; A3f2=0; A3f3=0
# 
H1m1=init*Nm1; Y1m1=0; Z1m1=0; A1m1=0
H1f1=init*Nf1; Y1f1=0; Z1f1=0; A1f1=0
H1m2=init*Nm2; Y1m2=0; Z1m2=0; A1m2=0
H1f2=init*Nf2; Y1f2=0; Z1f2=0; A1f2=0
H1m3=init*Nm3; Y1m3=0; Z1m3=0; A1m3=0
H1f3=init*Nf3; Y1f3=0; Z1f3=0; A1f3=0

S2f1=0; S2f2=0; S2f3=0; S2m1=0; S2m2=0; S2m3=0
S3m1=0; S3m2=0; S3m3=0; S4m1=0; S4m2=0; S4m3=0

S1f1=Nf1-S2f1-H1f1-Y1f1-Z1f1-A1f1-H2f1-Y2f1-Z2f1-A2f1-H3f1-Y3f1-Z3f1-A3f1
S1f2=Nf2-S2f2-H1f2-Y1f2-Z1f2-A1f2-H2f2-Y2f2-Z2f2-A2f2-H3f2-Y3f2-Z3f2-A3f2 
S1f3=Nf3-S2f3-H1f3-Y1f3-Z1f3-A1f3-H2f3-Y2f3-Z2f3-A2f3-H3f3-Y3f3-Z3f3-A3f3 

S1m1=Nm1-S2m1-S3m1-S4m1-H1m1-Y1m1-Z1m1-A1m1-H2m1-Y2m1-Z2m1-A2m1-H3m1-Y3m1-Z3m1-A3m1
S1m2=Nm2-S2m2-S3m2-S4m2-H1m2-Y1m2-Z1m2-A1m2-H2m2-Y2m2-Z2m2-A2m2-H3m2-Y3m2-Z3m2-A3m2 
S1m3=Nm3-S2m3-S3m3-S4m3-H1m3-Y1m3-Z1m3-A1m3-H2m3-Y2m3-Z2m3-A2m3-H3m3-Y3m3-Z3m3-A3m3 

########################
### INPUT PARAMETERS ###
########################

## Demographic
psi = 0.0019
mu = 0.0019 + (1.025^(1/12)-1)

## Epidemiologic
avg.nu1 = 3.25; avg.omega1 = 100; avg.epi1 = 14; avg.delta1 = 16
avg.nu2 = 98; avg.omega2 = 648; avg.epi2 = 30; avg.delta2 = 33


## Intervention parameters
ART_eff = 0.96 # reduced infectiousness
cond_eff = 0.1
HCT_cond = 1.1 # 10% increase in condom use (1=no increase) due to HCT
e_circ = 0.6 # relative susceptibility of circumcised men = 0.4 = reduction of 60%
prop_circ = prop_circ # proportion circumcised when entering the population
HCT_red = 0.9 # 10% reduction in number of sex acts (1=no reduction) due to HCT
AIDS_red = 0.1 # 90% reduction in number of sex acts due to untreated AIDS

## Intervention uptake
rHCT = 1/323 # HCT coverage increases over 323 months 
pHCT = pHCT # specified coverage level 2014 (proportion aware of HIV+ status)
# tau1 = (1+(-log(1-tau1)))^(1/12) - 1 # proportion initiating ART on average one year after testing
# tau2 = (1+(-log(1-tau2)))^(1/12) - 1
# tau3 = (1+(-log(1-tau3)))^(1/12) - 1
# tau4 = (1+(-log(1-tau4)))^(1/12) - 1
# pCIRC = (1+(-log(1-pCIRC)))^(1/12) - 1 # pCIRC = proportion accepting VMMC on average one year after testing
# 
# gamma = (1+(-log(1-0.2)))^(1/12) - 1 # monthly rate of ART dropout (20% of individuals no longer on treatment one year after starting)

ART_m = ART_m # ART uptake relative to females
HCT_m = HCT_m # HCT uptake relative to females

## PrEP
# pPREP = pPREP # prop initiating prep on average one year after testing
# # rPREP = 1/168
# pPREP = (1+(-log(1-pPREP)))^(1/12) - 1 # monthly rate of PrEP uptake; theta = proportion initiating PrEP on average one year after testing
# om = (1+(-log(1-0.2)))^(1/12) - 1 # monthly rate of PrEP discontinuation (20% of individuals no longer on PrEP one year after starting)
# om = 0
e_p = 0.67 # relative reduction in susceptibility on PrEP
prepF = 0.6 # prep adherence
prepM = 0.6
# prep efficacy = e_p*prepF
e_pF = e_p*prepF
e_pM = e_p*prepM

##### Intervention scale-up
# scaleART1 = (1+(-log(1-scaleART1)))^(1/12) - 1
# scaleART2 = (1+(-log(1-scaleART2)))^(1/12) - 1
# scaleART3 = (1+(-log(1-scaleART3)))^(1/12) - 1
# scaleART4 = (1+(-log(1-scaleART4)))^(1/12) - 1
# scaleCIRC = (1+(-log(1-scaleCIRC)))^(1/12) - 1 
scaleCOND1 = scaleCOND1
scaleCOND2 = scaleCOND2
scaleCOND3 = scaleCOND2
rHCT2 = 1/180
pHCT2 = pHCT2
# PhctS = (1+(-log(1-PhctS)))^(1/12) - 1 # yearly rate of accepting HCT for susceptibles
# PhctS = (1+(-log(1-0.9)))^(1/12) - 1

## Transmission
high_v1 = 26
high_v2 = 7
high_v3 = 1

## Behavioral
f_main = f_main # condom use
f_reg = f_reg
f_casual = f_casual

## Longer partnership = more sex acts
actsf11_main = actsf11_main
actsm11_main = actsf11_main
actsf12_main = actsf12_main
actsm21_main = actsf12_main
actsf13_main = actsf13_main
actsm31_main = actsf13_main
actsf21_main = actsf21_main
actsm12_main = actsf21_main
actsf22_main = actsf22_main
actsm22_main = actsf22_main
actsf23_main = actsf23_main
actsm32_main = actsf23_main
actsf31_main = actsf31_main
actsm13_main = actsf31_main
actsf32_main = actsf32_main
actsm23_main = actsf32_main
actsf33_main = actsf33_main
actsm33_main = actsf33_main
actsf22_reg = actsf22_reg
actsm22_reg = actsf22_reg
actsf23_reg = actsf23_reg
actsm32_reg = actsf23_reg
actsf32_reg = actsf32_reg
actsm23_reg = actsf32_reg
actsf33_reg = actsf33_reg
actsm33_reg = actsf33_reg
actsf33_casual = actsf33_casual
actsm33_casual = actsf33_casual

## Rate of partner acquisition (= average number of partners per month)
df3_casual = df3_casual
dm3_casual = dm3_casual
df2_reg = df2_reg
dm2_reg = dm2_reg
df3_reg = df3_reg
dm3_reg = dm3_reg
df1_main = df1_main
dm1_main = dm1_main
df2_main = df2_main
dm2_main = dm2_main
df3_main = df3_main
dm3_main = dm3_main

A_main = A_main # degree of assortativity for main partnerships
A_reg = A_reg


##################################
#### Transmission probabilities D

### Transmission parameters (depending on disease stage, PrEP use, ART, and circumcision status)
beta_f = beta_f # baseline (partner in chronic LV stage) assumed same for male/female infected partner
beta_m = beta_m
# Disease stages (1=initial, 2=preA, 3=AIDS)
beta1_f = high_v1*beta_f
beta1_m = high_v1*beta_m
beta2_f = high_v2*beta_f
beta2_m = high_v2*beta_m
beta3_f = 0 # no transmission in untreated AIDS stage
beta3_m = 0
# Infected partner on ART
beta4_f = beta_f*(1-ART_eff)
beta5_f = beta_f*(1-ART_eff)
beta6_f = beta_f*(1-ART_eff)
beta7_f = beta_f*(1-ART_eff)
beta4_m = beta_m*(1-ART_eff)
beta5_m = beta_m*(1-ART_eff)
beta6_m = beta_m*(1-ART_eff)
beta7_m = beta_m*(1-ART_eff)
# Susceptible on PrEP
beta8_f = beta_f*(1-e_pF)
beta9_f = beta1_f*(1-e_pF)
beta10_f = beta2_f*(1-e_pF)
beta11_f = beta3_f*(1-e_pF)
beta8_m = beta_m*(1-e_pM)
beta9_m = beta1_m*(1-e_pM)
beta10_m = beta2_m*(1-e_pM)
beta11_m = beta3_m*(1-e_pM)
# Circumcised susceptible males
beta12_m = beta_m*(1-e_circ)
beta13_m = beta1_m*(1-e_circ)
beta14_m = beta2_m*(1-e_circ)
beta15_m = beta3_m*(1-e_circ)
# Circumcised males on PrEP
beta28_m = beta12_m*(1-e_pM)
beta29_m = beta13_m*(1-e_pM)
beta30_m = beta14_m*(1-e_pM)
beta31_m = beta15_m*(1-e_pM)
# Infected partner on ART and susceptible on PrEP
beta16_f = beta4_f*(1-e_pF)
beta17_f = beta5_f*(1-e_pF)
beta18_f = beta6_f*(1-e_pF)
beta19_f = beta7_f*(1-e_pF)
beta16_m = beta4_m*(1-e_pM)
beta17_m = beta5_m*(1-e_pM)
beta18_m = beta6_m*(1-e_pM)
beta19_m = beta7_m*(1-e_pM)
# Infected partner on ART and susceptible male circumcised
beta20_m = beta12_m*(1-ART_eff)
beta21_m = beta12_m*(1-ART_eff)
beta22_m = beta12_m*(1-ART_eff)
beta23_m = beta12_m*(1-ART_eff)
# Infected partner on ART & susceptible male circumcised + PrEP
beta24_m = beta20_m*(1-e_pM)
beta25_m = beta21_m*(1-e_pM)
beta26_m = beta22_m*(1-e_pM)
beta27_m = beta23_m*(1-e_pM)


#######################
### ODE SYSTEM      ###
#######################

SHYZA_FM=function(t, state, parameters)
{
  
  with(as.list(c(state,parameters)),
       {
         ####################################################
         ### FORCE OF INFECTION --> FREQUENCY-DEPENDENT ! ###
         ####################################################
         
         Nf1 = S1f1 + S2f1 + H1f1 + H2f1 + Y1f1 + Y2f1 + Z1f1 + Z2f1 + A1f1 + A2f1 + H3f1 + Y3f1 + Z3f1 + A3f1
         Nf2 = S1f2 + S2f2 + H1f2 + H2f2 + Y1f2 + Y2f2 + Z1f2 + Z2f2 + A1f2 + A2f2 + H3f2 + Y3f2 + Z3f2 + A3f2
         Nf3 = S1f3 + S2f3 + H1f3 + H2f3 + Y1f3 + Y2f3 + Z1f3 + Z2f3 + A1f3 + A2f3 + H3f3 + Y3f3 + Z3f3 + A3f3
         Nm1 = S1m1 + S2m1 + S3m1 + S4m1 + H1m1 + H2m1 + Y1m1 + Y2m1 + Z1m1 + Z2m1 + A1m1 + A2m1 + H3m1 + Y3m1 + Z3m1 + A3m1
         Nm2 = S1m2 + S2m2 + S3m2 + S4m2 + H1m2 + H2m2 + Y1m2 + Y2m2 + Z1m2 + Z2m2 + A1m2 + A2m2 + H3m2 + Y3m2 + Z3m2 + A3m2
         Nm3 = S1m3 + S2m3 + S3m3 + S4m3 + H1m3 + H2m3 + Y1m3 + Y2m3 + Z1m3 + Z2m3 + A1m3 + A2m3 + H3m3 + Y3m3 + Z3m3 + A3m3
         N = Nf1 + Nf2 + Nf3 + Nm1 + Nm2 + Nm3  
         
         Stot = S1f1 + S2f1 + S1f2 + S2f2 + S1f3 + S2f3 + S1m1 + S2m1 + S3m1 + S4m1 + S1m2 + S2m2 + S3m2 + S4m2 + S1m3 + S2m3 + S3m3 + S4m3
         Htot = H1f1 + H2f1 + H3f1 + H1f2 + H2f2 + H3f2 + H1f3 + H2f3 + H3f3 + H1m1 + H2m1 + H3m1 + H1m2 + H2m2 + H3m2 + H1m3 + H2m3 + H3m3
         Ytot = Y1f1 + Y2f1 + Y3f1 + Y1f2 + Y2f2 + Y3f2 + Y1f3 + Y2f3 + Y3f3 + Y1m1 + Y2m1 + Y3m1 + Y1m2 + Y2m2 + Y3m2 + Y1m3 + Y2m3 + Y3m3
         Ztot = Z1f1 + Z2f1 + Z3f1 + Z1f2 + Z2f2 + Z3f2 + Z1f3 + Z2f3 + Z3f3 + Z1m1 + Z2m1 + Z3m1 + Z1m2 + Z2m2 + Z3m2 + Z1m3 + Z2m3 + Z3m3
         Atot = A1f1 + A2f1 + A3f1 + A1f2 + A2f2 + A3f2 + A1f3 + A2f3 + A3f3 + A1m1 + A2m1 + A3m1 + A1m2 + A2m2 + A3m2 + A1m3 + A2m3 + A3m3
         
         S.male = S1m1 + S2m1 + S3m1 + S4m1 + S1m2 + S2m2 + S3m2 + S4m2 + S1m3 + S2m3 + S3m3 + S4m3
         S.female = S1f1 + S2f1 + S1f2 + S2f2 + S1f3 + S2f3 
         
         
         
         ## Proportion aware of HIV+ status
         propaware = (H2f1+H3f1+H2f2+H3f2+H2f3+H3f3+H2m1+H3m1+H2m2+H3m2+H2m3+H3m3+
                        Y2f1+Y3f1+Y2f2+Y3f2+Y2f3+Y3f3+Y2m1+Y3m1+Y2m2+Y3m2+Y2m3+Y3m3+
                        Z2f1+Z3f1+Z2f2+Z3f2+Z2f3+Z3f3+Z2m1+Z3m1+Z2m2+Z3m2+Z2m3+Z3m3+
                        A2f1+A3f1+A2f2+A3f2+A2f3+A3f3+A2m1+A3m1+A2m2+A3m2+A2m3+A3m3)/(Htot+Ytot+Ztot+Atot)
         
         # Proportion of susceptibles on PrEP
         propPREP = (S2f1+S2f2+S2f3+S2m1+S2m2+S2m3+S4m1+S4m2+S4m3)/Stot
         propPREP.male = (S2m1+S2m2+S2m3+S4m1+S4m2+S4m3)/(S2m1+S2m2+S2m3+S4m1+S4m2+S4m3+S1m1+S1m2+S1m3+S3m1+S3m2+S3m3)
         propPREP.female = (S2f1+S2f2+S2f3)/(S2f1+S2f2+S2f3+S1f1+S1f2+S1f3)

         ## Prevalence
         PRf1 = (H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1)/Nf1
         PRf2 = (H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2)/Nf2
         PRf3 = (H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3)/Nf3
         PRm1 = (H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1)/Nm1
         PRm2 = (H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2)/Nm2
         PRm3 = (H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/Nm3
         
         PRtot = (( H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1
                    +H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2
                    +H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3
                    +H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1
                    +H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2
                    +H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/N)
         
         ## Treatment uptake current time
         # start ART in 2004
         if(t<=419) {tau3 = 0; tau4 = 0; gamma = 0}
         if(t<504) {tau1 = 0; tau2 = 0}
         # Uptake rates
         tauC = fun.circ(t,pCIRC,scaleCIRC-pCIRC)*HCT_m*fun.tauT(t,rHCT,pHCT,propaware,rHCT2,pHCT2)
         if(t>540){tauC = fun.circ(t,pCIRC,scaleCIRC-pCIRC)*HCT_m*fun.hct(t,PhctS)}
         tauA1 = fun.tauA12(t,tau1,scaleART1-tau1)
         tauA2 = fun.tauA12(t,tau2,scaleART2-tau2)
         tauA3 = fun.tauA34(t,tau3,scaleART3-tau3)
         tauA4 = fun.tauA34(t,tau4,scaleART4-tau4)
         tauT = fun.tauT(t,rHCT,pHCT,propaware,rHCT2,pHCT2)
         theta = fun.prep(t, pPREP)
         tauTS = fun.hct(t,PhctS)

         # ART incidence per 100 PY
         ARTinc = ((tauA1*H2f1 + tauA2*Y2f1 + tauA3*Z2f1 + tauA4*A2f1 
                    + tauA1*H2f2 + tauA2*Y2f2 + tauA3*Z2f2 + tauA4*A2f2 
                    + tauA1*H2f3 + tauA2*Y2f3 + tauA3*Z2f3 + tauA4*A2f3
                    + ART_m*tauA1*H2m1 + ART_m*tauA2*Y2m1 + ART_m*tauA3*Z2m1 + ART_m*tauA4*A2m1 
                    + ART_m*tauA1*H2m2 + ART_m*tauA2*Y2m2 + ART_m*tauA3*Z2m2 + ART_m*tauA4*A2m2 
                    + ART_m*tauA1*H2m3 + ART_m*tauA2*Y2m3 + ART_m*tauA3*Z2m3 + ART_m*tauA4*A2m3)/(Htot+Ytot+Ztot+Atot)*100)*12
         
         # ART dropout
         ART.drop = (gamma*H3f1 + gamma*H3f2 + gamma*H3f3 +
                       gamma*Y3f1 + gamma*Y3f2 + gamma*Y3f3 +
                       gamma*Z3f1 + gamma*Z3f2 + gamma*Z3f3 +
                       gamma*A3f1 + gamma*A3f2 + gamma*A3f3 +
                       gamma*H3m1 + gamma*H3m2 + gamma*H3m3 +
                       gamma*Y3m1 + gamma*Y3m2 + gamma*Y3m3 +
                       gamma*Z3m1 + gamma*Z3m2 + gamma*Z3m3 +
                       gamma*A3m1 + gamma*A3m2 + gamma*A3m3)
         
         # PrEP dropout
         PrEP.drop = (om*S2f1 + om*S2f2 + om*S2f3 + om*S2m1 + om*S2m2 + om*S2m3 +
                        om*S4m1 + om*S4m2 + om*S4m3)
         
         ### MIXING PARAMETERS
         
         ## MAIN PARTNERSHIPS
         # FEMALES
         rho_f1_1_main = (1 - A_main)* ((dm1_main*Nm1)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) + A_main 
         rho_f1_2_main = (1 - A_main)* ((dm2_main*Nm2)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
         rho_f1_3_main = (1 - A_main)* ((dm3_main*Nm3)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
         
         rho_f2_1_main = (1 - A_main)* ((dm1_main*Nm1)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
         rho_f2_2_main = (1 - A_main)* ((dm2_main*Nm2)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) + A_main
         rho_f2_3_main = (1 - A_main)* ((dm3_main*Nm3)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
         
         rho_f3_1_main = (1 - A_main)* ((dm1_main*Nm1)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
         rho_f3_2_main = (1 - A_main)* ((dm2_main*Nm2)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) 
         rho_f3_3_main = (1 - A_main)* ((dm3_main*Nm3)/(dm1_main*Nm1+dm2_main*Nm2+dm3_main*Nm3)) + A_main
         
         # MALES
         rho_m1_1_main = (1 - A_main)* ((df1_main*Nf1)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) + A_main 
         rho_m1_2_main = (1 - A_main)* ((df2_main*Nf2)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
         rho_m1_3_main = (1 - A_main)* ((df3_main*Nf3)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
         
         rho_m2_1_main = (1 - A_main)* ((df1_main*Nf1)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
         rho_m2_2_main = (1 - A_main)* ((df2_main*Nf2)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) + A_main
         rho_m2_3_main = (1 - A_main)* ((df3_main*Nf3)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
         
         rho_m3_1_main = (1 - A_main)* ((df1_main*Nf1)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
         rho_m3_2_main = (1 - A_main)* ((df2_main*Nf2)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) 
         rho_m3_3_main = (1 - A_main)* ((df3_main*Nf3)/(df1_main*Nf1+df2_main*Nf2+df3_main*Nf3)) + A_main
         
         ## Balancing the number of sexual partnerships
         B11_main = (dm1_main*rho_m1_1_main*Nm1)/(df1_main*rho_f1_1_main*Nf1)
         const1 = const1
         dm11_main = B11_main^(const1-1) * dm1_main
         df11_main = B11_main^const1 * df1_main
         
         B12_main = (dm2_main*rho_m2_1_main*Nm2)/(df1_main*rho_f1_2_main*Nf1)
         const1 = const1
         dm21_main = B12_main^(const1-1) * dm2_main
         df12_main = B12_main^const1 * df1_main
         
         B13_main = (dm3_main*rho_m3_1_main*Nm3)/(df1_main*rho_f1_3_main*Nf1)
         const1 = const1
         dm31_main = B13_main^(const1-1) * dm3_main
         df13_main = B13_main^const1 * df1_main
         
         B21_main = (dm1_main*rho_m1_2_main*Nm1)/(df2_main*rho_f2_1_main*Nf2)
         const1 = const1
         dm12_main = B21_main^(const1-1) * dm1_main
         df21_main = B21_main^const1 * df2_main
         
         B22_main = (dm2_main*rho_m2_2_main*Nm2)/(df2_main*rho_f2_2_main*Nf2)
         const1 = const1
         dm22_main = B22_main^(const1-1) * dm2_main
         df22_main = B22_main^const1 * df2_main
         
         B23_main = (dm3_main*rho_m3_2_main*Nm3)/(df2_main*rho_f2_3_main*Nf2)
         const1 = const1
         dm32_main = B23_main^(const1-1) * dm3_main
         df23_main = B23_main^const1 * df2_main
         
         B31_main = (dm1_main*rho_m1_3_main*Nm1)/(df3_main*rho_f3_1_main*Nf3)
         const1 = const1
         dm13_main = B31_main^(const1-1) * dm1_main
         df31_main = B31_main^const1 * df3_main
         
         B32_main = (dm2_main*rho_m2_3_main*Nm2)/(df3_main*rho_f3_2_main*Nf3)
         const1 = const1
         dm23_main = B32_main^(const1-1) * dm2_main
         df32_main = B32_main^const1 * df3_main
         
         B33_main = (dm3_main*rho_m3_3_main*Nm3)/(df3_main*rho_f3_3_main*Nf3)
         const1 = const1
         dm33_main = B33_main^(const1-1) * dm3_main
         df33_main = B33_main^const1 * df3_main
         
         ## REGULAR PARTNERSHIPS
         # FEMALES
         rho_f2_2_reg = (1 - A_reg)* ((dm2_reg*Nm2)/(dm2_reg*Nm2+dm3_reg*Nm3)) + A_reg 
         rho_f2_3_reg = (1 - A_reg)* ((dm3_reg*Nm3)/(dm2_reg*Nm2+dm3_reg*Nm3))
         
         rho_f3_2_reg = (1 - A_reg)* ((dm2_reg*Nm2)/(dm2_reg*Nm2+dm3_reg*Nm3))
         rho_f3_3_reg = (1 - A_reg)* ((dm3_reg*Nm3)/(dm2_reg*Nm2+dm3_reg*Nm3)) + A_reg
         
         # MALES
         rho_m2_2_reg = (1 - A_reg)* ((df2_reg*Nf2)/(df2_reg*Nf2+df3_reg*Nf3)) + A_reg 
         rho_m2_3_reg = (1 - A_reg)* ((df3_reg*Nf3)/(df2_reg*Nf2+df3_reg*Nf3))
         
         rho_m3_2_reg = (1 - A_reg)* ((df2_reg*Nf2)/(df2_reg*Nf2+df3_reg*Nf3))
         rho_m3_3_reg = (1 - A_reg)* ((df3_reg*Nf3)/(df2_reg*Nf2+df3_reg*Nf3)) + A_reg
         
         
         B22_reg = (dm2_reg*rho_m2_2_reg*Nm2)/(df2_reg*rho_f2_2_reg*Nf2)
         const2 = const2
         dm22_reg = B22_reg^(const2-1) * dm2_reg
         df22_reg = B22_reg^const2 * df2_reg
         
         B23_reg = (dm3_reg*rho_m3_2_reg*Nm3)/(df2_reg*rho_f2_3_reg*Nf2)
         const2 = const2
         dm32_reg = B23_reg^(const2-1) * dm3_reg
         df23_reg = B23_reg^const2 * df2_reg
         
         B32_reg = (dm2_reg*rho_m2_3_reg*Nm2)/(df3_reg*rho_f3_2_reg*Nf3)
         const2 = const2
         dm23_reg = B32_reg^(const2-1) * dm2_reg
         df32_reg = B32_reg^const2 * df3_reg
         
         B33_reg = (dm3_reg*rho_m3_3_reg*Nm3)/(df3_reg*rho_f3_3_reg*Nf3)
         const2 = const2
         dm33_reg = B33_reg^(const2-1) * dm3_reg
         df33_reg = B33_reg^const2 * df3_reg
         
         
         ## CASUAL PARTNERSHIPS (always 1)
         # FEMALES
         rho_f3_3_casual = 1
         
         # MALES
         rho_m3_3_casual = 1
         
         B33_casual = (dm3_casual*rho_m3_3_casual*Nm3)/(df3_casual*rho_f3_3_casual*Nf3)
         const3 = const3
         dm33_casual = B33_casual^(const3-1) * dm3_casual
         df33_casual = B33_casual^const3 * df3_casual
         
         
         ##### PROBABILITY OF TRANSMISSION
         
         ## Condom use at current time
         cond_main = fun.cond(t,f_main,scaleCOND1-f_main)
         cond_reg = fun.cond(t,f_reg,scaleCOND2-f_reg)
         cond_casual = fun.cond(t,f_casual,scaleCOND3-f_casual)
         
         ### Susceptibles not on PrEP/not circumcised ###
         
         ## Partner in LV stage, not on HCT/ART
         
         D1G1_f11_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf11_main) * (1-beta_f)^((1-cond_main)*actsf11_main))
         D1G1_f12_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf12_main) * (1-beta_f)^((1-cond_main)*actsf12_main))
         D1G1_f13_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf13_main) * (1-beta_f)^((1-cond_main)*actsf13_main))
         D1G1_f21_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf21_main) * (1-beta_f)^((1-cond_main)*actsf21_main))
         D1G1_f22_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf22_main) * (1-beta_f)^((1-cond_main)*actsf22_main))
         D1G1_f23_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf23_main) * (1-beta_f)^((1-cond_main)*actsf23_main))
         D1G1_f31_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf31_main) * (1-beta_f)^((1-cond_main)*actsf31_main))
         D1G1_f32_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf32_main) * (1-beta_f)^((1-cond_main)*actsf32_main))
         D1G1_f33_main = 1 - ((1-beta_f*cond_eff)^(cond_main*actsf33_main) * (1-beta_f)^((1-cond_main)*actsf33_main))
         D1G1_m11_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm11_main) * (1-beta_m)^((1-cond_main)*actsm11_main))
         D1G1_m12_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm12_main) * (1-beta_m)^((1-cond_main)*actsm12_main))
         D1G1_m13_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm13_main) * (1-beta_m)^((1-cond_main)*actsm13_main))
         D1G1_m21_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm21_main) * (1-beta_m)^((1-cond_main)*actsm21_main))
         D1G1_m22_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm22_main) * (1-beta_m)^((1-cond_main)*actsm22_main))
         D1G1_m23_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm23_main) * (1-beta_m)^((1-cond_main)*actsm23_main))
         D1G1_m31_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm31_main) * (1-beta_m)^((1-cond_main)*actsm31_main))
         D1G1_m32_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm32_main) * (1-beta_m)^((1-cond_main)*actsm32_main))
         D1G1_m33_main = 1 - ((1-beta_m*cond_eff)^(cond_main*actsm33_main) * (1-beta_m)^((1-cond_main)*actsm33_main))
         D1G1_f22_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta_f)^((1-cond_reg)*actsf22_reg))
         D1G1_f23_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta_f)^((1-cond_reg)*actsf23_reg))
         D1G1_f32_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta_f)^((1-cond_reg)*actsf32_reg))
         D1G1_f33_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta_f)^((1-cond_reg)*actsf33_reg))
         D1G1_m22_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta_m)^((1-cond_reg)*actsm22_reg))
         D1G1_m23_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta_m)^((1-cond_reg)*actsm23_reg))
         D1G1_m32_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta_m)^((1-cond_reg)*actsm32_reg))
         D1G1_m33_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta_m)^((1-cond_reg)*actsm33_reg))
         D1G1_f33_casual = 1 - ((1-beta_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta_f)^((1-cond_casual)*actsf33_casual))
         D1G1_m33_casual = 1 - ((1-beta_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in LV stage, on HCT
         
         D1G2_f11_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf11_main))
         D1G2_f12_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf12_main))
         D1G2_f13_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf13_main))
         D1G2_f21_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf21_main))
         D1G2_f22_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf22_main))
         D1G2_f23_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf23_main))
         D1G2_f31_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf31_main))
         D1G2_f32_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf32_main))
         D1G2_f33_main = 1 - ((1-beta_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta_f)^((1-cond_main*HCT_cond)*HCT_red*actsf33_main))
         D1G2_m11_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
         D1G2_m12_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
         D1G2_m13_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
         D1G2_m21_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
         D1G2_m22_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
         D1G2_m23_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
         D1G2_m31_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
         D1G2_m32_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
         D1G2_m33_main = 1 - ((1-beta_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
         D1G2_f22_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf22_reg))
         D1G2_f23_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf23_reg))
         D1G2_f32_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf32_reg))
         D1G2_f33_reg = 1 - ((1-beta_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf33_reg))
         D1G2_m22_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
         D1G2_m23_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
         D1G2_m32_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
         D1G2_m33_reg = 1 - ((1-beta_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
         D1G2_f33_casual = 1 - ((1-beta_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta_f)^((1-cond_casual*HCT_cond)*HCT_red*actsf33_casual))
         D1G2_m33_casual = 1 - ((1-beta_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))
         
         ## Partner in LV stage, on ART
         
         D1G3_f11_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf11_main) * (1-beta4_f)^((1-cond_main)*actsf11_main))
         D1G3_f12_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf12_main) * (1-beta4_f)^((1-cond_main)*actsf12_main))
         D1G3_f13_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf13_main) * (1-beta4_f)^((1-cond_main)*actsf13_main))
         D1G3_f21_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf21_main) * (1-beta4_f)^((1-cond_main)*actsf21_main))
         D1G3_f22_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf22_main) * (1-beta4_f)^((1-cond_main)*actsf22_main))
         D1G3_f23_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf23_main) * (1-beta4_f)^((1-cond_main)*actsf23_main))
         D1G3_f31_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf31_main) * (1-beta4_f)^((1-cond_main)*actsf31_main))
         D1G3_f32_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf32_main) * (1-beta4_f)^((1-cond_main)*actsf32_main))
         D1G3_f33_main = 1 - ((1-beta4_f*cond_eff)^(cond_main*actsf33_main) * (1-beta4_f)^((1-cond_main)*actsf33_main))
         D1G3_m11_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm11_main) * (1-beta4_m)^((1-cond_main)*actsm11_main))
         D1G3_m12_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm12_main) * (1-beta4_m)^((1-cond_main)*actsm12_main))
         D1G3_m13_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm13_main) * (1-beta4_m)^((1-cond_main)*actsm13_main))
         D1G3_m21_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm21_main) * (1-beta4_m)^((1-cond_main)*actsm21_main))
         D1G3_m22_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm22_main) * (1-beta4_m)^((1-cond_main)*actsm22_main))
         D1G3_m23_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm23_main) * (1-beta4_m)^((1-cond_main)*actsm23_main))
         D1G3_m31_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm31_main) * (1-beta4_m)^((1-cond_main)*actsm31_main))
         D1G3_m32_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm32_main) * (1-beta4_m)^((1-cond_main)*actsm32_main))
         D1G3_m33_main = 1 - ((1-beta4_m*cond_eff)^(cond_main*actsm33_main) * (1-beta4_m)^((1-cond_main)*actsm33_main))
         D1G3_f22_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta4_f)^((1-cond_reg)*actsf22_reg))
         D1G3_f23_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta4_f)^((1-cond_reg)*actsf23_reg))
         D1G3_f32_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta4_f)^((1-cond_reg)*actsf32_reg))
         D1G3_f33_reg = 1 - ((1-beta4_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta4_f)^((1-cond_reg)*actsf33_reg))
         D1G3_m22_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta4_m)^((1-cond_reg)*actsm22_reg))
         D1G3_m23_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta4_m)^((1-cond_reg)*actsm23_reg))
         D1G3_m32_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta4_m)^((1-cond_reg)*actsm32_reg))
         D1G3_m33_reg = 1 - ((1-beta4_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta4_m)^((1-cond_reg)*actsm33_reg))
         D1G3_f33_casual = 1 - ((1-beta4_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta4_f)^((1-cond_casual)*actsf33_casual))
         D1G3_m33_casual = 1 - ((1-beta4_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta4_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in initial HV stage, not on HCT/ART
         
         D1P1_f11_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf11_main) * (1-beta1_f)^((1-cond_main)*actsf11_main))
         D1P1_f12_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf12_main) * (1-beta1_f)^((1-cond_main)*actsf12_main))
         D1P1_f13_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf13_main) * (1-beta1_f)^((1-cond_main)*actsf13_main))
         D1P1_f21_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf21_main) * (1-beta1_f)^((1-cond_main)*actsf21_main))
         D1P1_f22_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf22_main) * (1-beta1_f)^((1-cond_main)*actsf22_main))
         D1P1_f23_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf23_main) * (1-beta1_f)^((1-cond_main)*actsf23_main))
         D1P1_f31_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf31_main) * (1-beta1_f)^((1-cond_main)*actsf31_main))
         D1P1_f32_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf32_main) * (1-beta1_f)^((1-cond_main)*actsf32_main))
         D1P1_f33_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*actsf33_main) * (1-beta1_f)^((1-cond_main)*actsf33_main))
         D1P1_m11_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm11_main) * (1-beta1_m)^((1-cond_main)*actsm11_main))
         D1P1_m12_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm12_main) * (1-beta1_m)^((1-cond_main)*actsm12_main))
         D1P1_m13_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm13_main) * (1-beta1_m)^((1-cond_main)*actsm13_main))
         D1P1_m21_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm21_main) * (1-beta1_m)^((1-cond_main)*actsm21_main))
         D1P1_m22_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm22_main) * (1-beta1_m)^((1-cond_main)*actsm22_main))
         D1P1_m23_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm23_main) * (1-beta1_m)^((1-cond_main)*actsm23_main))
         D1P1_m31_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm31_main) * (1-beta1_m)^((1-cond_main)*actsm31_main))
         D1P1_m32_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm32_main) * (1-beta1_m)^((1-cond_main)*actsm32_main))
         D1P1_m33_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*actsm33_main) * (1-beta1_m)^((1-cond_main)*actsm33_main))
         D1P1_f22_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta1_f)^((1-cond_reg)*actsf22_reg))
         D1P1_f23_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta1_f)^((1-cond_reg)*actsf23_reg))
         D1P1_f32_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta1_f)^((1-cond_reg)*actsf32_reg))
         D1P1_f33_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta1_f)^((1-cond_reg)*actsf33_reg))
         D1P1_m22_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta1_m)^((1-cond_reg)*actsm22_reg))
         D1P1_m23_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta1_m)^((1-cond_reg)*actsm23_reg))
         D1P1_m32_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta1_m)^((1-cond_reg)*actsm32_reg))
         D1P1_m33_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta1_m)^((1-cond_reg)*actsm33_reg))
         D1P1_f33_casual = 1 - ((1-beta1_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta1_f)^((1-cond_casual)*actsf33_casual))
         D1P1_m33_casual = 1 - ((1-beta1_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta1_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in initial HV stage, on HCT
         
         D1P2_f11_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf11_main))
         D1P2_f12_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf12_main))
         D1P2_f13_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf13_main))
         D1P2_f21_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf21_main))
         D1P2_f22_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf22_main))
         D1P2_f23_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf23_main))
         D1P2_f31_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf31_main))
         D1P2_f32_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf32_main))
         D1P2_f33_main = 1 - ((1-beta1_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta1_f)^((1-cond_main*HCT_cond)*HCT_red*actsf33_main))
         D1P2_m11_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
         D1P2_m12_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
         D1P2_m13_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
         D1P2_m21_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
         D1P2_m22_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
         D1P2_m23_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
         D1P2_m31_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
         D1P2_m32_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
         D1P2_m33_main = 1 - ((1-beta1_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta1_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
         D1P2_f22_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf22_reg))
         D1P2_f23_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf23_reg))
         D1P2_f32_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf32_reg))
         D1P2_f33_reg = 1 - ((1-beta1_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta1_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf33_reg))
         D1P2_m22_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
         D1P2_m23_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
         D1P2_m32_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
         D1P2_m33_reg = 1 - ((1-beta1_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta1_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
         D1P2_f33_casual = 1 - ((1-beta1_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta1_f)^((1-cond_casual*HCT_cond)*HCT_red*actsf33_casual))
         D1P2_m33_casual = 1 - ((1-beta1_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta1_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))
         
         ## Partner in initial HV stage, on ART
         
         D1P3_f11_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf11_main) * (1-beta5_f)^((1-cond_main)*actsf11_main))
         D1P3_f12_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf12_main) * (1-beta5_f)^((1-cond_main)*actsf12_main))
         D1P3_f13_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf13_main) * (1-beta5_f)^((1-cond_main)*actsf13_main))
         D1P3_f21_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf21_main) * (1-beta5_f)^((1-cond_main)*actsf21_main))
         D1P3_f22_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf22_main) * (1-beta5_f)^((1-cond_main)*actsf22_main))
         D1P3_f23_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf23_main) * (1-beta5_f)^((1-cond_main)*actsf23_main))
         D1P3_f31_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf31_main) * (1-beta5_f)^((1-cond_main)*actsf31_main))
         D1P3_f32_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf32_main) * (1-beta5_f)^((1-cond_main)*actsf32_main))
         D1P3_f33_main = 1 - ((1-beta5_f*cond_eff)^(cond_main*actsf33_main) * (1-beta5_f)^((1-cond_main)*actsf33_main))
         D1P3_m11_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm11_main) * (1-beta5_m)^((1-cond_main)*actsm11_main))
         D1P3_m12_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm12_main) * (1-beta5_m)^((1-cond_main)*actsm12_main))
         D1P3_m13_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm13_main) * (1-beta5_m)^((1-cond_main)*actsm13_main))
         D1P3_m21_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm21_main) * (1-beta5_m)^((1-cond_main)*actsm21_main))
         D1P3_m22_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm22_main) * (1-beta5_m)^((1-cond_main)*actsm22_main))
         D1P3_m23_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm23_main) * (1-beta5_m)^((1-cond_main)*actsm23_main))
         D1P3_m31_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm31_main) * (1-beta5_m)^((1-cond_main)*actsm31_main))
         D1P3_m32_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm32_main) * (1-beta5_m)^((1-cond_main)*actsm32_main))
         D1P3_m33_main = 1 - ((1-beta5_m*cond_eff)^(cond_main*actsm33_main) * (1-beta5_m)^((1-cond_main)*actsm33_main))
         D1P3_f22_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta5_f)^((1-cond_reg)*actsf22_reg))
         D1P3_f23_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta5_f)^((1-cond_reg)*actsf23_reg))
         D1P3_f32_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta5_f)^((1-cond_reg)*actsf32_reg))
         D1P3_f33_reg = 1 - ((1-beta5_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta5_f)^((1-cond_reg)*actsf33_reg))
         D1P3_m22_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta5_m)^((1-cond_reg)*actsm22_reg))
         D1P3_m23_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta5_m)^((1-cond_reg)*actsm23_reg))
         D1P3_m32_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta5_m)^((1-cond_reg)*actsm32_reg))
         D1P3_m33_reg = 1 - ((1-beta5_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta5_m)^((1-cond_reg)*actsm33_reg))
         D1P3_f33_casual = 1 - ((1-beta5_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta5_f)^((1-cond_casual)*actsf33_casual))
         D1P3_m33_casual = 1 - ((1-beta5_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta5_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in preAIDS HV stage, not on HCT/ART
         
         D1X1_f11_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf11_main) * (1-beta2_f)^((1-cond_main)*actsf11_main))
         D1X1_f12_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf12_main) * (1-beta2_f)^((1-cond_main)*actsf12_main))
         D1X1_f13_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf13_main) * (1-beta2_f)^((1-cond_main)*actsf13_main))
         D1X1_f21_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf21_main) * (1-beta2_f)^((1-cond_main)*actsf21_main))
         D1X1_f22_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf22_main) * (1-beta2_f)^((1-cond_main)*actsf22_main))
         D1X1_f23_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf23_main) * (1-beta2_f)^((1-cond_main)*actsf23_main))
         D1X1_f31_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf31_main) * (1-beta2_f)^((1-cond_main)*actsf31_main))
         D1X1_f32_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf32_main) * (1-beta2_f)^((1-cond_main)*actsf32_main))
         D1X1_f33_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*actsf33_main) * (1-beta2_f)^((1-cond_main)*actsf33_main))
         D1X1_m11_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm11_main) * (1-beta2_m)^((1-cond_main)*actsm11_main))
         D1X1_m12_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm12_main) * (1-beta2_m)^((1-cond_main)*actsm12_main))
         D1X1_m13_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm13_main) * (1-beta2_m)^((1-cond_main)*actsm13_main))
         D1X1_m21_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm21_main) * (1-beta2_m)^((1-cond_main)*actsm21_main))
         D1X1_m22_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm22_main) * (1-beta2_m)^((1-cond_main)*actsm22_main))
         D1X1_m23_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm23_main) * (1-beta2_m)^((1-cond_main)*actsm23_main))
         D1X1_m31_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm31_main) * (1-beta2_m)^((1-cond_main)*actsm31_main))
         D1X1_m32_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm32_main) * (1-beta2_m)^((1-cond_main)*actsm32_main))
         D1X1_m33_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*actsm33_main) * (1-beta2_m)^((1-cond_main)*actsm33_main))
         D1X1_f22_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta2_f)^((1-cond_reg)*actsf22_reg))
         D1X1_f23_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta2_f)^((1-cond_reg)*actsf23_reg))
         D1X1_f32_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta2_f)^((1-cond_reg)*actsf32_reg))
         D1X1_f33_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta2_f)^((1-cond_reg)*actsf33_reg))
         D1X1_m22_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta2_m)^((1-cond_reg)*actsm22_reg))
         D1X1_m23_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta2_m)^((1-cond_reg)*actsm23_reg))
         D1X1_m32_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta2_m)^((1-cond_reg)*actsm32_reg))
         D1X1_m33_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta2_m)^((1-cond_reg)*actsm33_reg))
         D1X1_f33_casual = 1 - ((1-beta2_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta2_f)^((1-cond_casual)*actsf33_casual))
         D1X1_m33_casual = 1 - ((1-beta2_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta2_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on HCT
         
         D1X2_f11_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf11_main))
         D1X2_f12_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf12_main))
         D1X2_f13_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf13_main))
         D1X2_f21_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf21_main))
         D1X2_f22_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf22_main))
         D1X2_f23_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf23_main))
         D1X2_f31_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf31_main))
         D1X2_f32_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf32_main))
         D1X2_f33_main = 1 - ((1-beta2_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta2_f)^((1-cond_main*HCT_cond)*HCT_red*actsf33_main))
         D1X2_m11_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
         D1X2_m12_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
         D1X2_m13_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
         D1X2_m21_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
         D1X2_m22_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
         D1X2_m23_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
         D1X2_m31_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
         D1X2_m32_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
         D1X2_m33_main = 1 - ((1-beta2_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta2_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
         D1X2_f22_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf22_reg))
         D1X2_f23_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf23_reg))
         D1X2_f32_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf32_reg))
         D1X2_f33_reg = 1 - ((1-beta2_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta2_f)^((1-cond_reg*HCT_cond)*HCT_red*actsf33_reg))
         D1X2_m22_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
         D1X2_m23_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
         D1X2_m32_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
         D1X2_m33_reg = 1 - ((1-beta2_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta2_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
         D1X2_f33_casual = 1 - ((1-beta2_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta2_f)^((1-cond_casual*HCT_cond)*HCT_red*actsf33_casual))
         D1X2_m33_casual = 1 - ((1-beta2_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta2_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on ART
         
         D1X3_f11_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf11_main) * (1-beta6_f)^((1-cond_main)*actsf11_main))
         D1X3_f12_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf12_main) * (1-beta6_f)^((1-cond_main)*actsf12_main))
         D1X3_f13_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf13_main) * (1-beta6_f)^((1-cond_main)*actsf13_main))
         D1X3_f21_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf21_main) * (1-beta6_f)^((1-cond_main)*actsf21_main))
         D1X3_f22_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf22_main) * (1-beta6_f)^((1-cond_main)*actsf22_main))
         D1X3_f23_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf23_main) * (1-beta6_f)^((1-cond_main)*actsf23_main))
         D1X3_f31_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf31_main) * (1-beta6_f)^((1-cond_main)*actsf31_main))
         D1X3_f32_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf32_main) * (1-beta6_f)^((1-cond_main)*actsf32_main))
         D1X3_f33_main = 1 - ((1-beta6_f*cond_eff)^(cond_main*actsf33_main) * (1-beta6_f)^((1-cond_main)*actsf33_main))
         D1X3_m11_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm11_main) * (1-beta6_m)^((1-cond_main)*actsm11_main))
         D1X3_m12_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm12_main) * (1-beta6_m)^((1-cond_main)*actsm12_main))
         D1X3_m13_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm13_main) * (1-beta6_m)^((1-cond_main)*actsm13_main))
         D1X3_m21_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm21_main) * (1-beta6_m)^((1-cond_main)*actsm21_main))
         D1X3_m22_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm22_main) * (1-beta6_m)^((1-cond_main)*actsm22_main))
         D1X3_m23_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm23_main) * (1-beta6_m)^((1-cond_main)*actsm23_main))
         D1X3_m31_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm31_main) * (1-beta6_m)^((1-cond_main)*actsm31_main))
         D1X3_m32_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm32_main) * (1-beta6_m)^((1-cond_main)*actsm32_main))
         D1X3_m33_main = 1 - ((1-beta6_m*cond_eff)^(cond_main*actsm33_main) * (1-beta6_m)^((1-cond_main)*actsm33_main))
         D1X3_f22_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta6_f)^((1-cond_reg)*actsf22_reg))
         D1X3_f23_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta6_f)^((1-cond_reg)*actsf23_reg))
         D1X3_f32_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta6_f)^((1-cond_reg)*actsf32_reg))
         D1X3_f33_reg = 1 - ((1-beta6_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta6_f)^((1-cond_reg)*actsf33_reg))
         D1X3_m22_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta6_m)^((1-cond_reg)*actsm22_reg))
         D1X3_m23_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta6_m)^((1-cond_reg)*actsm23_reg))
         D1X3_m32_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta6_m)^((1-cond_reg)*actsm32_reg))
         D1X3_m33_reg = 1 - ((1-beta6_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta6_m)^((1-cond_reg)*actsm33_reg))
         D1X3_f33_casual = 1 - ((1-beta6_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta6_f)^((1-cond_casual)*actsf33_casual))
         D1X3_m33_casual = 1 - ((1-beta6_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta6_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in AIDS stage, not on HCT/ART
         
         D1K1_f11_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf11_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf11_main))
         D1K1_f12_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf12_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf12_main))
         D1K1_f13_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf13_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf13_main))
         D1K1_f21_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf21_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf21_main))
         D1K1_f22_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf22_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf22_main))
         D1K1_f23_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf23_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf23_main))
         D1K1_f31_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf31_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf31_main))
         D1K1_f32_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf32_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf32_main))
         D1K1_f33_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*AIDS_red*actsf33_main) * (1-beta3_f)^((1-cond_main)*AIDS_red*actsf33_main))
         D1K1_m11_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm11_main))
         D1K1_m12_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm12_main))
         D1K1_m13_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm13_main))
         D1K1_m21_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm21_main))
         D1K1_m22_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm22_main))
         D1K1_m23_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm23_main))
         D1K1_m31_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm31_main))
         D1K1_m32_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm32_main))
         D1K1_m33_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta3_m)^((1-cond_main)*AIDS_red*actsm33_main))
         D1K1_f22_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf22_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf22_reg))
         D1K1_f23_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf23_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf23_reg))
         D1K1_f32_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf32_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf32_reg))
         D1K1_f33_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*AIDS_red*actsf33_reg) * (1-beta3_f)^((1-cond_reg)*AIDS_red*actsf33_reg))
         D1K1_m22_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
         D1K1_m23_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
         D1K1_m32_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
         D1K1_m33_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta3_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
         D1K1_f33_casual = 1 - ((1-beta3_f*cond_eff)^(cond_casual*AIDS_red*actsf33_casual) * (1-beta3_f)^((1-cond_casual)*AIDS_red*actsf33_casual))
         D1K1_m33_casual = 1 - ((1-beta3_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta3_m)^((1-cond_casual)*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on HCT
         
         D1K2_f11_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf11_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf11_main))
         D1K2_f12_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf12_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf12_main))
         D1K2_f13_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf13_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf13_main))
         D1K2_f21_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf21_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf21_main))
         D1K2_f22_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf22_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf22_main))
         D1K2_f23_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf23_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf23_main))
         D1K2_f31_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf31_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf31_main))
         D1K2_f32_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf32_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf32_main))
         D1K2_f33_main = 1 - ((1-beta3_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf33_main) * (1-beta3_f)^((1-cond_main*HCT_cond)*AIDS_red*actsf33_main))
         D1K2_m11_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm11_main))
         D1K2_m12_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm12_main))
         D1K2_m13_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm13_main))
         D1K2_m21_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm21_main))
         D1K2_m22_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm22_main))
         D1K2_m23_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm23_main))
         D1K2_m31_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm31_main))
         D1K2_m32_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm32_main))
         D1K2_m33_main = 1 - ((1-beta3_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta3_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm33_main))
         D1K2_f22_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf22_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf22_reg))
         D1K2_f23_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf23_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf23_reg))
         D1K2_f32_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf32_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf32_reg))
         D1K2_f33_reg = 1 - ((1-beta3_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf33_reg) * (1-beta3_f)^((1-cond_reg*HCT_cond)*AIDS_red*actsf33_reg))
         D1K2_m22_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm22_reg))
         D1K2_m23_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm23_reg))
         D1K2_m32_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm32_reg))
         D1K2_m33_reg = 1 - ((1-beta3_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta3_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm33_reg))
         D1K2_f33_casual = 1 - ((1-beta3_f*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsf33_casual) * (1-beta3_f)^((1-cond_casual*HCT_cond)*AIDS_red*actsf33_casual))
         D1K2_m33_casual = 1 - ((1-beta3_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta3_m)^((1-cond_casual*HCT_cond)*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on ART
         
         D1K3_f11_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf11_main) * (1-beta7_f)^((1-cond_main)*actsf11_main))
         D1K3_f12_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf12_main) * (1-beta7_f)^((1-cond_main)*actsf12_main))
         D1K3_f13_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf13_main) * (1-beta7_f)^((1-cond_main)*actsf13_main))
         D1K3_f21_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf21_main) * (1-beta7_f)^((1-cond_main)*actsf21_main))
         D1K3_f22_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf22_main) * (1-beta7_f)^((1-cond_main)*actsf22_main))
         D1K3_f23_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf23_main) * (1-beta7_f)^((1-cond_main)*actsf23_main))
         D1K3_f31_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf31_main) * (1-beta7_f)^((1-cond_main)*actsf31_main))
         D1K3_f32_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf32_main) * (1-beta7_f)^((1-cond_main)*actsf32_main))
         D1K3_f33_main = 1 - ((1-beta7_f*cond_eff)^(cond_main*actsf33_main) * (1-beta7_f)^((1-cond_main)*actsf33_main))
         D1K3_m11_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm11_main) * (1-beta7_m)^((1-cond_main)*actsm11_main))
         D1K3_m12_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm12_main) * (1-beta7_m)^((1-cond_main)*actsm12_main))
         D1K3_m13_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm13_main) * (1-beta7_m)^((1-cond_main)*actsm13_main))
         D1K3_m21_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm21_main) * (1-beta7_m)^((1-cond_main)*actsm21_main))
         D1K3_m22_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm22_main) * (1-beta7_m)^((1-cond_main)*actsm22_main))
         D1K3_m23_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm23_main) * (1-beta7_m)^((1-cond_main)*actsm23_main))
         D1K3_m31_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm31_main) * (1-beta7_m)^((1-cond_main)*actsm31_main))
         D1K3_m32_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm32_main) * (1-beta7_m)^((1-cond_main)*actsm32_main))
         D1K3_m33_main = 1 - ((1-beta7_m*cond_eff)^(cond_main*actsm33_main) * (1-beta7_m)^((1-cond_main)*actsm33_main))
         D1K3_f22_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta7_f)^((1-cond_reg)*actsf22_reg))
         D1K3_f23_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta7_f)^((1-cond_reg)*actsf23_reg))
         D1K3_f32_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta7_f)^((1-cond_reg)*actsf32_reg))
         D1K3_f33_reg = 1 - ((1-beta7_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta7_f)^((1-cond_reg)*actsf33_reg))
         D1K3_m22_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta7_m)^((1-cond_reg)*actsm22_reg))
         D1K3_m23_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta7_m)^((1-cond_reg)*actsm23_reg))
         D1K3_m32_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta7_m)^((1-cond_reg)*actsm32_reg))
         D1K3_m33_reg = 1 - ((1-beta7_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta7_m)^((1-cond_reg)*actsm33_reg))
         D1K3_f33_casual = 1 - ((1-beta7_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta7_f)^((1-cond_casual)*actsf33_casual))
         D1K3_m33_casual = 1 - ((1-beta7_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta7_m)^((1-cond_casual)*actsm33_casual))
         
         ### Susceptibles on PrEP/not circumcised ###
         
         ## Partner in LV stage, not on HCT/ART
         
         D2G1_f11_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf11_main) * (1-beta8_f)^((1-cond_main)*actsf11_main))
         D2G1_f12_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf12_main) * (1-beta8_f)^((1-cond_main)*actsf12_main))
         D2G1_f13_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf13_main) * (1-beta8_f)^((1-cond_main)*actsf13_main))
         D2G1_f21_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf21_main) * (1-beta8_f)^((1-cond_main)*actsf21_main))
         D2G1_f22_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf22_main) * (1-beta8_f)^((1-cond_main)*actsf22_main))
         D2G1_f23_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf23_main) * (1-beta8_f)^((1-cond_main)*actsf23_main))
         D2G1_f31_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf31_main) * (1-beta8_f)^((1-cond_main)*actsf31_main))
         D2G1_f32_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf32_main) * (1-beta8_f)^((1-cond_main)*actsf32_main))
         D2G1_f33_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*actsf33_main) * (1-beta8_f)^((1-cond_main)*actsf33_main))
         D2G1_m11_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm11_main) * (1-beta8_m)^((1-cond_main)*actsm11_main))
         D2G1_m12_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm12_main) * (1-beta8_m)^((1-cond_main)*actsm12_main))
         D2G1_m13_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm13_main) * (1-beta8_m)^((1-cond_main)*actsm13_main))
         D2G1_m21_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm21_main) * (1-beta8_m)^((1-cond_main)*actsm21_main))
         D2G1_m22_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm22_main) * (1-beta8_m)^((1-cond_main)*actsm22_main))
         D2G1_m23_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm23_main) * (1-beta8_m)^((1-cond_main)*actsm23_main))
         D2G1_m31_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm31_main) * (1-beta8_m)^((1-cond_main)*actsm31_main))
         D2G1_m32_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm32_main) * (1-beta8_m)^((1-cond_main)*actsm32_main))
         D2G1_m33_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*actsm33_main) * (1-beta8_m)^((1-cond_main)*actsm33_main))
         D2G1_f22_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta8_f)^((1-cond_reg)*actsf22_reg))
         D2G1_f23_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta8_f)^((1-cond_reg)*actsf23_reg))
         D2G1_f32_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta8_f)^((1-cond_reg)*actsf32_reg))
         D2G1_f33_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta8_f)^((1-cond_reg)*actsf33_reg))
         D2G1_m22_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta8_m)^((1-cond_reg)*actsm22_reg))
         D2G1_m23_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta8_m)^((1-cond_reg)*actsm23_reg))
         D2G1_m32_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta8_m)^((1-cond_reg)*actsm32_reg))
         D2G1_m33_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta8_m)^((1-cond_reg)*actsm33_reg))
         D2G1_f33_casual = 1 - ((1-beta8_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta8_f)^((1-cond_casual)*actsf33_casual))
         D2G1_m33_casual = 1 - ((1-beta8_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta8_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in LV stage, on HCT
         
         D2G2_f11_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf11_main))
         D2G2_f12_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf12_main))
         D2G2_f13_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf13_main))
         D2G2_f21_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf21_main))
         D2G2_f22_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf22_main))
         D2G2_f23_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf23_main))
         D2G2_f31_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf31_main))
         D2G2_f32_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf32_main))
         D2G2_f33_main = 1 - ((1-beta8_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta8_f)^((1-cond_main)*HCT_cond*HCT_red*actsf33_main))
         D2G2_m11_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
         D2G2_m12_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
         D2G2_m13_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
         D2G2_m21_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
         D2G2_m22_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
         D2G2_m23_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
         D2G2_m31_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
         D2G2_m32_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
         D2G2_m33_main = 1 - ((1-beta8_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta8_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
         D2G2_f22_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf22_reg))
         D2G2_f23_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf23_reg))
         D2G2_f32_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf32_reg))
         D2G2_f33_reg = 1 - ((1-beta8_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta8_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf33_reg))
         D2G2_m22_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
         D2G2_m23_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
         D2G2_m32_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
         D2G2_m33_reg = 1 - ((1-beta8_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta8_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
         D2G2_f33_casual = 1 - ((1-beta8_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta8_f)^((1-cond_casual)*HCT_cond*HCT_red*actsf33_casual))
         D2G2_m33_casual = 1 - ((1-beta8_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta8_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))
         
         ## Partner in LV stage, on ART
         
         D2G3_f11_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf11_main) * (1-beta16_f)^((1-cond_main)*actsf11_main))
         D2G3_f12_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf12_main) * (1-beta16_f)^((1-cond_main)*actsf12_main))
         D2G3_f13_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf13_main) * (1-beta16_f)^((1-cond_main)*actsf13_main))
         D2G3_f21_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf21_main) * (1-beta16_f)^((1-cond_main)*actsf21_main))
         D2G3_f22_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf22_main) * (1-beta16_f)^((1-cond_main)*actsf22_main))
         D2G3_f23_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf23_main) * (1-beta16_f)^((1-cond_main)*actsf23_main))
         D2G3_f31_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf31_main) * (1-beta16_f)^((1-cond_main)*actsf31_main))
         D2G3_f32_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf32_main) * (1-beta16_f)^((1-cond_main)*actsf32_main))
         D2G3_f33_main = 1 - ((1-beta16_f*cond_eff)^(cond_main*actsf33_main) * (1-beta16_f)^((1-cond_main)*actsf33_main))
         D2G3_m11_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm11_main) * (1-beta16_m)^((1-cond_main)*actsm11_main))
         D2G3_m12_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm12_main) * (1-beta16_m)^((1-cond_main)*actsm12_main))
         D2G3_m13_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm13_main) * (1-beta16_m)^((1-cond_main)*actsm13_main))
         D2G3_m21_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm21_main) * (1-beta16_m)^((1-cond_main)*actsm21_main))
         D2G3_m22_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm22_main) * (1-beta16_m)^((1-cond_main)*actsm22_main))
         D2G3_m23_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm23_main) * (1-beta16_m)^((1-cond_main)*actsm23_main))
         D2G3_m31_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm31_main) * (1-beta16_m)^((1-cond_main)*actsm31_main))
         D2G3_m32_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm32_main) * (1-beta16_m)^((1-cond_main)*actsm32_main))
         D2G3_m33_main = 1 - ((1-beta16_m*cond_eff)^(cond_main*actsm33_main) * (1-beta16_m)^((1-cond_main)*actsm33_main))
         D2G3_f22_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta16_f)^((1-cond_reg)*actsf22_reg))
         D2G3_f23_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta16_f)^((1-cond_reg)*actsf23_reg))
         D2G3_f32_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta16_f)^((1-cond_reg)*actsf32_reg))
         D2G3_f33_reg = 1 - ((1-beta16_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta16_f)^((1-cond_reg)*actsf33_reg))
         D2G3_m22_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta16_m)^((1-cond_reg)*actsm22_reg))
         D2G3_m23_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta16_m)^((1-cond_reg)*actsm23_reg))
         D2G3_m32_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta16_m)^((1-cond_reg)*actsm32_reg))
         D2G3_m33_reg = 1 - ((1-beta16_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta16_m)^((1-cond_reg)*actsm33_reg))
         D2G3_f33_casual = 1 - ((1-beta16_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta16_f)^((1-cond_casual)*actsf33_casual))
         D2G3_m33_casual = 1 - ((1-beta16_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta16_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in initial HV stage, not on HCT/ART
         
         D2P1_f11_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf11_main) * (1-beta9_f)^((1-cond_main)*actsf11_main))
         D2P1_f12_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf12_main) * (1-beta9_f)^((1-cond_main)*actsf12_main))
         D2P1_f13_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf13_main) * (1-beta9_f)^((1-cond_main)*actsf13_main))
         D2P1_f21_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf21_main) * (1-beta9_f)^((1-cond_main)*actsf21_main))
         D2P1_f22_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf22_main) * (1-beta9_f)^((1-cond_main)*actsf22_main))
         D2P1_f23_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf23_main) * (1-beta9_f)^((1-cond_main)*actsf23_main))
         D2P1_f31_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf31_main) * (1-beta9_f)^((1-cond_main)*actsf31_main))
         D2P1_f32_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf32_main) * (1-beta9_f)^((1-cond_main)*actsf32_main))
         D2P1_f33_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*actsf33_main) * (1-beta9_f)^((1-cond_main)*actsf33_main))
         D2P1_m11_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm11_main) * (1-beta9_m)^((1-cond_main)*actsm11_main))
         D2P1_m12_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm12_main) * (1-beta9_m)^((1-cond_main)*actsm12_main))
         D2P1_m13_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm13_main) * (1-beta9_m)^((1-cond_main)*actsm13_main))
         D2P1_m21_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm21_main) * (1-beta9_m)^((1-cond_main)*actsm21_main))
         D2P1_m22_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm22_main) * (1-beta9_m)^((1-cond_main)*actsm22_main))
         D2P1_m23_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm23_main) * (1-beta9_m)^((1-cond_main)*actsm23_main))
         D2P1_m31_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm31_main) * (1-beta9_m)^((1-cond_main)*actsm31_main))
         D2P1_m32_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm32_main) * (1-beta9_m)^((1-cond_main)*actsm32_main))
         D2P1_m33_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*actsm33_main) * (1-beta9_m)^((1-cond_main)*actsm33_main))
         D2P1_f22_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta9_f)^((1-cond_reg)*actsf22_reg))
         D2P1_f23_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta9_f)^((1-cond_reg)*actsf23_reg))
         D2P1_f32_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta9_f)^((1-cond_reg)*actsf32_reg))
         D2P1_f33_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta9_f)^((1-cond_reg)*actsf33_reg))
         D2P1_m22_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta9_m)^((1-cond_reg)*actsm22_reg))
         D2P1_m23_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta9_m)^((1-cond_reg)*actsm23_reg))
         D2P1_m32_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta9_m)^((1-cond_reg)*actsm32_reg))
         D2P1_m33_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta9_m)^((1-cond_reg)*actsm33_reg))
         D2P1_f33_casual = 1 - ((1-beta9_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta9_f)^((1-cond_casual)*actsf33_casual))
         D2P1_m33_casual = 1 - ((1-beta9_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta9_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in initial HV stage, on HCT
         
         D2P2_f11_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf11_main))
         D2P2_f12_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf12_main))
         D2P2_f13_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf13_main))
         D2P2_f21_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf21_main))
         D2P2_f22_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf22_main))
         D2P2_f23_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf23_main))
         D2P2_f31_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf31_main))
         D2P2_f32_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf32_main))
         D2P2_f33_main = 1 - ((1-beta9_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta9_f)^((1-cond_main)*HCT_cond*HCT_red*actsf33_main))
         D2P2_m11_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
         D2P2_m12_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
         D2P2_m13_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
         D2P2_m21_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
         D2P2_m22_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
         D2P2_m23_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
         D2P2_m31_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
         D2P2_m32_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
         D2P2_m33_main = 1 - ((1-beta9_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta9_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
         D2P2_f22_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf22_reg))
         D2P2_f23_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf23_reg))
         D2P2_f32_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf32_reg))
         D2P2_f33_reg = 1 - ((1-beta9_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta9_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf33_reg))
         D2P2_m22_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
         D2P2_m23_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
         D2P2_m32_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
         D2P2_m33_reg = 1 - ((1-beta9_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta9_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
         D2P2_f33_casual = 1 - ((1-beta9_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta9_f)^((1-cond_casual)*HCT_cond*HCT_red*actsf33_casual))
         D2P2_m33_casual = 1 - ((1-beta9_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta9_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))
         
         ## Partner in initial HV stage, on ART
         
         D2P3_f11_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf11_main) * (1-beta17_f)^((1-cond_main)*actsf11_main))
         D2P3_f12_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf12_main) * (1-beta17_f)^((1-cond_main)*actsf12_main))
         D2P3_f13_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf13_main) * (1-beta17_f)^((1-cond_main)*actsf13_main))
         D2P3_f21_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf21_main) * (1-beta17_f)^((1-cond_main)*actsf21_main))
         D2P3_f22_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf22_main) * (1-beta17_f)^((1-cond_main)*actsf22_main))
         D2P3_f23_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf23_main) * (1-beta17_f)^((1-cond_main)*actsf23_main))
         D2P3_f31_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf31_main) * (1-beta17_f)^((1-cond_main)*actsf31_main))
         D2P3_f32_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf32_main) * (1-beta17_f)^((1-cond_main)*actsf32_main))
         D2P3_f33_main = 1 - ((1-beta17_f*cond_eff)^(cond_main*actsf33_main) * (1-beta17_f)^((1-cond_main)*actsf33_main))
         D2P3_m11_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm11_main) * (1-beta17_m)^((1-cond_main)*actsm11_main))
         D2P3_m12_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm12_main) * (1-beta17_m)^((1-cond_main)*actsm12_main))
         D2P3_m13_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm13_main) * (1-beta17_m)^((1-cond_main)*actsm13_main))
         D2P3_m21_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm21_main) * (1-beta17_m)^((1-cond_main)*actsm21_main))
         D2P3_m22_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm22_main) * (1-beta17_m)^((1-cond_main)*actsm22_main))
         D2P3_m23_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm23_main) * (1-beta17_m)^((1-cond_main)*actsm23_main))
         D2P3_m31_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm31_main) * (1-beta17_m)^((1-cond_main)*actsm31_main))
         D2P3_m32_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm32_main) * (1-beta17_m)^((1-cond_main)*actsm32_main))
         D2P3_m33_main = 1 - ((1-beta17_m*cond_eff)^(cond_main*actsm33_main) * (1-beta17_m)^((1-cond_main)*actsm33_main))
         D2P3_f22_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta17_f)^((1-cond_reg)*actsf22_reg))
         D2P3_f23_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta17_f)^((1-cond_reg)*actsf23_reg))
         D2P3_f32_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta17_f)^((1-cond_reg)*actsf32_reg))
         D2P3_f33_reg = 1 - ((1-beta17_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta17_f)^((1-cond_reg)*actsf33_reg))
         D2P3_m22_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta17_m)^((1-cond_reg)*actsm22_reg))
         D2P3_m23_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta17_m)^((1-cond_reg)*actsm23_reg))
         D2P3_m32_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta17_m)^((1-cond_reg)*actsm32_reg))
         D2P3_m33_reg = 1 - ((1-beta17_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta17_m)^((1-cond_reg)*actsm33_reg))
         D2P3_f33_casual = 1 - ((1-beta17_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta17_f)^((1-cond_casual)*actsf33_casual))
         D2P3_m33_casual = 1 - ((1-beta17_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta17_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in preAIDS HV stage, not on HCT/ART
         
         D2X1_f11_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf11_main) * (1-beta10_f)^((1-cond_main)*actsf11_main))
         D2X1_f12_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf12_main) * (1-beta10_f)^((1-cond_main)*actsf12_main))
         D2X1_f13_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf13_main) * (1-beta10_f)^((1-cond_main)*actsf13_main))
         D2X1_f21_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf21_main) * (1-beta10_f)^((1-cond_main)*actsf21_main))
         D2X1_f22_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf22_main) * (1-beta10_f)^((1-cond_main)*actsf22_main))
         D2X1_f23_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf23_main) * (1-beta10_f)^((1-cond_main)*actsf23_main))
         D2X1_f31_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf31_main) * (1-beta10_f)^((1-cond_main)*actsf31_main))
         D2X1_f32_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf32_main) * (1-beta10_f)^((1-cond_main)*actsf32_main))
         D2X1_f33_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*actsf33_main) * (1-beta10_f)^((1-cond_main)*actsf33_main))
         D2X1_m11_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm11_main) * (1-beta10_m)^((1-cond_main)*actsm11_main))
         D2X1_m12_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm12_main) * (1-beta10_m)^((1-cond_main)*actsm12_main))
         D2X1_m13_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm13_main) * (1-beta10_m)^((1-cond_main)*actsm13_main))
         D2X1_m21_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm21_main) * (1-beta10_m)^((1-cond_main)*actsm21_main))
         D2X1_m22_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm22_main) * (1-beta10_m)^((1-cond_main)*actsm22_main))
         D2X1_m23_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm23_main) * (1-beta10_m)^((1-cond_main)*actsm23_main))
         D2X1_m31_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm31_main) * (1-beta10_m)^((1-cond_main)*actsm31_main))
         D2X1_m32_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm32_main) * (1-beta10_m)^((1-cond_main)*actsm32_main))
         D2X1_m33_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*actsm33_main) * (1-beta10_m)^((1-cond_main)*actsm33_main))
         D2X1_f22_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta10_f)^((1-cond_reg)*actsf22_reg))
         D2X1_f23_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta10_f)^((1-cond_reg)*actsf23_reg))
         D2X1_f32_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta10_f)^((1-cond_reg)*actsf32_reg))
         D2X1_f33_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta10_f)^((1-cond_reg)*actsf33_reg))
         D2X1_m22_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta10_m)^((1-cond_reg)*actsm22_reg))
         D2X1_m23_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta10_m)^((1-cond_reg)*actsm23_reg))
         D2X1_m32_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta10_m)^((1-cond_reg)*actsm32_reg))
         D2X1_m33_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta10_m)^((1-cond_reg)*actsm33_reg))
         D2X1_f33_casual = 1 - ((1-beta10_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta10_f)^((1-cond_casual)*actsf33_casual))
         D2X1_m33_casual = 1 - ((1-beta10_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta10_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on HCT
         
         D2X2_f11_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf11_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf11_main))
         D2X2_f12_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf12_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf12_main))
         D2X2_f13_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf13_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf13_main))
         D2X2_f21_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf21_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf21_main))
         D2X2_f22_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf22_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf22_main))
         D2X2_f23_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf23_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf23_main))
         D2X2_f31_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf31_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf31_main))
         D2X2_f32_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf32_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf32_main))
         D2X2_f33_main = 1 - ((1-beta10_f*cond_eff)^(cond_main*HCT_cond*HCT_red*actsf33_main) * (1-beta10_f)^((1-cond_main)*HCT_cond*HCT_red*actsf33_main))
         D2X2_m11_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
         D2X2_m12_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
         D2X2_m13_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
         D2X2_m21_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
         D2X2_m22_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
         D2X2_m23_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
         D2X2_m31_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
         D2X2_m32_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
         D2X2_m33_main = 1 - ((1-beta10_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta10_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
         D2X2_f22_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf22_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf22_reg))
         D2X2_f23_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf23_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf23_reg))
         D2X2_f32_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf32_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf32_reg))
         D2X2_f33_reg = 1 - ((1-beta10_f*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsf33_reg) * (1-beta10_f)^((1-cond_reg)*HCT_cond*HCT_red*actsf33_reg))
         D2X2_m22_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
         D2X2_m23_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
         D2X2_m32_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
         D2X2_m33_reg = 1 - ((1-beta10_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta10_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
         D2X2_f33_casual = 1 - ((1-beta10_f*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsf33_casual) * (1-beta10_f)^((1-cond_casual)*HCT_cond*HCT_red*actsf33_casual))
         D2X2_m33_casual = 1 - ((1-beta10_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta10_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on ART
         
         D2X3_f11_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf11_main) * (1-beta18_f)^((1-cond_main)*actsf11_main))
         D2X3_f12_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf12_main) * (1-beta18_f)^((1-cond_main)*actsf12_main))
         D2X3_f13_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf13_main) * (1-beta18_f)^((1-cond_main)*actsf13_main))
         D2X3_f21_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf21_main) * (1-beta18_f)^((1-cond_main)*actsf21_main))
         D2X3_f22_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf22_main) * (1-beta18_f)^((1-cond_main)*actsf22_main))
         D2X3_f23_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf23_main) * (1-beta18_f)^((1-cond_main)*actsf23_main))
         D2X3_f31_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf31_main) * (1-beta18_f)^((1-cond_main)*actsf31_main))
         D2X3_f32_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf32_main) * (1-beta18_f)^((1-cond_main)*actsf32_main))
         D2X3_f33_main = 1 - ((1-beta18_f*cond_eff)^(cond_main*actsf33_main) * (1-beta18_f)^((1-cond_main)*actsf33_main))
         D2X3_m11_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm11_main) * (1-beta18_m)^((1-cond_main)*actsm11_main))
         D2X3_m12_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm12_main) * (1-beta18_m)^((1-cond_main)*actsm12_main))
         D2X3_m13_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm13_main) * (1-beta18_m)^((1-cond_main)*actsm13_main))
         D2X3_m21_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm21_main) * (1-beta18_m)^((1-cond_main)*actsm21_main))
         D2X3_m22_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm22_main) * (1-beta18_m)^((1-cond_main)*actsm22_main))
         D2X3_m23_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm23_main) * (1-beta18_m)^((1-cond_main)*actsm23_main))
         D2X3_m31_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm31_main) * (1-beta18_m)^((1-cond_main)*actsm31_main))
         D2X3_m32_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm32_main) * (1-beta18_m)^((1-cond_main)*actsm32_main))
         D2X3_m33_main = 1 - ((1-beta18_m*cond_eff)^(cond_main*actsm33_main) * (1-beta18_m)^((1-cond_main)*actsm33_main))
         D2X3_f22_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta18_f)^((1-cond_reg)*actsf22_reg))
         D2X3_f23_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta18_f)^((1-cond_reg)*actsf23_reg))
         D2X3_f32_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta18_f)^((1-cond_reg)*actsf32_reg))
         D2X3_f33_reg = 1 - ((1-beta18_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta18_f)^((1-cond_reg)*actsf33_reg))
         D2X3_m22_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta18_m)^((1-cond_reg)*actsm22_reg))
         D2X3_m23_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta18_m)^((1-cond_reg)*actsm23_reg))
         D2X3_m32_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta18_m)^((1-cond_reg)*actsm32_reg))
         D2X3_m33_reg = 1 - ((1-beta18_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta18_m)^((1-cond_reg)*actsm33_reg))
         D2X3_f33_casual = 1 - ((1-beta18_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta18_f)^((1-cond_casual)*actsf33_casual))
         D2X3_m33_casual = 1 - ((1-beta18_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta18_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in AIDS stage, not on HCT/ART
         
         D2K1_f11_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf11_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf11_main))
         D2K1_f12_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf12_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf12_main))
         D2K1_f13_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf13_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf13_main))
         D2K1_f21_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf21_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf21_main))
         D2K1_f22_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf22_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf22_main))
         D2K1_f23_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf23_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf23_main))
         D2K1_f31_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf31_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf31_main))
         D2K1_f32_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf32_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf32_main))
         D2K1_f33_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*AIDS_red*actsf33_main) * (1-beta11_f)^((1-cond_main)*AIDS_red*actsf33_main))
         D2K1_m11_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm11_main))
         D2K1_m12_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm12_main))
         D2K1_m13_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm13_main))
         D2K1_m21_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm21_main))
         D2K1_m22_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm22_main))
         D2K1_m23_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm23_main))
         D2K1_m31_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm31_main))
         D2K1_m32_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm32_main))
         D2K1_m33_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta11_m)^((1-cond_main)*AIDS_red*actsm33_main))
         D2K1_f22_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf22_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf22_reg))
         D2K1_f23_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf23_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf23_reg))
         D2K1_f32_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf32_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf32_reg))
         D2K1_f33_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*AIDS_red*actsf33_reg) * (1-beta11_f)^((1-cond_reg)*AIDS_red*actsf33_reg))
         D2K1_m22_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
         D2K1_m23_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
         D2K1_m32_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
         D2K1_m33_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta11_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
         D2K1_f33_casual = 1 - ((1-beta11_f*cond_eff)^(cond_casual*AIDS_red*actsf33_casual) * (1-beta11_f)^((1-cond_casual)*AIDS_red*actsf33_casual))
         D2K1_m33_casual = 1 - ((1-beta11_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta11_m)^((1-cond_casual)*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on HCT
         
         D2K2_f11_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf11_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf11_main))
         D2K2_f12_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf12_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf12_main))
         D2K2_f13_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf13_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf13_main))
         D2K2_f21_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf21_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf21_main))
         D2K2_f22_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf22_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf22_main))
         D2K2_f23_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf23_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf23_main))
         D2K2_f31_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf31_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf31_main))
         D2K2_f32_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf32_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf32_main))
         D2K2_f33_main = 1 - ((1-beta11_f*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsf33_main) * (1-beta11_f)^((1-cond_main)*HCT_cond*AIDS_red*actsf33_main))
         D2K2_m11_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm11_main))
         D2K2_m12_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm12_main))
         D2K2_m13_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm13_main))
         D2K2_m21_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm21_main))
         D2K2_m22_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm22_main))
         D2K2_m23_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm23_main))
         D2K2_m31_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm31_main))
         D2K2_m32_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm32_main))
         D2K2_m33_main = 1 - ((1-beta11_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta11_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm33_main))
         D2K2_f22_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf22_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf22_reg))
         D2K2_f23_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf23_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf23_reg))
         D2K2_f32_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf32_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf32_reg))
         D2K2_f33_reg = 1 - ((1-beta11_f*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsf33_reg) * (1-beta11_f)^((1-cond_reg)*HCT_cond*AIDS_red*actsf33_reg))
         D2K2_m22_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm22_reg))
         D2K2_m23_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm23_reg))
         D2K2_m32_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm32_reg))
         D2K2_m33_reg = 1 - ((1-beta11_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta11_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm33_reg))
         D2K2_f33_casual = 1 - ((1-beta11_f*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsf33_casual) * (1-beta11_f)^((1-cond_casual)*HCT_cond*AIDS_red*actsf33_casual))
         D2K2_m33_casual = 1 - ((1-beta11_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta11_m)^((1-cond_casual)*HCT_cond*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on ART
         
         D2K3_f11_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf11_main) * (1-beta19_f)^((1-cond_main)*actsf11_main))
         D2K3_f12_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf12_main) * (1-beta19_f)^((1-cond_main)*actsf12_main))
         D2K3_f13_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf13_main) * (1-beta19_f)^((1-cond_main)*actsf13_main))
         D2K3_f21_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf21_main) * (1-beta19_f)^((1-cond_main)*actsf21_main))
         D2K3_f22_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf22_main) * (1-beta19_f)^((1-cond_main)*actsf22_main))
         D2K3_f23_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf23_main) * (1-beta19_f)^((1-cond_main)*actsf23_main))
         D2K3_f31_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf31_main) * (1-beta19_f)^((1-cond_main)*actsf31_main))
         D2K3_f32_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf32_main) * (1-beta19_f)^((1-cond_main)*actsf32_main))
         D2K3_f33_main = 1 - ((1-beta19_f*cond_eff)^(cond_main*actsf33_main) * (1-beta19_f)^((1-cond_main)*actsf33_main))
         D2K3_m11_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm11_main) * (1-beta19_m)^((1-cond_main)*actsm11_main))
         D2K3_m12_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm12_main) * (1-beta19_m)^((1-cond_main)*actsm12_main))
         D2K3_m13_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm13_main) * (1-beta19_m)^((1-cond_main)*actsm13_main))
         D2K3_m21_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm21_main) * (1-beta19_m)^((1-cond_main)*actsm21_main))
         D2K3_m22_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm22_main) * (1-beta19_m)^((1-cond_main)*actsm22_main))
         D2K3_m23_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm23_main) * (1-beta19_m)^((1-cond_main)*actsm23_main))
         D2K3_m31_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm31_main) * (1-beta19_m)^((1-cond_main)*actsm31_main))
         D2K3_m32_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm32_main) * (1-beta19_m)^((1-cond_main)*actsm32_main))
         D2K3_m33_main = 1 - ((1-beta19_m*cond_eff)^(cond_main*actsm33_main) * (1-beta19_m)^((1-cond_main)*actsm33_main))
         D2K3_f22_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf22_reg) * (1-beta19_f)^((1-cond_reg)*actsf22_reg))
         D2K3_f23_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf23_reg) * (1-beta19_f)^((1-cond_reg)*actsf23_reg))
         D2K3_f32_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf32_reg) * (1-beta19_f)^((1-cond_reg)*actsf32_reg))
         D2K3_f33_reg = 1 - ((1-beta19_f*cond_eff)^(cond_reg*actsf33_reg) * (1-beta19_f)^((1-cond_reg)*actsf33_reg))
         D2K3_m22_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta19_m)^((1-cond_reg)*actsm22_reg))
         D2K3_m23_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta19_m)^((1-cond_reg)*actsm23_reg))
         D2K3_m32_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta19_m)^((1-cond_reg)*actsm32_reg))
         D2K3_m33_reg = 1 - ((1-beta19_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta19_m)^((1-cond_reg)*actsm33_reg))
         D2K3_f33_casual = 1 - ((1-beta19_f*cond_eff)^(cond_casual*actsf33_casual) * (1-beta19_f)^((1-cond_casual)*actsf33_casual))
         D2K3_m33_casual = 1 - ((1-beta19_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta19_m)^((1-cond_casual)*actsm33_casual))
         
         ### Circumcised males not on PrEP ###
         
         ## Partner in LV stage, not on HCT/ART
         
         D3G1_m11_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm11_main) * (1-beta12_m)^((1-cond_main)*actsm11_main))
         D3G1_m12_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm12_main) * (1-beta12_m)^((1-cond_main)*actsm12_main))
         D3G1_m13_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm13_main) * (1-beta12_m)^((1-cond_main)*actsm13_main))
         D3G1_m21_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm21_main) * (1-beta12_m)^((1-cond_main)*actsm21_main))
         D3G1_m22_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm22_main) * (1-beta12_m)^((1-cond_main)*actsm22_main))
         D3G1_m23_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm23_main) * (1-beta12_m)^((1-cond_main)*actsm23_main))
         D3G1_m31_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm31_main) * (1-beta12_m)^((1-cond_main)*actsm31_main))
         D3G1_m32_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm32_main) * (1-beta12_m)^((1-cond_main)*actsm32_main))
         D3G1_m33_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*actsm33_main) * (1-beta12_m)^((1-cond_main)*actsm33_main))
         D3G1_m22_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta12_m)^((1-cond_reg)*actsm22_reg))
         D3G1_m23_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta12_m)^((1-cond_reg)*actsm23_reg))
         D3G1_m32_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta12_m)^((1-cond_reg)*actsm32_reg))
         D3G1_m33_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta12_m)^((1-cond_reg)*actsm33_reg))
         D3G1_m33_casual = 1 - ((1-beta12_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta12_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in LV stage, on HCT
         
         D3G2_m11_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
         D3G2_m12_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
         D3G2_m13_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
         D3G2_m21_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
         D3G2_m22_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
         D3G2_m23_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
         D3G2_m31_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
         D3G2_m32_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
         D3G2_m33_main = 1 - ((1-beta12_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta12_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
         D3G2_m22_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
         D3G2_m23_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
         D3G2_m32_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
         D3G2_m33_reg = 1 - ((1-beta12_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta12_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
         D3G2_m33_casual = 1 - ((1-beta12_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta12_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))
         
         ## Partner in LV stage, on ART
         
         D3G3_m11_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm11_main) * (1-beta20_m)^((1-cond_main)*actsm11_main))
         D3G3_m12_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm12_main) * (1-beta20_m)^((1-cond_main)*actsm12_main))
         D3G3_m13_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm13_main) * (1-beta20_m)^((1-cond_main)*actsm13_main))
         D3G3_m21_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm21_main) * (1-beta20_m)^((1-cond_main)*actsm21_main))
         D3G3_m22_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm22_main) * (1-beta20_m)^((1-cond_main)*actsm22_main))
         D3G3_m23_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm23_main) * (1-beta20_m)^((1-cond_main)*actsm23_main))
         D3G3_m31_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm31_main) * (1-beta20_m)^((1-cond_main)*actsm31_main))
         D3G3_m32_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm32_main) * (1-beta20_m)^((1-cond_main)*actsm32_main))
         D3G3_m33_main = 1 - ((1-beta20_m*cond_eff)^(cond_main*actsm33_main) * (1-beta20_m)^((1-cond_main)*actsm33_main))
         D3G3_m22_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta20_m)^((1-cond_reg)*actsm22_reg))
         D3G3_m23_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta20_m)^((1-cond_reg)*actsm23_reg))
         D3G3_m32_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta20_m)^((1-cond_reg)*actsm32_reg))
         D3G3_m33_reg = 1 - ((1-beta20_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta20_m)^((1-cond_reg)*actsm33_reg))
         D3G3_m33_casual = 1 - ((1-beta20_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta20_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in initial HV stage, not on HCT/ART
         
         D3P1_m11_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm11_main) * (1-beta13_m)^((1-cond_main)*actsm11_main))
         D3P1_m12_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm12_main) * (1-beta13_m)^((1-cond_main)*actsm12_main))
         D3P1_m13_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm13_main) * (1-beta13_m)^((1-cond_main)*actsm13_main))
         D3P1_m21_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm21_main) * (1-beta13_m)^((1-cond_main)*actsm21_main))
         D3P1_m22_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm22_main) * (1-beta13_m)^((1-cond_main)*actsm22_main))
         D3P1_m23_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm23_main) * (1-beta13_m)^((1-cond_main)*actsm23_main))
         D3P1_m31_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm31_main) * (1-beta13_m)^((1-cond_main)*actsm31_main))
         D3P1_m32_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm32_main) * (1-beta13_m)^((1-cond_main)*actsm32_main))
         D3P1_m33_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*actsm33_main) * (1-beta13_m)^((1-cond_main)*actsm33_main))
         D3P1_m22_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta13_m)^((1-cond_reg)*actsm22_reg))
         D3P1_m23_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta13_m)^((1-cond_reg)*actsm23_reg))
         D3P1_m32_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta13_m)^((1-cond_reg)*actsm32_reg))
         D3P1_m33_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta13_m)^((1-cond_reg)*actsm33_reg))
         D3P1_m33_casual = 1 - ((1-beta13_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta13_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in initial HV stage, on HCT
         
         D3P2_m11_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
         D3P2_m12_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
         D3P2_m13_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
         D3P2_m21_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
         D3P2_m22_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
         D3P2_m23_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
         D3P2_m31_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
         D3P2_m32_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
         D3P2_m33_main = 1 - ((1-beta13_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta13_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
         D3P2_m22_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
         D3P2_m23_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
         D3P2_m32_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
         D3P2_m33_reg = 1 - ((1-beta13_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta13_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
         D3P2_m33_casual = 1 - ((1-beta13_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta13_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))
         
         ## Partner in initial HV stage, on ART
         
         D3P3_m11_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm11_main) * (1-beta21_m)^((1-cond_main)*actsm11_main))
         D3P3_m12_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm12_main) * (1-beta21_m)^((1-cond_main)*actsm12_main))
         D3P3_m13_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm13_main) * (1-beta21_m)^((1-cond_main)*actsm13_main))
         D3P3_m21_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm21_main) * (1-beta21_m)^((1-cond_main)*actsm21_main))
         D3P3_m22_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm22_main) * (1-beta21_m)^((1-cond_main)*actsm22_main))
         D3P3_m23_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm23_main) * (1-beta21_m)^((1-cond_main)*actsm23_main))
         D3P3_m31_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm31_main) * (1-beta21_m)^((1-cond_main)*actsm31_main))
         D3P3_m32_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm32_main) * (1-beta21_m)^((1-cond_main)*actsm32_main))
         D3P3_m33_main = 1 - ((1-beta21_m*cond_eff)^(cond_main*actsm33_main) * (1-beta21_m)^((1-cond_main)*actsm33_main))
         D3P3_m22_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta21_m)^((1-cond_reg)*actsm22_reg))
         D3P3_m23_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta21_m)^((1-cond_reg)*actsm23_reg))
         D3P3_m32_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta21_m)^((1-cond_reg)*actsm32_reg))
         D3P3_m33_reg = 1 - ((1-beta21_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta21_m)^((1-cond_reg)*actsm33_reg))
         D3P3_m33_casual = 1 - ((1-beta21_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta21_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in preAIDS HV stage, not on HCT/ART
         
         D3X1_m11_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm11_main) * (1-beta14_m)^((1-cond_main)*actsm11_main))
         D3X1_m12_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm12_main) * (1-beta14_m)^((1-cond_main)*actsm12_main))
         D3X1_m13_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm13_main) * (1-beta14_m)^((1-cond_main)*actsm13_main))
         D3X1_m21_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm21_main) * (1-beta14_m)^((1-cond_main)*actsm21_main))
         D3X1_m22_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm22_main) * (1-beta14_m)^((1-cond_main)*actsm22_main))
         D3X1_m23_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm23_main) * (1-beta14_m)^((1-cond_main)*actsm23_main))
         D3X1_m31_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm31_main) * (1-beta14_m)^((1-cond_main)*actsm31_main))
         D3X1_m32_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm32_main) * (1-beta14_m)^((1-cond_main)*actsm32_main))
         D3X1_m33_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*actsm33_main) * (1-beta14_m)^((1-cond_main)*actsm33_main))
         D3X1_m22_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta14_m)^((1-cond_reg)*actsm22_reg))
         D3X1_m23_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta14_m)^((1-cond_reg)*actsm23_reg))
         D3X1_m32_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta14_m)^((1-cond_reg)*actsm32_reg))
         D3X1_m33_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta14_m)^((1-cond_reg)*actsm33_reg))
         D3X1_m33_casual = 1 - ((1-beta14_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta14_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on HCT
         
         D3X2_m11_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm11_main))
         D3X2_m12_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm12_main))
         D3X2_m13_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm13_main))
         D3X2_m21_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm21_main))
         D3X2_m22_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm22_main))
         D3X2_m23_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm23_main))
         D3X2_m31_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm31_main))
         D3X2_m32_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm32_main))
         D3X2_m33_main = 1 - ((1-beta14_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta14_m)^((1-cond_main*HCT_cond)*HCT_red*actsm33_main))
         D3X2_m22_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm22_reg))
         D3X2_m23_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm23_reg))
         D3X2_m32_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm32_reg))
         D3X2_m33_reg = 1 - ((1-beta14_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta14_m)^((1-cond_reg*HCT_cond)*HCT_red*actsm33_reg))
         D3X2_m33_casual = 1 - ((1-beta14_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta14_m)^((1-cond_casual*HCT_cond)*HCT_red*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on ART
         
         D3X3_m11_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm11_main) * (1-beta22_m)^((1-cond_main)*actsm11_main))
         D3X3_m12_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm12_main) * (1-beta22_m)^((1-cond_main)*actsm12_main))
         D3X3_m13_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm13_main) * (1-beta22_m)^((1-cond_main)*actsm13_main))
         D3X3_m21_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm21_main) * (1-beta22_m)^((1-cond_main)*actsm21_main))
         D3X3_m22_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm22_main) * (1-beta22_m)^((1-cond_main)*actsm22_main))
         D3X3_m23_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm23_main) * (1-beta22_m)^((1-cond_main)*actsm23_main))
         D3X3_m31_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm31_main) * (1-beta22_m)^((1-cond_main)*actsm31_main))
         D3X3_m32_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm32_main) * (1-beta22_m)^((1-cond_main)*actsm32_main))
         D3X3_m33_main = 1 - ((1-beta22_m*cond_eff)^(cond_main*actsm33_main) * (1-beta22_m)^((1-cond_main)*actsm33_main))
         D3X3_m22_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta22_m)^((1-cond_reg)*actsm22_reg))
         D3X3_m23_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta22_m)^((1-cond_reg)*actsm23_reg))
         D3X3_m32_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta22_m)^((1-cond_reg)*actsm32_reg))
         D3X3_m33_reg = 1 - ((1-beta22_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta22_m)^((1-cond_reg)*actsm33_reg))
         D3X3_m33_casual = 1 - ((1-beta22_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta22_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in AIDS stage, not on HCT/ART
         
         D3K1_m11_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm11_main))
         D3K1_m12_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm12_main))
         D3K1_m13_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm13_main))
         D3K1_m21_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm21_main))
         D3K1_m22_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm22_main))
         D3K1_m23_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm23_main))
         D3K1_m31_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm31_main))
         D3K1_m32_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm32_main))
         D3K1_m33_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta15_m)^((1-cond_main)*AIDS_red*actsm33_main))
         D3K1_m22_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
         D3K1_m23_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
         D3K1_m32_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
         D3K1_m33_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta15_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
         D3K1_m33_casual = 1 - ((1-beta15_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta15_m)^((1-cond_casual)*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on HCT
         
         D3K2_m11_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm11_main))
         D3K2_m12_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm12_main))
         D3K2_m13_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm13_main))
         D3K2_m21_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm21_main))
         D3K2_m22_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm22_main))
         D3K2_m23_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm23_main))
         D3K2_m31_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm31_main))
         D3K2_m32_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm32_main))
         D3K2_m33_main = 1 - ((1-beta15_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta15_m)^((1-cond_main*HCT_cond)*AIDS_red*actsm33_main))
         D3K2_m22_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm22_reg))
         D3K2_m23_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm23_reg))
         D3K2_m32_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm32_reg))
         D3K2_m33_reg = 1 - ((1-beta15_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta15_m)^((1-cond_reg*HCT_cond)*AIDS_red*actsm33_reg))
         D3K2_m33_casual = 1 - ((1-beta15_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta15_m)^((1-cond_casual*HCT_cond)*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on ART
         
         D3K3_m11_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm11_main) * (1-beta23_m)^((1-cond_main)*actsm11_main))
         D3K3_m12_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm12_main) * (1-beta23_m)^((1-cond_main)*actsm12_main))
         D3K3_m13_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm13_main) * (1-beta23_m)^((1-cond_main)*actsm13_main))
         D3K3_m21_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm21_main) * (1-beta23_m)^((1-cond_main)*actsm21_main))
         D3K3_m22_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm22_main) * (1-beta23_m)^((1-cond_main)*actsm22_main))
         D3K3_m23_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm23_main) * (1-beta23_m)^((1-cond_main)*actsm23_main))
         D3K3_m31_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm31_main) * (1-beta23_m)^((1-cond_main)*actsm31_main))
         D3K3_m32_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm32_main) * (1-beta23_m)^((1-cond_main)*actsm32_main))
         D3K3_m33_main = 1 - ((1-beta23_m*cond_eff)^(cond_main*actsm33_main) * (1-beta23_m)^((1-cond_main)*actsm33_main))
         D3K3_m22_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta23_m)^((1-cond_reg)*actsm22_reg))
         D3K3_m23_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta23_m)^((1-cond_reg)*actsm23_reg))
         D3K3_m32_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta23_m)^((1-cond_reg)*actsm32_reg))
         D3K3_m33_reg = 1 - ((1-beta23_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta23_m)^((1-cond_reg)*actsm33_reg))
         D3K3_m33_casual = 1 - ((1-beta23_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta23_m)^((1-cond_casual)*actsm33_casual))
         
         ### Circumcised males on PrEP ###
         
         ## Partner in LV stage, not on HCT/ART
         
         D4G1_m11_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm11_main) * (1-beta28_m)^((1-cond_main)*actsm11_main))
         D4G1_m12_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm12_main) * (1-beta28_m)^((1-cond_main)*actsm12_main))
         D4G1_m13_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm13_main) * (1-beta28_m)^((1-cond_main)*actsm13_main))
         D4G1_m21_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm21_main) * (1-beta28_m)^((1-cond_main)*actsm21_main))
         D4G1_m22_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm22_main) * (1-beta28_m)^((1-cond_main)*actsm22_main))
         D4G1_m23_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm23_main) * (1-beta28_m)^((1-cond_main)*actsm23_main))
         D4G1_m31_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm31_main) * (1-beta28_m)^((1-cond_main)*actsm31_main))
         D4G1_m32_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm32_main) * (1-beta28_m)^((1-cond_main)*actsm32_main))
         D4G1_m33_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*actsm33_main) * (1-beta28_m)^((1-cond_main)*actsm33_main))
         D4G1_m22_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta28_m)^((1-cond_reg)*actsm22_reg))
         D4G1_m23_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta28_m)^((1-cond_reg)*actsm23_reg))
         D4G1_m32_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta28_m)^((1-cond_reg)*actsm32_reg))
         D4G1_m33_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta28_m)^((1-cond_reg)*actsm33_reg))
         D4G1_m33_casual = 1 - ((1-beta28_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta28_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in LV stage, on HCT
         
         D4G2_m11_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
         D4G2_m12_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
         D4G2_m13_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
         D4G2_m21_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
         D4G2_m22_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
         D4G2_m23_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
         D4G2_m31_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
         D4G2_m32_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
         D4G2_m33_main = 1 - ((1-beta28_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta28_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
         D4G2_m22_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
         D4G2_m23_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
         D4G2_m32_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
         D4G2_m33_reg = 1 - ((1-beta28_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta28_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
         D4G2_m33_casual = 1 - ((1-beta28_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta28_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))
         
         ## Partner in LV stage, on ART
         
         D4G3_m11_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm11_main) * (1-beta24_m)^((1-cond_main)*actsm11_main))
         D4G3_m12_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm12_main) * (1-beta24_m)^((1-cond_main)*actsm12_main))
         D4G3_m13_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm13_main) * (1-beta24_m)^((1-cond_main)*actsm13_main))
         D4G3_m21_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm21_main) * (1-beta24_m)^((1-cond_main)*actsm21_main))
         D4G3_m22_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm22_main) * (1-beta24_m)^((1-cond_main)*actsm22_main))
         D4G3_m23_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm23_main) * (1-beta24_m)^((1-cond_main)*actsm23_main))
         D4G3_m31_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm31_main) * (1-beta24_m)^((1-cond_main)*actsm31_main))
         D4G3_m32_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm32_main) * (1-beta24_m)^((1-cond_main)*actsm32_main))
         D4G3_m33_main = 1 - ((1-beta24_m*cond_eff)^(cond_main*actsm33_main) * (1-beta24_m)^((1-cond_main)*actsm33_main))
         D4G3_m22_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta24_m)^((1-cond_reg)*actsm22_reg))
         D4G3_m23_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta24_m)^((1-cond_reg)*actsm23_reg))
         D4G3_m32_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta24_m)^((1-cond_reg)*actsm32_reg))
         D4G3_m33_reg = 1 - ((1-beta24_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta24_m)^((1-cond_reg)*actsm33_reg))
         D4G3_m33_casual = 1 - ((1-beta24_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta24_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in initial HV stage, not on HCT/ART
         
         D4P1_m11_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm11_main) * (1-beta29_m)^((1-cond_main)*actsm11_main))
         D4P1_m12_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm12_main) * (1-beta29_m)^((1-cond_main)*actsm12_main))
         D4P1_m13_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm13_main) * (1-beta29_m)^((1-cond_main)*actsm13_main))
         D4P1_m21_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm21_main) * (1-beta29_m)^((1-cond_main)*actsm21_main))
         D4P1_m22_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm22_main) * (1-beta29_m)^((1-cond_main)*actsm22_main))
         D4P1_m23_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm23_main) * (1-beta29_m)^((1-cond_main)*actsm23_main))
         D4P1_m31_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm31_main) * (1-beta29_m)^((1-cond_main)*actsm31_main))
         D4P1_m32_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm32_main) * (1-beta29_m)^((1-cond_main)*actsm32_main))
         D4P1_m33_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*actsm33_main) * (1-beta29_m)^((1-cond_main)*actsm33_main))
         D4P1_m22_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta29_m)^((1-cond_reg)*actsm22_reg))
         D4P1_m23_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta29_m)^((1-cond_reg)*actsm23_reg))
         D4P1_m32_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta29_m)^((1-cond_reg)*actsm32_reg))
         D4P1_m33_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta29_m)^((1-cond_reg)*actsm33_reg))
         D4P1_m33_casual = 1 - ((1-beta29_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta29_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in initial HV stage, on HCT
         
         D4P2_m11_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
         D4P2_m12_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
         D4P2_m13_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
         D4P2_m21_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
         D4P2_m22_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
         D4P2_m23_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
         D4P2_m31_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
         D4P2_m32_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
         D4P2_m33_main = 1 - ((1-beta29_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta29_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
         D4P2_m22_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
         D4P2_m23_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
         D4P2_m32_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
         D4P2_m33_reg = 1 - ((1-beta29_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta29_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
         D4P2_m33_casual = 1 - ((1-beta29_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta29_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))
         
         ## Partner in initial HV stage, on ART
         
         D4P3_m11_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm11_main) * (1-beta25_m)^((1-cond_main)*actsm11_main))
         D4P3_m12_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm12_main) * (1-beta25_m)^((1-cond_main)*actsm12_main))
         D4P3_m13_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm13_main) * (1-beta25_m)^((1-cond_main)*actsm13_main))
         D4P3_m21_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm21_main) * (1-beta25_m)^((1-cond_main)*actsm21_main))
         D4P3_m22_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm22_main) * (1-beta25_m)^((1-cond_main)*actsm22_main))
         D4P3_m23_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm23_main) * (1-beta25_m)^((1-cond_main)*actsm23_main))
         D4P3_m31_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm31_main) * (1-beta25_m)^((1-cond_main)*actsm31_main))
         D4P3_m32_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm32_main) * (1-beta25_m)^((1-cond_main)*actsm32_main))
         D4P3_m33_main = 1 - ((1-beta25_m*cond_eff)^(cond_main*actsm33_main) * (1-beta25_m)^((1-cond_main)*actsm33_main))
         D4P3_m22_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta25_m)^((1-cond_reg)*actsm22_reg))
         D4P3_m23_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta25_m)^((1-cond_reg)*actsm23_reg))
         D4P3_m32_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta25_m)^((1-cond_reg)*actsm32_reg))
         D4P3_m33_reg = 1 - ((1-beta25_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta25_m)^((1-cond_reg)*actsm33_reg))
         D4P3_m33_casual = 1 - ((1-beta25_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta25_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in preAIDS HV stage, not on HCT/ART
         
         D4X1_m11_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm11_main) * (1-beta30_m)^((1-cond_main)*actsm11_main))
         D4X1_m12_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm12_main) * (1-beta30_m)^((1-cond_main)*actsm12_main))
         D4X1_m13_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm13_main) * (1-beta30_m)^((1-cond_main)*actsm13_main))
         D4X1_m21_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm21_main) * (1-beta30_m)^((1-cond_main)*actsm21_main))
         D4X1_m22_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm22_main) * (1-beta30_m)^((1-cond_main)*actsm22_main))
         D4X1_m23_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm23_main) * (1-beta30_m)^((1-cond_main)*actsm23_main))
         D4X1_m31_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm31_main) * (1-beta30_m)^((1-cond_main)*actsm31_main))
         D4X1_m32_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm32_main) * (1-beta30_m)^((1-cond_main)*actsm32_main))
         D4X1_m33_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*actsm33_main) * (1-beta30_m)^((1-cond_main)*actsm33_main))
         D4X1_m22_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta30_m)^((1-cond_reg)*actsm22_reg))
         D4X1_m23_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta30_m)^((1-cond_reg)*actsm23_reg))
         D4X1_m32_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta30_m)^((1-cond_reg)*actsm32_reg))
         D4X1_m33_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta30_m)^((1-cond_reg)*actsm33_reg))
         D4X1_m33_casual = 1 - ((1-beta30_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta30_m)^((1-cond_casual)*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on HCT
         
         D4X2_m11_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm11_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm11_main))
         D4X2_m12_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm12_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm12_main))
         D4X2_m13_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm13_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm13_main))
         D4X2_m21_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm21_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm21_main))
         D4X2_m22_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm22_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm22_main))
         D4X2_m23_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm23_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm23_main))
         D4X2_m31_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm31_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm31_main))
         D4X2_m32_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm32_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm32_main))
         D4X2_m33_main = 1 - ((1-beta30_m*cond_eff)^(cond_main*HCT_cond*HCT_red*actsm33_main) * (1-beta30_m)^((1-cond_main)*HCT_cond*HCT_red*actsm33_main))
         D4X2_m22_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm22_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm22_reg))
         D4X2_m23_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm23_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm23_reg))
         D4X2_m32_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm32_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm32_reg))
         D4X2_m33_reg = 1 - ((1-beta30_m*cond_eff)^(cond_reg*HCT_cond*HCT_red*actsm33_reg) * (1-beta30_m)^((1-cond_reg)*HCT_cond*HCT_red*actsm33_reg))
         D4X2_m33_casual = 1 - ((1-beta30_m*cond_eff)^(cond_casual*HCT_cond*HCT_red*actsm33_casual) * (1-beta30_m)^((1-cond_casual)*HCT_cond*HCT_red*actsm33_casual))
         
         ## Partner in preAIDS HV stage, on ART
         
         D4X3_m11_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm11_main) * (1-beta26_m)^((1-cond_main)*actsm11_main))
         D4X3_m12_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm12_main) * (1-beta26_m)^((1-cond_main)*actsm12_main))
         D4X3_m13_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm13_main) * (1-beta26_m)^((1-cond_main)*actsm13_main))
         D4X3_m21_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm21_main) * (1-beta26_m)^((1-cond_main)*actsm21_main))
         D4X3_m22_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm22_main) * (1-beta26_m)^((1-cond_main)*actsm22_main))
         D4X3_m23_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm23_main) * (1-beta26_m)^((1-cond_main)*actsm23_main))
         D4X3_m31_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm31_main) * (1-beta26_m)^((1-cond_main)*actsm31_main))
         D4X3_m32_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm32_main) * (1-beta26_m)^((1-cond_main)*actsm32_main))
         D4X3_m33_main = 1 - ((1-beta26_m*cond_eff)^(cond_main*actsm33_main) * (1-beta26_m)^((1-cond_main)*actsm33_main))
         D4X3_m22_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta26_m)^((1-cond_reg)*actsm22_reg))
         D4X3_m23_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta26_m)^((1-cond_reg)*actsm23_reg))
         D4X3_m32_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta26_m)^((1-cond_reg)*actsm32_reg))
         D4X3_m33_reg = 1 - ((1-beta26_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta26_m)^((1-cond_reg)*actsm33_reg))
         D4X3_m33_casual = 1 - ((1-beta26_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta26_m)^((1-cond_casual)*actsm33_casual))
         
         #####
         ## Partner in AIDS stage, not on HCT/ART
         
         D4K1_m11_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm11_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm11_main))
         D4K1_m12_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm12_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm12_main))
         D4K1_m13_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm13_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm13_main))
         D4K1_m21_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm21_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm21_main))
         D4K1_m22_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm22_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm22_main))
         D4K1_m23_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm23_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm23_main))
         D4K1_m31_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm31_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm31_main))
         D4K1_m32_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm32_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm32_main))
         D4K1_m33_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*AIDS_red*actsm33_main) * (1-beta31_m)^((1-cond_main)*AIDS_red*actsm33_main))
         D4K1_m22_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm22_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm22_reg))
         D4K1_m23_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm23_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm23_reg))
         D4K1_m32_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm32_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm32_reg))
         D4K1_m33_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*AIDS_red*actsm33_reg) * (1-beta31_m)^((1-cond_reg)*AIDS_red*actsm33_reg))
         D4K1_m33_casual = 1 - ((1-beta31_m*cond_eff)^(cond_casual*AIDS_red*actsm33_casual) * (1-beta31_m)^((1-cond_casual)*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on HCT
         
         D4K2_m11_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm11_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm11_main))
         D4K2_m12_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm12_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm12_main))
         D4K2_m13_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm13_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm13_main))
         D4K2_m21_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm21_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm21_main))
         D4K2_m22_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm22_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm22_main))
         D4K2_m23_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm23_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm23_main))
         D4K2_m31_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm31_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm31_main))
         D4K2_m32_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm32_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm32_main))
         D4K2_m33_main = 1 - ((1-beta31_m*cond_eff)^(cond_main*HCT_cond*AIDS_red*actsm33_main) * (1-beta31_m)^((1-cond_main)*HCT_cond*AIDS_red*actsm33_main))
         D4K2_m22_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm22_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm22_reg))
         D4K2_m23_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm23_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm23_reg))
         D4K2_m32_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm32_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm32_reg))
         D4K2_m33_reg = 1 - ((1-beta31_m*cond_eff)^(cond_reg*HCT_cond*AIDS_red*actsm33_reg) * (1-beta31_m)^((1-cond_reg)*HCT_cond*AIDS_red*actsm33_reg))
         D4K2_m33_casual = 1 - ((1-beta31_m*cond_eff)^(cond_casual*HCT_cond*AIDS_red*actsm33_casual) * (1-beta31_m)^((1-cond_casual)*HCT_cond*AIDS_red*actsm33_casual))
         
         ## Partner in AIDS stage, on ART
         
         D4K3_m11_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm11_main) * (1-beta27_m)^((1-cond_main)*actsm11_main))
         D4K3_m12_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm12_main) * (1-beta27_m)^((1-cond_main)*actsm12_main))
         D4K3_m13_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm13_main) * (1-beta27_m)^((1-cond_main)*actsm13_main))
         D4K3_m21_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm21_main) * (1-beta27_m)^((1-cond_main)*actsm21_main))
         D4K3_m22_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm22_main) * (1-beta27_m)^((1-cond_main)*actsm22_main))
         D4K3_m23_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm23_main) * (1-beta27_m)^((1-cond_main)*actsm23_main))
         D4K3_m31_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm31_main) * (1-beta27_m)^((1-cond_main)*actsm31_main))
         D4K3_m32_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm32_main) * (1-beta27_m)^((1-cond_main)*actsm32_main))
         D4K3_m33_main = 1 - ((1-beta27_m*cond_eff)^(cond_main*actsm33_main) * (1-beta27_m)^((1-cond_main)*actsm33_main))
         D4K3_m22_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm22_reg) * (1-beta27_m)^((1-cond_reg)*actsm22_reg))
         D4K3_m23_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm23_reg) * (1-beta27_m)^((1-cond_reg)*actsm23_reg))
         D4K3_m32_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm32_reg) * (1-beta27_m)^((1-cond_reg)*actsm32_reg))
         D4K3_m33_reg = 1 - ((1-beta27_m*cond_eff)^(cond_reg*actsm33_reg) * (1-beta27_m)^((1-cond_reg)*actsm33_reg))
         D4K3_m33_casual = 1 - ((1-beta27_m*cond_eff)^(cond_casual*actsm33_casual) * (1-beta27_m)^((1-cond_casual)*actsm33_casual))
         
         ##########################
         ### FORCE OF INFECTION ###
         
         # Females not on PrEP
         pi1_f1 = (
           # low-low
           (df11_main*rho_f1_1_main)*(
             # infected partner not on HCT/ART
             D1P1_f11_main*(H1m1/Nm1) + D1G1_f11_main*(Y1m1/Nm1) + D1X1_f11_main*(Z1m1/Nm1) + D1K1_f11_main*(A1m1/Nm1) # = (D1P1_f11_main*H1m1 + D1G1_f11_main*Y1m1 + D1X1_f11_main*Z1m1 + D1K1_f11_main*A1m1)/Nm1
             # infected partner on HCT
             + D1P2_f11_main*(H2m1/Nm1) + D1G2_f11_main*(Y2m1/Nm1) + D1X2_f11_main*(Z2m1/Nm1) + D1K2_f11_main*(A2m1/Nm1)
             # infected partner on ART
             + D1P3_f11_main*(H3m1/Nm1) + D1G3_f11_main*(Y3m1/Nm1) + D1X3_f11_main*(Z3m1/Nm1) + D1K3_f11_main*(A3m1/Nm1)
           )
           # low-medium
           +(df12_main*rho_f1_2_main)*(
             # infected partner not on HCT/ART
             D1P1_f12_main*(H1m2/Nm2) + D1G1_f12_main*(Y1m2/Nm2) + D1X1_f12_main*(Z1m2/Nm2) + D1K1_f12_main*(A1m2/Nm2)
             # infected partner on HCT
             + D1P2_f12_main*(H2m2/Nm2) + D1G2_f12_main*(Y2m2/Nm2) + D1X2_f12_main*(Z2m2/Nm2) + D1K2_f12_main*(A2m2/Nm2)
             # infected partner on ART
             + D1P3_f12_main*(H3m2/Nm2) + D1G3_f12_main*(Y3m2/Nm2) + D1X3_f12_main*(Z3m2/Nm2) + D1K3_f12_main*(A3m2/Nm2)
           )
           # low-high
           +(df13_main*rho_f1_3_main)*(
             # infected partner not on HCT/ART
             D1P1_f13_main*(H1m3/Nm3) + D1G1_f13_main*(Y1m3/Nm3) + D1X1_f13_main*(Z1m3/Nm3) + D1K1_f13_main*(A1m3/Nm3)
             # infected partner on HCT
             + D1P2_f13_main*(H2m3/Nm3) + D1G2_f13_main*(Y2m3/Nm3) + D1X2_f13_main*(Z2m3/Nm3) + D1K2_f13_main*(A2m3/Nm3)
             # infected partner on ART
             + D1P3_f13_main*(H3m3/Nm3) + D1G3_f13_main*(Y3m3/Nm3) + D1X3_f13_main*(Z3m3/Nm3) + D1K3_f13_main*(A3m3/Nm3)
           )
         )
         
         
         pi1_f2 = (
           # main
           (df21_main*rho_f2_1_main)*(
             # infected partner not on HCT/ART
             D1P1_f21_main*(H1m1/Nm1) + D1G1_f21_main*(Y1m1/Nm1) + D1X1_f21_main*(Z1m1/Nm1) + D1K1_f21_main*(A1m1/Nm1)
             # infected partner on HCT
             + D1P2_f21_main*(H2m1/Nm1) + D1G2_f21_main*(Y2m1/Nm1) + D1X2_f21_main*(Z2m1/Nm1) + D1K2_f21_main*(A2m1/Nm1)
             # infected partner on ART
             + D1P3_f21_main*(H3m1/Nm1) + D1G3_f21_main*(Y3m1/Nm1) + D1X3_f21_main*(Z3m1/Nm1) + D1K3_f21_main*(A3m1/Nm1)
           )
           +(df22_main*rho_f2_2_main)*(
             # infected partner not on HCT/ART
             D1P1_f22_main*(H1m2/Nm2) + D1G1_f22_main*(Y1m2/Nm2) + D1X1_f22_main*(Z1m2/Nm2) + D1K1_f22_main*(A1m2/Nm2)
             # infected partner on HCT
             + D1P2_f22_main*(H2m2/Nm2) + D1G2_f22_main*(Y2m2/Nm2) + D1X2_f22_main*(Z2m2/Nm2) + D1K2_f22_main*(A2m2/Nm2)
             # infected partner on ART
             + D1P3_f22_main*(H3m2/Nm2) + D1G3_f22_main*(Y3m2/Nm2) + D1X3_f22_main*(Z3m2/Nm2) + D1K3_f22_main*(A3m2/Nm2)
           )
           +(df23_main*rho_f2_3_main)*(
             # infected partner not on HCT/ART
             D1P1_f23_main*(H1m3/Nm3) + D1G1_f23_main*(Y1m3/Nm3) + D1X1_f23_main*(Z1m3/Nm3) + D1K1_f23_main*(A1m3/Nm3)
             # infected partner on HCT
             + D1P2_f23_main*(H2m3/Nm3) + D1G2_f23_main*(Y2m3/Nm3) + D1X2_f23_main*(Z2m3/Nm3) + D1K2_f23_main*(A2m3/Nm3)
             # infected partner on ART
             + D1P3_f23_main*(H3m3/Nm3) + D1G3_f23_main*(Y3m3/Nm3) + D1X3_f23_main*(Z3m3/Nm3) + D1K3_f23_main*(A3m3/Nm3)
           )
           # regular
           +(df22_reg*rho_f2_2_reg)*(
             # infected partner not on HCT/ART
             D1P1_f22_reg*(H1m2/Nm2) + D1G1_f22_reg*(Y1m2/Nm2) + D1X1_f22_reg*(Z1m2/Nm2) + D1K1_f22_reg*(A1m2/Nm2)
             # infected partner on HCT
             + D1P2_f22_reg*(H2m2/Nm2) + D1G2_f22_reg*(Y2m2/Nm2) + D1X2_f22_reg*(Z2m2/Nm2) + D1K2_f22_reg*(A2m2/Nm2)
             # infected partner on ART
             + D1P3_f22_reg*(H3m2/Nm2) + D1G3_f22_reg*(Y3m2/Nm2) + D1X3_f22_reg*(Z3m2/Nm2) + D1K3_f22_reg*(A3m2/Nm2)
           )
           +(df23_reg*rho_f2_3_reg)*(
             # infected partner not on HCT/ART
             D1P1_f23_reg*(H1m3/Nm3) + D1G1_f23_reg*(Y1m3/Nm3) + D1X1_f23_reg*(Z1m3/Nm3) + D1K1_f23_reg*(A1m3/Nm3)
             # infected partner on HCT
             + D1P2_f23_reg*(H2m3/Nm3) + D1G2_f23_reg*(Y2m3/Nm3) + D1X2_f23_reg*(Z2m3/Nm3) + D1K2_f23_reg*(A2m3/Nm3)
             # infected partner on ART
             + D1P3_f23_reg*(H3m3/Nm3) + D1G3_f23_reg*(Y3m3/Nm3) + D1X3_f23_reg*(Z3m3/Nm3) + D1K3_f23_reg*(A3m3/Nm3)
           )
         )
         
         # pitest = (
         #   (1 - (((D1P1_f21_main*(H1m1/Nm1))^(df21_main*rho_f2_1_main)) * ((D1G1_f21_main*(Y1m1/Nm1))^(df21_main*rho_f2_1_main)) * ((D1X1_f21_main*(Z1m1/Nm1))^(df21_main*rho_f2_1_main)) * ((D1K1_f11_main*(A1m1/Nm1))^(df21_main*rho_f2_1_main))))
         # )
         
         pi1_f3 = (
           # main
           (df31_main*rho_f3_1_main)*(
             # infected partner not on HCT/ART
             D1P1_f31_main*(H1m1/Nm1) + D1G1_f31_main*(Y1m1/Nm1) + D1X1_f31_main*(Z1m1/Nm1) + D1K1_f31_main*(A1m1/Nm1)
             # infected partner on HCT
             + D1P2_f31_main*(H2m1/Nm1) + D1G2_f31_main*(Y2m1/Nm1) + D1X2_f31_main*(Z2m1/Nm1) + D1K2_f31_main*(A2m1/Nm1)
             # infected partner on ART
             + D1P3_f31_main*(H3m1/Nm1) + D1G3_f31_main*(Y3m1/Nm1) + D1X3_f31_main*(Z3m1/Nm1) + D1K3_f31_main*(A3m1/Nm1)
           )
           +(df32_main*rho_f3_2_main)*(
             # infected partner not on HCT/ART
             D1P1_f32_main*(H1m2/Nm2) + D1G1_f32_main*(Y1m2/Nm2) + D1X1_f32_main*(Z1m2/Nm2) + D1K1_f32_main*(A1m2/Nm2)
             # infected partner on HCT
             + D1P2_f32_main*(H2m2/Nm2) + D1G2_f32_main*(Y2m2/Nm2) + D1X2_f32_main*(Z2m2/Nm2) + D1K2_f32_main*(A2m2/Nm2)
             # infected partner on ART
             + D1P3_f32_main*(H3m2/Nm2) + D1G3_f32_main*(Y3m2/Nm2) + D1X3_f32_main*(Z3m2/Nm2) + D1K3_f32_main*(A3m2/Nm2)
           )
           +(df33_main*rho_f3_3_main)*(
             # infected partner not on HCT/ART
             D1P1_f33_main*(H1m3/Nm3) + D1G1_f33_main*(Y1m3/Nm3) + D1X1_f33_main*(Z1m3/Nm3) + D1K1_f33_main*(A1m3/Nm3)
             # infected partner on HCT
             + D1P2_f33_main*(H2m3/Nm3) + D1G2_f33_main*(Y2m3/Nm3) + D1X2_f33_main*(Z2m3/Nm3) + D1K2_f33_main*(A2m3/Nm3)
             # infected partner on ART
             + D1P3_f33_main*(H3m3/Nm3) + D1G3_f33_main*(Y3m3/Nm3) + D1X3_f33_main*(Z3m3/Nm3) + D1K3_f33_main*(A3m3/Nm3)
           )
           # regular
           +(df32_reg*rho_f3_2_reg)*(
             # infected partner not on HCT/ART
             D1P1_f32_reg*(H1m2/Nm2) + D1G1_f32_reg*(Y1m2/Nm2) + D1X1_f32_reg*(Z1m2/Nm2) + D1K1_f32_reg*(A1m2/Nm2)
             # infected partner on HCT
             + D1P2_f32_reg*(H2m2/Nm2) + D1G2_f32_reg*(Y2m2/Nm2) + D1X2_f32_reg*(Z2m2/Nm2) + D1K2_f32_reg*(A2m2/Nm2)
             # infected partner on ART
             + D1P3_f32_reg*(H3m2/Nm2) + D1G3_f32_reg*(Y3m2/Nm2) + D1X3_f32_reg*(Z3m2/Nm2) + D1K3_f32_reg*(A3m2/Nm2)
           )
           +(df33_reg*rho_f3_3_reg)*(
             # infected partner not on HCT/ART
             D1P1_f33_reg*(H1m3/Nm3) + D1G1_f33_reg*(Y1m3/Nm3) + D1X1_f33_reg*(Z1m3/Nm3) + D1K1_f33_reg*(A1m3/Nm3)
             # infected partner on HCT
             + D1P2_f33_reg*(H2m3/Nm3) + D1G2_f33_reg*(Y2m3/Nm3) + D1X2_f33_reg*(Z2m3/Nm3) + D1K2_f33_reg*(A2m3/Nm3)
             # infected partner on ART
             + D1P3_f33_reg*(H3m3/Nm3) + D1G3_f33_reg*(Y3m3/Nm3) + D1X3_f33_reg*(Z3m3/Nm3) + D1K3_f33_reg*(A3m3/Nm3)
           )
           # casual
           +(df33_casual*rho_f3_3_casual)*(
             # infected partner not on HCT/ART
             D1P1_f33_casual*(H1m3/Nm3) + D1G1_f33_casual*(Y1m3/Nm3) + D1X1_f33_casual*(Z1m3/Nm3) + D1K1_f33_casual*(A1m3/Nm3)
             # infected partner on HCT
             + D1P2_f33_casual*(H2m3/Nm3) + D1G2_f33_casual*(Y2m3/Nm3) + D1X2_f33_casual*(Z2m3/Nm3) + D1K2_f33_casual*(A2m3/Nm3)
             # infected partner on ART
             + D1P3_f33_casual*(H3m3/Nm3) + D1G3_f33_casual*(Y3m3/Nm3) + D1X3_f33_casual*(Z3m3/Nm3) + D1K3_f33_casual*(A3m3/Nm3)
           )
         )
         
         
         
         # Uncircumcised not on PrEP
         pi1_m1 = (
           # low-low
           (dm11_main*rho_m1_1_main)*(
             # infected partner not on HCT/ART
             D1P1_m11_main*(H1f1/Nf1) + D1G1_m11_main*(Y1f1/Nf1) + D1X1_m11_main*(Z1f1/Nf1) + D1K1_m11_main*(A1f1/Nf1)
             # infected partner on HCT
             + D1P2_m11_main*(H2f1/Nf1) + D1G2_m11_main*(Y2f1/Nf1) + D1X2_m11_main*(Z2f1/Nf1) + D1K2_m11_main*(A2f1/Nf1)
             # infected partner on ART
             + D1P3_m11_main*(H3f1/Nf1) + D1G3_m11_main*(Y3f1/Nf1) + D1X3_m11_main*(Z3f1/Nf1) + D1K3_m11_main*(A3f1/Nf1)
           )
           # low-medium
           +(dm12_main*rho_m1_2_main)*(
             # infected partner not on HCT/ART
             D1P1_m12_main*(H1f2/Nf2) + D1G1_m12_main*(Y1f2/Nf2) + D1X1_m12_main*(Z1f2/Nf2) + D1K1_m12_main*(A1f2/Nf2)
             # infected partner on HCT
             + D1P2_m12_main*(H2f2/Nf2) + D1G2_m12_main*(Y2f2/Nf2) + D1X2_m12_main*(Z2f2/Nf2) + D1K2_m12_main*(A2f2/Nf2)
             # infected partner on ART
             + D1P3_m12_main*(H3f2/Nf2) + D1G3_m12_main*(Y3f2/Nf2) + D1X3_m12_main*(Z3f2/Nf2) + D1K3_m12_main*(A3f2/Nf2)
           )
           # low-high
           +(dm13_main*rho_m1_3_main)*(
             # infected partner not on HCT/ART
             D1P1_m13_main*(H1f3/Nf3) + D1G1_m13_main*(Y1f3/Nf3) + D1X1_m13_main*(Z1f3/Nf3) + D1K1_m13_main*(A1f3/Nf3)
             # infected partner on HCT
             + D1P2_m13_main*(H2f3/Nf3) + D1G2_m13_main*(Y2f3/Nf3) + D1X2_m13_main*(Z2f3/Nf3) + D1K2_m13_main*(A2f3/Nf3)
             # infected partner on ART
             + D1P3_m13_main*(H3f3/Nf3) + D1G3_m13_main*(Y3f3/Nf3) + D1X3_m13_main*(Z3f3/Nf3) + D1K3_m13_main*(A3f3/Nf3)
           )
         )
         
         pi1_m2 = (
           # main
           (dm21_main*rho_m2_1_main)*(
             # infected partner not on HCT/ART
             D1P1_m21_main*(H1f1/Nf1) + D1G1_m21_main*(Y1f1/Nf1) + D1X1_m21_main*(Z1f1/Nf1) + D1K1_m21_main*(A1f1/Nf1)
             # infected partner on HCT
             + D1P2_m21_main*(H2f1/Nf1) + D1G2_m21_main*(Y2f1/Nf1) + D1X2_m21_main*(Z2f1/Nf1) + D1K2_m21_main*(A2f1/Nf1)
             # infected partner on ART
             + D1P3_m21_main*(H3f1/Nf1) + D1G3_m21_main*(Y3f1/Nf1) + D1X3_m21_main*(Z3f1/Nf1) + D1K3_m21_main*(A3f1/Nf1)
           )
           +(dm22_main*rho_m2_2_main)*(
             # infected partner not on HCT/ART
             D1P1_m22_main*(H1f2/Nf2) + D1G1_m22_main*(Y1f2/Nf2) + D1X1_m22_main*(Z1f2/Nf2) + D1K1_m22_main*(A1f2/Nf2)
             # infected partner on HCT
             + D1P2_m22_main*(H2f2/Nf2) + D1G2_m22_main*(Y2f2/Nf2) + D1X2_m22_main*(Z2f2/Nf2) + D1K2_m22_main*(A2f2/Nf2)
             # infected partner on ART
             + D1P3_m22_main*(H3f2/Nf2) + D1G3_m22_main*(Y3f2/Nf2) + D1X3_m22_main*(Z3f2/Nf2) + D1K3_m22_main*(A3f2/Nf2)
           )
           +(dm23_main*rho_m2_3_main)*(
             # infected partner not on HCT/ART
             D1P1_m23_main*(H1f3/Nf3) + D1G1_m23_main*(Y1f3/Nf3) + D1X1_m23_main*(Z1f3/Nf3) + D1K1_m23_main*(A1f3/Nf3)
             # infected partner on HCT
             + D1P2_m23_main*(H2f3/Nf3) + D1G2_m23_main*(Y2f3/Nf3) + D1X2_m23_main*(Z2f3/Nf3) + D1K2_m23_main*(A2f3/Nf3)
             # infected partner on ART
             + D1P3_m23_main*(H3f3/Nf3) + D1G3_m23_main*(Y3f3/Nf3) + D1X3_m23_main*(Z3f3/Nf3) + D1K3_m23_main*(A3f3/Nf3)
           )
           # regular
           +(dm22_reg*rho_m2_2_reg)*(
             # infected partner not on HCT/ART
             D1P1_m22_reg*(H1f2/Nf2) + D1G1_m22_reg*(Y1f2/Nf2) + D1X1_m22_reg*(Z1f2/Nf2) + D1K1_m22_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D1P2_m22_reg*(H2f2/Nf2) + D1G2_m22_reg*(Y2f2/Nf2) + D1X2_m22_reg*(Z2f2/Nf2) + D1K2_m22_reg*(A2f2/Nf2)
             # infected partner on ART
             + D1P3_m22_reg*(H3f2/Nf2) + D1G3_m22_reg*(Y3f2/Nf2) + D1X3_m22_reg*(Z3f2/Nf2) + D1K3_m22_reg*(A3f2/Nf2)
           )
           +(dm23_reg*rho_m2_3_reg)*(
             # infected partner not on HCT/ART
             D1P1_m23_reg*(H1f3/Nf3) + D1G1_m23_reg*(Y1f3/Nf3) + D1X1_m23_reg*(Z1f3/Nf3) + D1K1_m23_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D1P2_m23_reg*(H2f3/Nf3) + D1G2_m23_reg*(Y2f3/Nf3) + D1X2_m23_reg*(Z2f3/Nf3) + D1K2_m23_reg*(A2f3/Nf3)
             # infected partner on ART
             + D1P3_m23_reg*(H3f3/Nf3) + D1G3_m23_reg*(Y3f3/Nf3) + D1X3_m23_reg*(Z3f3/Nf3) + D1K3_m23_reg*(A3f3/Nf3)
           )
         )
         
         pi1_m3 = (
           # main
           (dm31_main*rho_m3_1_main)*(
             # infected partner not on HCT/ART
             D1P1_m31_main*(H1f1/Nf1) + D1G1_m31_main*(Y1f1/Nf1) + D1X1_m31_main*(Z1f1/Nf1) + D1K1_m31_main*(A1f1/Nf1)
             # infected partner on HCT
             + D1P2_m31_main*(H2f1/Nf1) + D1G2_m31_main*(Y2f1/Nf1) + D1X2_m31_main*(Z2f1/Nf1) + D1K2_m31_main*(A2f1/Nf1)
             # infected partner on ART
             + D1P3_m31_main*(H3f1/Nf1) + D1G3_m31_main*(Y3f1/Nf1) + D1X3_m31_main*(Z3f1/Nf1) + D1K3_m31_main*(A3f1/Nf1)
           )
           +(dm32_main*rho_m3_2_main)*(
             # infected partner not on HCT/ART
             D1P1_m32_main*(H1f2/Nf2) + D1G1_m32_main*(Y1f2/Nf2) + D1X1_m32_main*(Z1f2/Nf2) + D1K1_m32_main*(A1f2/Nf2)
             # infected partner on HCT
             + D1P2_m32_main*(H2f2/Nf2) + D1G2_m32_main*(Y2f2/Nf2) + D1X2_m32_main*(Z2f2/Nf2) + D1K2_m32_main*(A2f2/Nf2)
             # infected partner on ART
             + D1P3_m32_main*(H3f2/Nf2) + D1G3_m32_main*(Y3f2/Nf2) + D1X3_m32_main*(Z3f2/Nf2) + D1K3_m32_main*(A3f2/Nf2)
           )
           +(dm33_main*rho_m3_3_main)*(
             # infected partner not on HCT/ART
             D1P1_m33_main*(H1f3/Nf3) + D1G1_m33_main*(Y1f3/Nf3) + D1X1_m33_main*(Z1f3/Nf3) + D1K1_m33_main*(A1f3/Nf3)
             # infected partner on HCT
             + D1P2_m33_main*(H2f3/Nf3) + D1G2_m33_main*(Y2f3/Nf3) + D1X2_m33_main*(Z2f3/Nf3) + D1K2_m33_main*(A2f3/Nf3)
             # infected partner on ART
             + D1P3_m33_main*(H3f3/Nf3) + D1G3_m33_main*(Y3f3/Nf3) + D1X3_m33_main*(Z3f3/Nf3) + D1K3_m33_main*(A3f3/Nf3)
           )
           # regular
           +(dm32_reg*rho_m3_2_reg)*(
             # infected partner not on HCT/ART
             D1P1_m32_reg*(H1f2/Nf2) + D1G1_m32_reg*(Y1f2/Nf2) + D1X1_m32_reg*(Z1f2/Nf2) + D1K1_m32_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D1P2_m32_reg*(H2f2/Nf2) + D1G2_m32_reg*(Y2f2/Nf2) + D1X2_m32_reg*(Z2f2/Nf2) + D1K2_m32_reg*(A2f2/Nf2)
             # infected partner on ART
             + D1P3_m32_reg*(H3f2/Nf2) + D1G3_m32_reg*(Y3f2/Nf2) + D1X3_m32_reg*(Z3f2/Nf2) + D1K3_m32_reg*(A3f2/Nf2)
           )
           +(dm33_reg*rho_m3_3_reg)*(
             # infected partner not on HCT/ART
             D1P1_m33_reg*(H1f3/Nf3) + D1G1_m33_reg*(Y1f3/Nf3) + D1X1_m33_reg*(Z1f3/Nf3) + D1K1_m33_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D1P2_m33_reg*(H2f3/Nf3) + D1G2_m33_reg*(Y2f3/Nf3) + D1X2_m33_reg*(Z2f3/Nf3) + D1K2_m33_reg*(A2f3/Nf3)
             # infected partner on ART
             + D1P3_m33_reg*(H3f3/Nf3) + D1G3_m33_reg*(Y3f3/Nf3) + D1X3_m33_reg*(Z3f3/Nf3) + D1K3_m33_reg*(A3f3/Nf3)
           )
           # casual
           +(dm33_casual*rho_m3_3_casual)*(
             # infected partner not on HCT/ART
             D1P1_m33_casual*(H1f3/Nf3) + D1G1_m33_casual*(Y1f3/Nf3) + D1X1_m33_casual*(Z1f3/Nf3) + D1K1_m33_casual*(A1f3/Nf3)
             # infected partner on HCT
             + D1P2_m33_casual*(H2f3/Nf3) + D1G2_m33_casual*(Y2f3/Nf3) + D1X2_m33_casual*(Z2f3/Nf3) + D1K2_m33_casual*(A2f3/Nf3)
             # infected partner on ART
             + D1P3_m33_casual*(H3f3/Nf3) + D1G3_m33_casual*(Y3f3/Nf3) + D1X3_m33_casual*(Z3f3/Nf3) + D1K3_m33_casual*(A3f3/Nf3)
           )
         )
         
         
         # Females on PrEP
         pi2_f1 = (
           # low-low
           (df11_main*rho_f1_1_main)*(
             # infected partner not on HCT/ART
             D2P1_f11_main*(H1m1/Nm1) + D2G1_f11_main*(Y1m1/Nm1) + D2X1_f11_main*(Z1m1/Nm1) + D2K1_f11_main*(A1m1/Nm1)
             # infected partner on HCT
             + D2P2_f11_main*(H2m1/Nm1) + D2G2_f11_main*(Y2m1/Nm1) + D2X2_f11_main*(Z2m1/Nm1) + D2K2_f11_main*(A2m1/Nm1)
             # infected partner on ART
             + D2P3_f11_main*(H3m1/Nm1) + D2G3_f11_main*(Y3m1/Nm1) + D2X3_f11_main*(Z3m1/Nm1) + D2K3_f11_main*(A3m1/Nm1)
           )
           # low-medium
           +(df12_main*rho_f1_2_main)*(
             # infected partner not on HCT/ART
             D2P1_f12_main*(H1m2/Nm2) + D2G1_f12_main*(Y1m2/Nm2) + D2X1_f12_main*(Z1m2/Nm2) + D2K1_f12_main*(A1m2/Nm2)
             # infected partner on HCT
             + D2P2_f12_main*(H2m2/Nm2) + D2G2_f12_main*(Y2m2/Nm2) + D2X2_f12_main*(Z2m2/Nm2) + D2K2_f12_main*(A2m2/Nm2)
             # infected partner on ART
             + D2P3_f12_main*(H3m2/Nm2) + D2G3_f12_main*(Y3m2/Nm2) + D2X3_f12_main*(Z3m2/Nm2) + D2K3_f12_main*(A3m2/Nm2)
           )
           # low-high
           +(df13_main*rho_f1_3_main)*(
             # infected partner not on HCT/ART
             D2P1_f13_main*(H1m3/Nm3) + D2G1_f13_main*(Y1m3/Nm3) + D2X1_f13_main*(Z1m3/Nm3) + D2K1_f13_main*(A1m3/Nm3)
             # infected partner on HCT
             + D2P2_f13_main*(H2m3/Nm3) + D2G2_f13_main*(Y2m3/Nm3) + D2X2_f13_main*(Z2m3/Nm3) + D2K2_f13_main*(A2m3/Nm3)
             # infected partner on ART
             + D2P3_f13_main*(H3m3/Nm3) + D2G3_f13_main*(Y3m3/Nm3) + D2X3_f13_main*(Z3m3/Nm3) + D2K3_f13_main*(A3m3/Nm3)
           )
         )
         
         pi2_f2 = (
           # main
           (df21_main*rho_f2_1_main)*(
             # infected partner not on HCT/ART
             D2P1_f21_main*(H1m1/Nm1) + D2G1_f21_main*(Y1m1/Nm1) + D2X1_f21_main*(Z1m1/Nm1) + D2K1_f21_main*(A1m1/Nm1)
             # infected partner on HCT
             + D2P2_f21_main*(H2m1/Nm1) + D2G2_f21_main*(Y2m1/Nm1) + D2X2_f21_main*(Z2m1/Nm1) + D2K2_f21_main*(A2m1/Nm1)
             # infected partner on ART
             + D2P3_f21_main*(H3m1/Nm1) + D2G3_f21_main*(Y3m1/Nm1) + D2X3_f21_main*(Z3m1/Nm1) + D2K3_f21_main*(A3m1/Nm1)
           )
           +(df22_main*rho_f2_2_main)*(
             # infected partner not on HCT/ART
             D2P1_f22_main*(H1m2/Nm2) + D2G1_f22_main*(Y1m2/Nm2) + D2X1_f22_main*(Z1m2/Nm2) + D2K1_f22_main*(A1m2/Nm2)
             # infected partner on HCT
             + D2P2_f22_main*(H2m2/Nm2) + D2G2_f22_main*(Y2m2/Nm2) + D2X2_f22_main*(Z2m2/Nm2) + D2K2_f22_main*(A2m2/Nm2)
             # infected partner on ART
             + D2P3_f22_main*(H3m2/Nm2) + D2G3_f22_main*(Y3m2/Nm2) + D2X3_f22_main*(Z3m2/Nm2) + D2K3_f22_main*(A3m2/Nm2)
           )
           +(df23_main*rho_f2_3_main)*(
             # infected partner not on HCT/ART
             D2P1_f23_main*(H1m3/Nm3) + D2G1_f23_main*(Y1m3/Nm3) + D2X1_f23_main*(Z1m3/Nm3) + D2K1_f23_main*(A1m3/Nm3)
             # infected partner on HCT
             + D2P2_f23_main*(H2m3/Nm3) + D2G2_f23_main*(Y2m3/Nm3) + D2X2_f23_main*(Z2m3/Nm3) + D2K2_f23_main*(A2m3/Nm3)
             # infected partner on ART
             + D2P3_f23_main*(H3m3/Nm3) + D2G3_f23_main*(Y3m3/Nm3) + D2X3_f23_main*(Z3m3/Nm3) + D2K3_f23_main*(A3m3/Nm3)
           )
           # regular
           +(df22_reg*rho_f2_2_reg)*(
             # infected partner not on HCT/ART
             D2P1_f22_reg*(H1m2/Nm2) + D2G1_f22_reg*(Y1m2/Nm2) + D2X1_f22_reg*(Z1m2/Nm2) + D2K1_f22_reg*(A1m2/Nm2)
             # infected partner on HCT
             + D2P2_f22_reg*(H2m2/Nm2) + D2G2_f22_reg*(Y2m2/Nm2) + D2X2_f22_reg*(Z2m2/Nm2) + D2K2_f22_reg*(A2m2/Nm2)
             # infected partner on ART
             + D2P3_f22_reg*(H3m2/Nm2) + D2G3_f22_reg*(Y3m2/Nm2) + D2X3_f22_reg*(Z3m2/Nm2) + D2K3_f22_reg*(A3m2/Nm2)
           )
           +(df23_reg*rho_f2_3_reg)*(
             # infected partner not on HCT/ART
             D2P1_f23_reg*(H1m3/Nm3) + D2G1_f23_reg*(Y1m3/Nm3) + D2X1_f23_reg*(Z1m3/Nm3) + D2K1_f23_reg*(A1m3/Nm3)
             # infected partner on HCT
             + D2P2_f23_reg*(H2m3/Nm3) + D2G2_f23_reg*(Y2m3/Nm3) + D2X2_f23_reg*(Z2m3/Nm3) + D2K2_f23_reg*(A2m3/Nm3)
             # infected partner on ART
             + D2P3_f23_reg*(H3m3/Nm3) + D2G3_f23_reg*(Y3m3/Nm3) + D2X3_f23_reg*(Z3m3/Nm3) + D2K3_f23_reg*(A3m3/Nm3)
           )
         )
         
         pi2_f3 = (
           # main
           (df31_main*rho_f3_1_main)*(
             # infected partner not on HCT/ART
             D2P1_f31_main*(H1m1/Nm1) + D2G1_f31_main*(Y1m1/Nm1) + D2X1_f31_main*(Z1m1/Nm1) + D2K1_f31_main*(A1m1/Nm1)
             # infected partner on HCT
             + D2P2_f31_main*(H2m1/Nm1) + D2G2_f31_main*(Y2m1/Nm1) + D2X2_f31_main*(Z2m1/Nm1) + D2K2_f31_main*(A2m1/Nm1)
             # infected partner on ART
             + D2P3_f31_main*(H3m1/Nm1) + D2G3_f31_main*(Y3m1/Nm1) + D2X3_f31_main*(Z3m1/Nm1) + D2K3_f31_main*(A3m1/Nm1)
           )
           +(df32_main*rho_f3_2_main)*(
             # infected partner not on HCT/ART
             D2P1_f32_main*(H1m2/Nm2) + D2G1_f32_main*(Y1m2/Nm2) + D2X1_f32_main*(Z1m2/Nm2) + D2K1_f32_main*(A1m2/Nm2)
             # infected partner on HCT
             + D2P2_f32_main*(H2m2/Nm2) + D2G2_f32_main*(Y2m2/Nm2) + D2X2_f32_main*(Z2m2/Nm2) + D2K2_f32_main*(A2m2/Nm2)
             # infected partner on ART
             + D2P3_f32_main*(H3m2/Nm2) + D2G3_f32_main*(Y3m2/Nm2) + D2X3_f32_main*(Z3m2/Nm2) + D2K3_f32_main*(A3m2/Nm2)
           )
           +(df33_main*rho_f3_3_main)*(
             # infected partner not on HCT/ART
             D2P1_f33_main*(H1m3/Nm3) + D2G1_f33_main*(Y1m3/Nm3) + D2X1_f33_main*(Z1m3/Nm3) + D2K1_f33_main*(A1m3/Nm3)
             # infected partner on HCT
             + D2P2_f33_main*(H2m3/Nm3) + D2G2_f33_main*(Y2m3/Nm3) + D2X2_f33_main*(Z2m3/Nm3) + D2K2_f33_main*(A2m3/Nm3)
             # infected partner on ART
             + D2P3_f33_main*(H3m3/Nm3) + D2G3_f33_main*(Y3m3/Nm3) + D2X3_f33_main*(Z3m3/Nm3) + D2K3_f33_main*(A3m3/Nm3)
           )
           # regular
           +(df32_reg*rho_f3_2_reg)*(
             # infected partner not on HCT/ART
             D2P1_f32_reg*(H1m2/Nm2) + D2G1_f32_reg*(Y1m2/Nm2) + D2X1_f32_reg*(Z1m2/Nm2) + D2K1_f32_reg*(A1m2/Nm2)
             # infected partner on HCT
             + D2P2_f32_reg*(H2m2/Nm2) + D2G2_f32_reg*(Y2m2/Nm2) + D2X2_f32_reg*(Z2m2/Nm2) + D2K2_f32_reg*(A2m2/Nm2)
             # infected partner on ART
             + D2P3_f32_reg*(H3m2/Nm2) + D2G3_f32_reg*(Y3m2/Nm2) + D2X3_f32_reg*(Z3m2/Nm2) + D2K3_f32_reg*(A3m2/Nm2)
           )
           +(df33_reg*rho_f3_3_reg)*(
             # infected partner not on HCT/ART
             D2P1_f33_reg*(H1m3/Nm3) + D2G1_f33_reg*(Y1m3/Nm3) + D2X1_f33_reg*(Z1m3/Nm3) + D2K1_f33_reg*(A1m3/Nm3)
             # infected partner on HCT
             + D2P2_f33_reg*(H2m3/Nm3) + D2G2_f33_reg*(Y2m3/Nm3) + D2X2_f33_reg*(Z2m3/Nm3) + D2K2_f33_reg*(A2m3/Nm3)
             # infected partner on ART
             + D2P3_f33_reg*(H3m3/Nm3) + D2G3_f33_reg*(Y3m3/Nm3) + D2X3_f33_reg*(Z3m3/Nm3) + D2K3_f33_reg*(A3m3/Nm3)
           )
           # casual
           +(df33_casual*rho_f3_3_casual)*(
             # infected partner not on HCT/ART
             D2P1_f33_casual*(H1m3/Nm3) + D2G1_f33_casual*(Y1m3/Nm3) + D2X1_f33_casual*(Z1m3/Nm3) + D2K1_f33_casual*(A1m3/Nm3)
             # infected partner on HCT
             + D2P2_f33_casual*(H2m3/Nm3) + D2G2_f33_casual*(Y2m3/Nm3) + D2X2_f33_casual*(Z2m3/Nm3) + D2K2_f33_casual*(A2m3/Nm3)
             # infected partner on ART
             + D2P3_f33_casual*(H3m3/Nm3) + D2G3_f33_casual*(Y3m3/Nm3) + D2X3_f33_casual*(Z3m3/Nm3) + D2K3_f33_casual*(A3m3/Nm3)
           )
         )
         
         ###
         # Uncircumcised males on PrEP
         pi2_m1 = (
           # low-low
           (dm11_main*rho_m1_1_main)*(
             # infected partner not on HCT/ART
             D2P1_m11_main*(H1f1/Nf1) + D2G1_m11_main*(Y1f1/Nf1) + D2X1_m11_main*(Z1f1/Nf1) + D2K1_m11_main*(A1f1/Nf1)
             # infected partner on HCT
             + D2P2_m11_main*(H2f1/Nf1) + D2G2_m11_main*(Y2f1/Nf1) + D2X2_m11_main*(Z2f1/Nf1) + D2K2_m11_main*(A2f1/Nf1)
             # infected partner on ART
             + D2P3_m11_main*(H3f1/Nf1) + D2G3_m11_main*(Y3f1/Nf1) + D2X3_m11_main*(Z3f1/Nf1) + D2K3_m11_main*(A3f1/Nf1)
           )
           # low-medium
           +(dm12_main*rho_m1_2_main)*(
             # infected partner not on HCT/ART
             D2P1_m12_main*(H1f2/Nf2) + D2G1_m12_main*(Y1f2/Nf2) + D2X1_m12_main*(Z1f2/Nf2) + D2K1_m12_main*(A1f2/Nf2)
             # infected partner on HCT
             + D2P2_m12_main*(H2f2/Nf2) + D2G2_m12_main*(Y2f2/Nf2) + D2X2_m12_main*(Z2f2/Nf2) + D2K2_m12_main*(A2f2/Nf2)
             # infected partner on ART
             + D2P3_m12_main*(H3f2/Nf2) + D2G3_m12_main*(Y3f2/Nf2) + D2X3_m12_main*(Z3f2/Nf2) + D2K3_m12_main*(A3f2/Nf2)
           )
           # low-high
           +(dm13_main*rho_m1_3_main)*(
             # infected partner not on HCT/ART
             D2P1_m13_main*(H1f3/Nf3) + D2G1_m13_main*(Y1f3/Nf3) + D2X1_m13_main*(Z1f3/Nf3) + D2K1_m13_main*(A1f3/Nf3)
             # infected partner on HCT
             + D2P2_m13_main*(H2f3/Nf3) + D2G2_m13_main*(Y2f3/Nf3) + D2X2_m13_main*(Z2f3/Nf3) + D2K2_m13_main*(A2f3/Nf3)
             # infected partner on ART
             + D2P3_m13_main*(H3f3/Nf3) + D2G3_m13_main*(Y3f3/Nf3) + D2X3_m13_main*(Z3f3/Nf3) + D2K3_m13_main*(A3f3/Nf3)
           )
         )
         
         pi2_m2 = (
           # main
           (dm21_main*rho_m2_1_main)*(
             # infected partner not on HCT/ART
             D2P1_m21_main*(H1f1/Nf1) + D2G1_m21_main*(Y1f1/Nf1) + D2X1_m21_main*(Z1f1/Nf1) + D2K1_m21_main*(A1f1/Nf1)
             # infected partner on HCT
             + D2P2_m21_main*(H2f1/Nf1) + D2G2_m21_main*(Y2f1/Nf1) + D2X2_m21_main*(Z2f1/Nf1) + D2K2_m21_main*(A2f1/Nf1)
             # infected partner on ART
             + D2P3_m21_main*(H3f1/Nf1) + D2G3_m21_main*(Y3f1/Nf1) + D2X3_m21_main*(Z3f1/Nf1) + D2K3_m21_main*(A3f1/Nf1)
           )
           +(dm22_main*rho_m2_2_main)*(
             # infected partner not on HCT/ART
             D2P1_m22_main*(H1f2/Nf2) + D2G1_m22_main*(Y1f2/Nf2) + D2X1_m22_main*(Z1f2/Nf2) + D2K1_m22_main*(A1f2/Nf2)
             # infected partner on HCT
             + D2P2_m22_main*(H2f2/Nf2) + D2G2_m22_main*(Y2f2/Nf2) + D2X2_m22_main*(Z2f2/Nf2) + D2K2_m22_main*(A2f2/Nf2)
             # infected partner on ART
             + D2P3_m22_main*(H3f2/Nf2) + D2G3_m22_main*(Y3f2/Nf2) + D2X3_m22_main*(Z3f2/Nf2) + D2K3_m22_main*(A3f2/Nf2)
           )
           +(dm23_main*rho_m2_3_main)*(
             # infected partner not on HCT/ART
             D2P1_m23_main*(H1f3/Nf3) + D2G1_m23_main*(Y1f3/Nf3) + D2X1_m23_main*(Z1f3/Nf3) + D2K1_m23_main*(A1f3/Nf3)
             # infected partner on HCT
             + D2P2_m23_main*(H2f3/Nf3) + D2G2_m23_main*(Y2f3/Nf3) + D2X2_m23_main*(Z2f3/Nf3) + D2K2_m23_main*(A2f3/Nf3)
             # infected partner on ART
             + D2P3_m23_main*(H3f3/Nf3) + D2G3_m23_main*(Y3f3/Nf3) + D2X3_m23_main*(Z3f3/Nf3) + D2K3_m23_main*(A3f3/Nf3)
           )
           # regular
           +(dm22_reg*rho_m2_2_reg)*(
             # infected partner not on HCT/ART
             D2P1_m22_reg*(H1f2/Nf2) + D2G1_m22_reg*(Y1f2/Nf2) + D2X1_m22_reg*(Z1f2/Nf2) + D2K1_m22_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D2P2_m22_reg*(H2f2/Nf2) + D2G2_m22_reg*(Y2f2/Nf2) + D2X2_m22_reg*(Z2f2/Nf2) + D2K2_m22_reg*(A2f2/Nf2)
             # infected partner on ART
             + D2P3_m22_reg*(H3f2/Nf2) + D2G3_m22_reg*(Y3f2/Nf2) + D2X3_m22_reg*(Z3f2/Nf2) + D2K3_m22_reg*(A3f2/Nf2)
           )
           +(dm23_reg*rho_m2_3_reg)*(
             # infected partner not on HCT/ART
             D2P1_m23_reg*(H1f3/Nf3) + D2G1_m23_reg*(Y1f3/Nf3) + D2X1_m23_reg*(Z1f3/Nf3) + D2K1_m23_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D2P2_m23_reg*(H2f3/Nf3) + D2G2_m23_reg*(Y2f3/Nf3) + D2X2_m23_reg*(Z2f3/Nf3) + D2K2_m23_reg*(A2f3/Nf3)
             # infected partner on ART
             + D2P3_m23_reg*(H3f3/Nf3) + D2G3_m23_reg*(Y3f3/Nf3) + D2X3_m23_reg*(Z3f3/Nf3) + D2K3_m23_reg*(A3f3/Nf3)
           )
         )
         
         pi2_m3 = (
           # main
           (dm31_main*rho_m3_1_main)*(
             # infected partner not on HCT/ART
             D2P1_m31_main*(H1f1/Nf1) + D2G1_m31_main*(Y1f1/Nf1) + D2X1_m31_main*(Z1f1/Nf1) + D2K1_m31_main*(A1f1/Nf1)
             # infected partner on HCT
             + D2P2_m31_main*(H2f1/Nf1) + D2G2_m31_main*(Y2f1/Nf1) + D2X2_m31_main*(Z2f1/Nf1) + D2K2_m31_main*(A2f1/Nf1)
             # infected partner on ART
             + D2P3_m31_main*(H3f1/Nf1) + D2G3_m31_main*(Y3f1/Nf1) + D2X3_m31_main*(Z3f1/Nf1) + D2K3_m31_main*(A3f1/Nf1)
           )
           +(dm32_main*rho_m3_2_main)*(
             # infected partner not on HCT/ART
             D2P1_m32_main*(H1f2/Nf2) + D2G1_m32_main*(Y1f2/Nf2) + D2X1_m32_main*(Z1f2/Nf2) + D2K1_m32_main*(A1f2/Nf2)
             # infected partner on HCT
             + D2P2_m32_main*(H2f2/Nf2) + D2G2_m32_main*(Y2f2/Nf2) + D2X2_m32_main*(Z2f2/Nf2) + D2K2_m32_main*(A2f2/Nf2)
             # infected partner on ART
             + D2P3_m32_main*(H3f2/Nf2) + D2G3_m32_main*(Y3f2/Nf2) + D2X3_m32_main*(Z3f2/Nf2) + D2K3_m32_main*(A3f2/Nf2)
           )
           +(dm33_main*rho_m3_3_main)*(
             # infected partner not on HCT/ART
             D2P1_m33_main*(H1f3/Nf3) + D2G1_m33_main*(Y1f3/Nf3) + D2X1_m33_main*(Z1f3/Nf3) + D2K1_m33_main*(A1f3/Nf3)
             # infected partner on HCT
             + D2P2_m33_main*(H2f3/Nf3) + D2G2_m33_main*(Y2f3/Nf3) + D2X2_m33_main*(Z2f3/Nf3) + D2K2_m33_main*(A2f3/Nf3)
             # infected partner on ART
             + D2P3_m33_main*(H3f3/Nf3) + D2G3_m33_main*(Y3f3/Nf3) + D2X3_m33_main*(Z3f3/Nf3) + D2K3_m33_main*(A3f3/Nf3)
           )
           # regular
           +(dm32_reg*rho_m3_2_reg)*(
             # infected partner not on HCT/ART
             D2P1_m32_reg*(H1f2/Nf2) + D2G1_m32_reg*(Y1f2/Nf2) + D2X1_m32_reg*(Z1f2/Nf2) + D2K1_m32_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D2P2_m32_reg*(H2f2/Nf2) + D2G2_m32_reg*(Y2f2/Nf2) + D2X2_m32_reg*(Z2f2/Nf2) + D2K2_m32_reg*(A2f2/Nf2)
             # infected partner on ART
             + D2P3_m32_reg*(H3f2/Nf2) + D2G3_m32_reg*(Y3f2/Nf2) + D2X3_m32_reg*(Z3f2/Nf2) + D2K3_m32_reg*(A3f2/Nf2)
           )
           +(dm33_reg*rho_m3_3_reg)*(
             # infected partner not on HCT/ART
             D2P1_m33_reg*(H1f3/Nf3) + D2G1_m33_reg*(Y1f3/Nf3) + D2X1_m33_reg*(Z1f3/Nf3) + D2K1_m33_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D2P2_m33_reg*(H2f3/Nf3) + D2G2_m33_reg*(Y2f3/Nf3) + D2X2_m33_reg*(Z2f3/Nf3) + D2K2_m33_reg*(A2f3/Nf3)
             # infected partner on ART
             + D2P3_m33_reg*(H3f3/Nf3) + D2G3_m33_reg*(Y3f3/Nf3) + D2X3_m33_reg*(Z3f3/Nf3) + D2K3_m33_reg*(A3f3/Nf3)
           )
           # casual
           +(dm33_casual*rho_m3_3_casual)*(
             # infected partner not on HCT/ART
             D2P1_m33_casual*(H1f3/Nf3) + D2G1_m33_casual*(Y1f3/Nf3) + D2X1_m33_casual*(Z1f3/Nf3) + D2K1_m33_casual*(A1f3/Nf3)
             # infected partner on HCT
             + D2P2_m33_casual*(H2f3/Nf3) + D2G2_m33_casual*(Y2f3/Nf3) + D2X2_m33_casual*(Z2f3/Nf3) + D2K2_m33_casual*(A2f3/Nf3)
             # infected partner on ART
             + D2P3_m33_casual*(H3f3/Nf3) + D2G3_m33_casual*(Y3f3/Nf3) + D2X3_m33_casual*(Z3f3/Nf3) + D2K3_m33_casual*(A3f3/Nf3)
           )
         )
         
         ###
         # Circumcised males not on PrEP
         pi3_m1 = (
           # low-low
           (dm11_main*rho_m1_1_main)*(
             # infected partner not on HCT/ART
             D3P1_m11_main*(H1f1/Nf1) + D3G1_m11_main*(Y1f1/Nf1) + D3X1_m11_main*(Z1f1/Nf1) + D3K1_m11_main*(A1f1/Nf1)
             # infected partner on HCT
             + D3P2_m11_main*(H2f1/Nf1) + D3G2_m11_main*(Y2f1/Nf1) + D3X2_m11_main*(Z2f1/Nf1) + D3K2_m11_main*(A2f1/Nf1)
             # infected partner on ART
             + D3P3_m11_main*(H3f1/Nf1) + D3G3_m11_main*(Y3f1/Nf1) + D3X3_m11_main*(Z3f1/Nf1) + D3K3_m11_main*(A3f1/Nf1)
           )
           # low-medium
           +(dm12_main*rho_m1_2_main)*(
             # infected partner not on HCT/ART
             D3P1_m12_main*(H1f2/Nf2) + D3G1_m12_main*(Y1f2/Nf2) + D3X1_m12_main*(Z1f2/Nf2) + D3K1_m12_main*(A1f2/Nf2)
             # infected partner on HCT
             + D3P2_m12_main*(H2f2/Nf2) + D3G2_m12_main*(Y2f2/Nf2) + D3X2_m12_main*(Z2f2/Nf2) + D3K2_m12_main*(A2f2/Nf2)
             # infected partner on ART
             + D3P3_m12_main*(H3f2/Nf2) + D3G3_m12_main*(Y3f2/Nf2) + D3X3_m12_main*(Z3f2/Nf2) + D3K3_m12_main*(A3f2/Nf2)
           )
           # low-high
           +(dm13_main*rho_m1_3_main)*(
             # infected partner not on HCT/ART
             D3P1_m13_main*(H1f3/Nf3) + D3G1_m13_main*(Y1f3/Nf3) + D3X1_m13_main*(Z1f3/Nf3) + D3K1_m13_main*(A1f3/Nf3)
             # infected partner on HCT
             + D3P2_m13_main*(H2f3/Nf3) + D3G2_m13_main*(Y2f3/Nf3) + D3X2_m13_main*(Z2f3/Nf3) + D3K2_m13_main*(A2f3/Nf3)
             # infected partner on ART
             + D3P3_m13_main*(H3f3/Nf3) + D3G3_m13_main*(Y3f3/Nf3) + D3X3_m13_main*(Z3f3/Nf3) + D3K3_m13_main*(A3f3/Nf3)
           )
         )
         
         pi3_m2 = (
           # main
           (dm21_main*rho_m2_1_main)*(
             # infected partner not on HCT/ART
             D3P1_m21_main*(H1f1/Nf1) + D3G1_m21_main*(Y1f1/Nf1) + D3X1_m21_main*(Z1f1/Nf1) + D3K1_m21_main*(A1f1/Nf1)
             # infected partner on HCT
             + D3P2_m21_main*(H2f1/Nf1) + D3G2_m21_main*(Y2f1/Nf1) + D3X2_m21_main*(Z2f1/Nf1) + D3K2_m21_main*(A2f1/Nf1)
             # infected partner on ART
             + D3P3_m21_main*(H3f1/Nf1) + D3G3_m21_main*(Y3f1/Nf1) + D3X3_m21_main*(Z3f1/Nf1) + D3K3_m21_main*(A3f1/Nf1)
           )
           +(dm22_main*rho_m2_2_main)*(
             # infected partner not on HCT/ART
             D3P1_m22_main*(H1f2/Nf2) + D3G1_m22_main*(Y1f2/Nf2) + D3X1_m22_main*(Z1f2/Nf2) + D3K1_m22_main*(A1f2/Nf2)
             # infected partner on HCT
             + D3P2_m22_main*(H2f2/Nf2) + D3G2_m22_main*(Y2f2/Nf2) + D3X2_m22_main*(Z2f2/Nf2) + D3K2_m22_main*(A2f2/Nf2)
             # infected partner on ART
             + D3P3_m22_main*(H3f2/Nf2) + D3G3_m22_main*(Y3f2/Nf2) + D3X3_m22_main*(Z3f2/Nf2) + D3K3_m22_main*(A3f2/Nf2)
           )
           +(dm23_main*rho_m2_3_main)*(
             # infected partner not on HCT/ART
             D3P1_m23_main*(H1f3/Nf3) + D3G1_m23_main*(Y1f3/Nf3) + D3X1_m23_main*(Z1f3/Nf3) + D3K1_m23_main*(A1f3/Nf3)
             # infected partner on HCT
             + D3P2_m23_main*(H2f3/Nf3) + D3G2_m23_main*(Y2f3/Nf3) + D3X2_m23_main*(Z2f3/Nf3) + D3K2_m23_main*(A2f3/Nf3)
             # infected partner on ART
             + D3P3_m23_main*(H3f3/Nf3) + D3G3_m23_main*(Y3f3/Nf3) + D3X3_m23_main*(Z3f3/Nf3) + D3K3_m23_main*(A3f3/Nf3)
           )
           # regular
           +(dm22_reg*rho_m2_2_reg)*(
             # infected partner not on HCT/ART
             D3P1_m22_reg*(H1f2/Nf2) + D3G1_m22_reg*(Y1f2/Nf2) + D3X1_m22_reg*(Z1f2/Nf2) + D3K1_m22_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D3P2_m22_reg*(H2f2/Nf2) + D3G2_m22_reg*(Y2f2/Nf2) + D3X2_m22_reg*(Z2f2/Nf2) + D3K2_m22_reg*(A2f2/Nf2)
             # infected partner on ART
             + D3P3_m22_reg*(H3f2/Nf2) + D3G3_m22_reg*(Y3f2/Nf2) + D3X3_m22_reg*(Z3f2/Nf2) + D3K3_m22_reg*(A3f2/Nf2)
           )
           +(dm23_reg*rho_m2_3_reg)*(
             # infected partner not on HCT/ART
             D3P1_m23_reg*(H1f3/Nf3) + D3G1_m23_reg*(Y1f3/Nf3) + D3X1_m23_reg*(Z1f3/Nf3) + D3K1_m23_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D3P2_m23_reg*(H2f3/Nf3) + D3G2_m23_reg*(Y2f3/Nf3) + D3X2_m23_reg*(Z2f3/Nf3) + D3K2_m23_reg*(A2f3/Nf3)
             # infected partner on ART
             + D3P3_m23_reg*(H3f3/Nf3) + D3G3_m23_reg*(Y3f3/Nf3) + D3X3_m23_reg*(Z3f3/Nf3) + D3K3_m23_reg*(A3f3/Nf3)
           )
         )
         
         pi3_m3 = (
           # main
           (dm31_main*rho_m3_1_main)*(
             # infected partner not on HCT/ART
             D3P1_m31_main*(H1f1/Nf1) + D3G1_m31_main*(Y1f1/Nf1) + D3X1_m31_main*(Z1f1/Nf1) + D3K1_m31_main*(A1f1/Nf1)
             # infected partner on HCT
             + D3P2_m31_main*(H2f1/Nf1) + D3G2_m31_main*(Y2f1/Nf1) + D3X2_m31_main*(Z2f1/Nf1) + D3K2_m31_main*(A2f1/Nf1)
             # infected partner on ART
             + D3P3_m31_main*(H3f1/Nf1) + D3G3_m31_main*(Y3f1/Nf1) + D3X3_m31_main*(Z3f1/Nf1) + D3K3_m31_main*(A3f1/Nf1)
           )
           +(dm32_main*rho_m3_2_main)*(
             # infected partner not on HCT/ART
             D3P1_m32_main*(H1f2/Nf2) + D3G1_m32_main*(Y1f2/Nf2) + D3X1_m32_main*(Z1f2/Nf2) + D3K1_m32_main*(A1f2/Nf2)
             # infected partner on HCT
             + D3P2_m32_main*(H2f2/Nf2) + D3G2_m32_main*(Y2f2/Nf2) + D3X2_m32_main*(Z2f2/Nf2) + D3K2_m32_main*(A2f2/Nf2)
             # infected partner on ART
             + D3P3_m32_main*(H3f2/Nf2) + D3G3_m32_main*(Y3f2/Nf2) + D3X3_m32_main*(Z3f2/Nf2) + D3K3_m32_main*(A3f2/Nf2)
           )
           +(dm33_main*rho_m3_3_main)*(
             # infected partner not on HCT/ART
             D3P1_m33_main*(H1f3/Nf3) + D3G1_m33_main*(Y1f3/Nf3) + D3X1_m33_main*(Z1f3/Nf3) + D3K1_m33_main*(A1f3/Nf3)
             # infected partner on HCT
             + D3P2_m33_main*(H2f3/Nf3) + D3G2_m33_main*(Y2f3/Nf3) + D3X2_m33_main*(Z2f3/Nf3) + D3K2_m33_main*(A2f3/Nf3)
             # infected partner on ART
             + D3P3_m33_main*(H3f3/Nf3) + D3G3_m33_main*(Y3f3/Nf3) + D3X3_m33_main*(Z3f3/Nf3) + D3K3_m33_main*(A3f3/Nf3)
           )
           # regular
           +(dm32_reg*rho_m3_2_reg)*(
             # infected partner not on HCT/ART
             D3P1_m32_reg*(H1f2/Nf2) + D3G1_m32_reg*(Y1f2/Nf2) + D3X1_m32_reg*(Z1f2/Nf2) + D3K1_m32_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D3P2_m32_reg*(H2f2/Nf2) + D3G2_m32_reg*(Y2f2/Nf2) + D3X2_m32_reg*(Z2f2/Nf2) + D3K2_m32_reg*(A2f2/Nf2)
             # infected partner on ART
             + D3P3_m32_reg*(H3f2/Nf2) + D3G3_m32_reg*(Y3f2/Nf2) + D3X3_m32_reg*(Z3f2/Nf2) + D3K3_m32_reg*(A3f2/Nf2)
           )
           +(dm33_reg*rho_m3_3_reg)*(
             # infected partner not on HCT/ART
             D3P1_m33_reg*(H1f3/Nf3) + D3G1_m33_reg*(Y1f3/Nf3) + D3X1_m33_reg*(Z1f3/Nf3) + D3K1_m33_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D3P2_m33_reg*(H2f3/Nf3) + D3G2_m33_reg*(Y2f3/Nf3) + D3X2_m33_reg*(Z2f3/Nf3) + D3K2_m33_reg*(A2f3/Nf3)
             # infected partner on ART
             + D3P3_m33_reg*(H3f3/Nf3) + D3G3_m33_reg*(Y3f3/Nf3) + D3X3_m33_reg*(Z3f3/Nf3) + D3K3_m33_reg*(A3f3/Nf3)
           )
           # casual
           +(dm33_casual*rho_m3_3_casual)*(
             # infected partner not on HCT/ART
             D3P1_m33_casual*(H1f3/Nf3) + D3G1_m33_casual*(Y1f3/Nf3) + D3X1_m33_casual*(Z1f3/Nf3) + D3K1_m33_casual*(A1f3/Nf3)
             # infected partner on HCT
             + D3P2_m33_casual*(H2f3/Nf3) + D3G2_m33_casual*(Y2f3/Nf3) + D3X2_m33_casual*(Z2f3/Nf3) + D3K2_m33_casual*(A2f3/Nf3)
             # infected partner on ART
             + D3P3_m33_casual*(H3f3/Nf3) + D3G3_m33_casual*(Y3f3/Nf3) + D3X3_m33_casual*(Z3f3/Nf3) + D3K3_m33_casual*(A3f3/Nf3)
           )
         )
         
         
         ###
         # Circumcised males on PrEP
         pi4_m1 = (
           # low-low
           (dm11_main*rho_m1_1_main)*(
             # infected partner not on HCT/ART
             D4P1_m11_main*(H1f1/Nf1) + D4G1_m11_main*(Y1f1/Nf1) + D4X1_m11_main*(Z1f1/Nf1) + D4K1_m11_main*(A1f1/Nf1)
             # infected partner on HCT
             + D4P2_m11_main*(H2f1/Nf1) + D4G2_m11_main*(Y2f1/Nf1) + D4X2_m11_main*(Z2f1/Nf1) + D4K2_m11_main*(A2f1/Nf1)
             # infected partner on ART
             + D4P3_m11_main*(H3f1/Nf1) + D4G3_m11_main*(Y3f1/Nf1) + D4X3_m11_main*(Z3f1/Nf1) + D4K3_m11_main*(A3f1/Nf1)
           )
           # low-medium
           +(dm12_main*rho_m1_2_main)*(
             # infected partner not on HCT/ART
             D4P1_m12_main*(H1f2/Nf2) + D4G1_m12_main*(Y1f2/Nf2) + D4X1_m12_main*(Z1f2/Nf2) + D4K1_m12_main*(A1f2/Nf2)
             # infected partner on HCT
             + D4P2_m12_main*(H2f2/Nf2) + D4G2_m12_main*(Y2f2/Nf2) + D4X2_m12_main*(Z2f2/Nf2) + D4K2_m12_main*(A2f2/Nf2)
             # infected partner on ART
             + D4P3_m12_main*(H3f2/Nf2) + D4G3_m12_main*(Y3f2/Nf2) + D4X3_m12_main*(Z3f2/Nf2) + D4K3_m12_main*(A3f2/Nf2)
           )
           # low-high
           +(dm13_main*rho_m1_3_main)*(
             # infected partner not on HCT/ART
             D4P1_m13_main*(H1f3/Nf3) + D4G1_m13_main*(Y1f3/Nf3) + D4X1_m13_main*(Z1f3/Nf3) + D4K1_m13_main*(A1f3/Nf3)
             # infected partner on HCT
             + D4P2_m13_main*(H2f3/Nf3) + D4G2_m13_main*(Y2f3/Nf3) + D4X2_m13_main*(Z2f3/Nf3) + D4K2_m13_main*(A2f3/Nf3)
             # infected partner on ART
             + D4P3_m13_main*(H3f3/Nf3) + D4G3_m13_main*(Y3f3/Nf3) + D4X3_m13_main*(Z3f3/Nf3) + D4K3_m13_main*(A3f3/Nf3)
           )
         )
         
         pi4_m2 = (
           # main
           (dm21_main*rho_m2_1_main)*(
             # infected partner not on HCT/ART
             D4P1_m21_main*(H1f1/Nf1) + D4G1_m21_main*(Y1f1/Nf1) + D4X1_m21_main*(Z1f1/Nf1) + D4K1_m21_main*(A1f1/Nf1)
             # infected partner on HCT
             + D4P2_m21_main*(H2f1/Nf1) + D4G2_m21_main*(Y2f1/Nf1) + D4X2_m21_main*(Z2f1/Nf1) + D4K2_m21_main*(A2f1/Nf1)
             # infected partner on ART
             + D4P3_m21_main*(H3f1/Nf1) + D4G3_m21_main*(Y3f1/Nf1) + D4X3_m21_main*(Z3f1/Nf1) + D4K3_m21_main*(A3f1/Nf1)
           )
           +(dm22_main*rho_m2_2_main)*(
             # infected partner not on HCT/ART
             D4P1_m22_main*(H1f2/Nf2) + D4G1_m22_main*(Y1f2/Nf2) + D4X1_m22_main*(Z1f2/Nf2) + D4K1_m22_main*(A1f2/Nf2)
             # infected partner on HCT
             + D4P2_m22_main*(H2f2/Nf2) + D4G2_m22_main*(Y2f2/Nf2) + D4X2_m22_main*(Z2f2/Nf2) + D4K2_m22_main*(A2f2/Nf2)
             # infected partner on ART
             + D4P3_m22_main*(H3f2/Nf2) + D4G3_m22_main*(Y3f2/Nf2) + D4X3_m22_main*(Z3f2/Nf2) + D4K3_m22_main*(A3f2/Nf2)
           )
           +(dm23_main*rho_m2_3_main)*(
             # infected partner not on HCT/ART
             D4P1_m23_main*(H1f3/Nf3) + D4G1_m23_main*(Y1f3/Nf3) + D4X1_m23_main*(Z1f3/Nf3) + D4K1_m23_main*(A1f3/Nf3)
             # infected partner on HCT
             + D4P2_m23_main*(H2f3/Nf3) + D4G2_m23_main*(Y2f3/Nf3) + D4X2_m23_main*(Z2f3/Nf3) + D4K2_m23_main*(A2f3/Nf3)
             # infected partner on ART
             + D4P3_m23_main*(H3f3/Nf3) + D4G3_m23_main*(Y3f3/Nf3) + D4X3_m23_main*(Z3f3/Nf3) + D4K3_m23_main*(A3f3/Nf3)
           )
           # regular
           +(dm22_reg*rho_m2_2_reg)*(
             # infected partner not on HCT/ART
             D4P1_m22_reg*(H1f2/Nf2) + D4G1_m22_reg*(Y1f2/Nf2) + D4X1_m22_reg*(Z1f2/Nf2) + D4K1_m22_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D4P2_m22_reg*(H2f2/Nf2) + D4G2_m22_reg*(Y2f2/Nf2) + D4X2_m22_reg*(Z2f2/Nf2) + D4K2_m22_reg*(A2f2/Nf2)
             # infected partner on ART
             + D4P3_m22_reg*(H3f2/Nf2) + D4G3_m22_reg*(Y3f2/Nf2) + D4X3_m22_reg*(Z3f2/Nf2) + D4K3_m22_reg*(A3f2/Nf2)
           )
           +(dm23_reg*rho_m2_3_reg)*(
             # infected partner not on HCT/ART
             D4P1_m23_reg*(H1f3/Nf3) + D4G1_m23_reg*(Y1f3/Nf3) + D4X1_m23_reg*(Z1f3/Nf3) + D4K1_m23_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D4P2_m23_reg*(H2f3/Nf3) + D4G2_m23_reg*(Y2f3/Nf3) + D4X2_m23_reg*(Z2f3/Nf3) + D4K2_m23_reg*(A2f3/Nf3)
             # infected partner on ART
             + D4P3_m23_reg*(H3f3/Nf3) + D4G3_m23_reg*(Y3f3/Nf3) + D4X3_m23_reg*(Z3f3/Nf3) + D4K3_m23_reg*(A3f3/Nf3)
           )
         )
         
         pi4_m3 = (
           # main
           (dm31_main*rho_m3_1_main)*(
             # infected partner not on HCT/ART
             D4P1_m31_main*(H1f1/Nf1) + D4G1_m31_main*(Y1f1/Nf1) + D4X1_m31_main*(Z1f1/Nf1) + D4K1_m31_main*(A1f1/Nf1)
             # infected partner on HCT
             + D4P2_m31_main*(H2f1/Nf1) + D4G2_m31_main*(Y2f1/Nf1) + D4X2_m31_main*(Z2f1/Nf1) + D4K2_m31_main*(A2f1/Nf1)
             # infected partner on ART
             + D4P3_m31_main*(H3f1/Nf1) + D4G3_m31_main*(Y3f1/Nf1) + D4X3_m31_main*(Z3f1/Nf1) + D4K3_m31_main*(A3f1/Nf1)
           )
           +(dm32_main*rho_m3_2_main)*(
             # infected partner not on HCT/ART
             D4P1_m32_main*(H1f2/Nf2) + D4G1_m32_main*(Y1f2/Nf2) + D4X1_m32_main*(Z1f2/Nf2) + D4K1_m32_main*(A1f2/Nf2)
             # infected partner on HCT
             + D4P2_m32_main*(H2f2/Nf2) + D4G2_m32_main*(Y2f2/Nf2) + D4X2_m32_main*(Z2f2/Nf2) + D4K2_m32_main*(A2f2/Nf2)
             # infected partner on ART
             + D4P3_m32_main*(H3f2/Nf2) + D4G3_m32_main*(Y3f2/Nf2) + D4X3_m32_main*(Z3f2/Nf2) + D4K3_m32_main*(A3f2/Nf2)
           )
           +(dm33_main*rho_m3_3_main)*(
             # infected partner not on HCT/ART
             D4P1_m33_main*(H1f3/Nf3) + D4G1_m33_main*(Y1f3/Nf3) + D4X1_m33_main*(Z1f3/Nf3) + D4K1_m33_main*(A1f3/Nf3)
             # infected partner on HCT
             + D4P2_m33_main*(H2f3/Nf3) + D4G2_m33_main*(Y2f3/Nf3) + D4X2_m33_main*(Z2f3/Nf3) + D4K2_m33_main*(A2f3/Nf3)
             # infected partner on ART
             + D4P3_m33_main*(H3f3/Nf3) + D4G3_m33_main*(Y3f3/Nf3) + D4X3_m33_main*(Z3f3/Nf3) + D4K3_m33_main*(A3f3/Nf3)
           )
           # regular
           +(dm32_reg*rho_m3_2_reg)*(
             # infected partner not on HCT/ART
             D4P1_m32_reg*(H1f2/Nf2) + D4G1_m32_reg*(Y1f2/Nf2) + D4X1_m32_reg*(Z1f2/Nf2) + D4K1_m32_reg*(A1f2/Nf2)
             # infected partner on HCT
             + D4P2_m32_reg*(H2f2/Nf2) + D4G2_m32_reg*(Y2f2/Nf2) + D4X2_m32_reg*(Z2f2/Nf2) + D4K2_m32_reg*(A2f2/Nf2)
             # infected partner on ART
             + D4P3_m32_reg*(H3f2/Nf2) + D4G3_m32_reg*(Y3f2/Nf2) + D4X3_m32_reg*(Z3f2/Nf2) + D4K3_m32_reg*(A3f2/Nf2)
           )
           +(dm33_reg*rho_m3_3_reg)*(
             # infected partner not on HCT/ART
             D4P1_m33_reg*(H1f3/Nf3) + D4G1_m33_reg*(Y1f3/Nf3) + D4X1_m33_reg*(Z1f3/Nf3) + D4K1_m33_reg*(A1f3/Nf3)
             # infected partner on HCT
             + D4P2_m33_reg*(H2f3/Nf3) + D4G2_m33_reg*(Y2f3/Nf3) + D4X2_m33_reg*(Z2f3/Nf3) + D4K2_m33_reg*(A2f3/Nf3)
             # infected partner on ART
             + D4P3_m33_reg*(H3f3/Nf3) + D4G3_m33_reg*(Y3f3/Nf3) + D4X3_m33_reg*(Z3f3/Nf3) + D4K3_m33_reg*(A3f3/Nf3)
           )
           # casual
           +(dm33_casual*rho_m3_3_casual)*(
             # infected partner not on HCT/ART
             D4P1_m33_casual*(H1f3/Nf3) + D4G1_m33_casual*(Y1f3/Nf3) + D4X1_m33_casual*(Z1f3/Nf3) + D4K1_m33_casual*(A1f3/Nf3)
             # infected partner on HCT
             + D4P2_m33_casual*(H2f3/Nf3) + D4G2_m33_casual*(Y2f3/Nf3) + D4X2_m33_casual*(Z2f3/Nf3) + D4K2_m33_casual*(A2f3/Nf3)
             # infected partner on ART
             + D4P3_m33_casual*(H3f3/Nf3) + D4G3_m33_casual*(Y3f3/Nf3) + D4X3_m33_casual*(Z3f3/Nf3) + D4K3_m33_casual*(A3f3/Nf3)
           )
         )
         
         
         ## Incidence per 100 PY
         IRf1 = ((pi1_f1*S1f1 + pi2_f1*S2f1)/(S1f1+S2f1)*100)*12
         IRf2 = ((pi1_f2*S1f2 + pi2_f2*S2f2)/(S1f2+S2f2)*100)*12
         IRf3 = ((pi1_f3*S1f3 + pi2_f3*S2f3)/(S1f3+S2f3)*100)*12
         IRm1 = ((pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1)/(S1m1+S2m1+S3m1+S4m1)*100)*12
         IRm2 = ((pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2)/(S1m2+S2m2+S3m2+S4m2)*100)*12
         IRm3 = ((pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)/(S1m3+S2m3+S3m3+S4m3)*100)*12
         
         person.time = S1f1+S2f1+S1m1+S2m1+S3m1+S4m1
         +S1f2+S2f2+S1m2+S2m2+S3m2+S4m2
         +S1f3+S2f3+S1m3+S2m3+S3m3+S4m3
         
         IRtot = ((pi1_f1*S1f1 + pi2_f1*S2f1 + pi1_f2*S1f2 + pi2_f2*S2f2 + pi1_f3*S1f3 + pi2_f3*S2f3
                   + pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1 + pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2
                   + pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)/person.time*100)*12
         
         ## New infections
         NItot = (pi1_f1*S1f1 + pi2_f1*S2f1 + pi1_f2*S1f2 + pi2_f2*S2f2 + pi1_f3*S1f3 + pi2_f3*S2f3
                  + pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1 + pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2
                  + pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)
         
         NIf1 = (pi1_f1*S1f1 + pi2_f1*S2f1)
         NIf2 = (pi1_f2*S1f2 + pi2_f2*S2f2)
         NIf3 = (pi1_f3*S1f3 + pi2_f3*S2f3)
         NIm1 = (pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1)         
         NIm2 = (pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2)         
         NIm3 = (pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3)         

                  
         # AIDS-related mortality per 100 PY
         AMtot = ((delta1*A1f1 + delta1*A2f1 + delta2*A3f1 + delta1*A1f2 + delta1*A2f2 + delta2*A3f2 + delta1*A1f3 + delta1*A2f3 + delta2*A3f3
                   + delta1*A1m1 + delta1*A2m1 + delta2*A3m1 + delta1*A1m2 + delta1*A2m2 + delta2*A3m2 + delta1*A1m3 + delta1*A2m3 + delta2*A3m3)/N*100)*12
         
         AMf1 = ((delta1*A1f1 + delta1*A2f1 + delta2*A3f1)/Nf1*100)*12
         AMf2 = ((delta1*A1f2 + delta1*A2f2 + delta2*A3f2)/Nf2*100)*12
         AMf3 = ((delta1*A1f3 + delta1*A2f3 + delta2*A3f3)/Nf3*100)*12
         AMm1 = ((delta1*A1m1 + delta1*A2m1 + delta2*A3m1)/Nm1*100)*12
         AMm2 = ((delta1*A1m2 + delta1*A2m2 + delta2*A3m2)/Nm2*100)*12
         AMm3 = ((delta1*A1m3 + delta1*A2m3 + delta2*A3m3)/Nm3*100)*12
         
         
         ###########
         ### ODE ###
         ###########
         
         ### Females       
         dS1f1 = mu*Nf1 + om*S2f1 - (pi1_f1 + tauTS*theta + psi)*S1f1
         dS2f1 = tauTS*theta*S1f1 - (pi2_f1 + om + psi)*S2f1
         dH1f1 = pi1_f1*S1f1 + pi2_f1*S2f1 - (nu1 + tauT + psi)*H1f1
         dH2f1 = tauT*H1f1 + gamma*H3f1 - (nu1 + tauA1 + psi)*H2f1
         dH3f1 = tauA1*H2f1 - (nu2 + gamma + psi)*H3f1
         dY1f1 = nu1*H1f1 - (omega1 + tauT + psi)*Y1f1
         dY2f1 = nu1*H2f1 + tauT*Y1f1 + gamma*Y3f1 - (omega1 + tauA2 + psi)*Y2f1
         dY3f1 = nu2*H3f1 + tauA2*Y2f1 - (omega2 + gamma + psi)*Y3f1
         dZ1f1 = omega1*Y1f1 - (epi1 + tauT + psi)*Z1f1
         dZ2f1 = omega1*Y2f1 + tauT*Z1f1 + gamma*Z3f1 - (epi1 + tauA3 + psi)*Z2f1
         dZ3f1 = omega2*Y3f1 + tauA3*Z2f1 - (epi2 + gamma + psi)*Z3f1
         dA1f1 = epi1*Z1f1 - (delta1 + tauT + psi)*A1f1
         dA2f1 = epi1*Z2f1 + tauT*A1f1 + gamma*A3f1 - (delta1 + tauA4 + psi)*A2f1
         dA3f1 = epi2*Z3f1 + tauA4*A2f1 - (delta2 + gamma + psi)*A3f1
         
         dS1f2 = mu*Nf2 + om*S2f2 - (pi1_f2 + tauTS*theta + psi)*S1f2
         dS2f2 = tauTS*theta*S1f2 - (pi2_f2 + om + psi)*S2f2
         dH1f2 = pi1_f2*S1f2 + pi2_f2*S2f2 - (nu1 + tauT + psi)*H1f2
         dH2f2 = tauT*H1f2 + gamma*H3f2 - (nu1 + tauA1 + psi)*H2f2
         dH3f2 = tauA1*H2f2 - (nu2 + gamma + psi)*H3f2
         dY1f2 = nu1*H1f2 - (omega1 + tauT + psi)*Y1f2
         dY2f2 = nu1*H2f2 + tauT*Y1f2 + gamma*Y3f2 - (omega1 + tauA2 + psi)*Y2f2
         dY3f2 = nu2*H3f2 + tauA2*Y2f2 - (omega2 + gamma + psi)*Y3f2
         dZ1f2 = omega1*Y1f2 - (epi1 + tauT + psi)*Z1f2
         dZ2f2 = omega1*Y2f2 + tauT*Z1f2 + gamma*Z3f2 - (epi1 + tauA3 + psi)*Z2f2
         dZ3f2 = omega2*Y3f2 + tauA3*Z2f2 - (epi2 + gamma + psi)*Z3f2
         dA1f2 = epi1*Z1f2 - (delta1 + tauT + psi)*A1f2
         dA2f2 = epi1*Z2f2 + tauT*A1f2 + gamma*A3f2 - (delta1 + tauA4 + psi)*A2f2
         dA3f2 = epi2*Z3f2 + tauA4*A2f2 - (delta2 + gamma + psi)*A3f2
         
         dS1f3 = mu*Nf3 + om*S2f3 - (pi1_f3 + tauTS*theta + psi)*S1f3
         dS2f3 = tauTS*theta*S1f3 - (pi2_f3 + om + psi)*S2f3
         dH1f3 = pi1_f3*S1f3 + pi2_f3*S2f3 - (nu1 + tauT + psi)*H1f3
         dH2f3 = tauT*H1f3 + gamma*H3f3 - (nu1 + tauA1 + psi)*H2f3
         dH3f3 = tauA1*H2f3 - (nu2 + gamma + psi)*H3f3
         dY1f3 = nu1*H1f3 - (omega1 + tauT + psi)*Y1f3
         dY2f3 = nu1*H2f3 + tauT*Y1f3 + gamma*Y3f3 - (omega1 + tauA2 + psi)*Y2f3
         dY3f3 = nu2*H3f3 + tauA2*Y2f3 - (omega2 + gamma + psi)*Y3f3
         dZ1f3 = omega1*Y1f3 - (epi1 + tauT + psi)*Z1f3
         dZ2f3 = omega1*Y2f3 + tauT*Z1f3 + gamma*Z3f3 - (epi1 + tauA3 + psi)*Z2f3
         dZ3f3 = omega2*Y3f3 + tauA3*Z2f3 - (epi2 + gamma + psi)*Z3f3
         dA1f3 = epi1*Z1f3 - (delta1 + tauT + psi)*A1f3
         dA2f3 = epi1*Z2f3 + tauT*A1f3 + gamma*A3f3 - (delta1 + tauA4 + psi)*A2f3
         dA3f3 = epi2*Z3f3 + tauA4*A2f3 - (delta2 + gamma + psi)*A3f3
         
         ### Males
         
         dS1m1 = (1-prop_circ)*mu*Nm1 + om*S2m1 - (pi1_m1 + tauC + HCT_m*tauTS*theta + psi)*S1m1
         dS2m1 = HCT_m*tauTS*theta*S1m1 - (pi2_m1 + tauC + om + psi)*S2m1
         dS3m1 = prop_circ*mu*Nm1 + tauC*S1m1 + om*S4m1 - (pi3_m1 + HCT_m*tauTS*theta + psi)*S3m1
         dS4m1 = HCT_m*tauTS*theta*S3m1 + tauC*S2m1 - (pi4_m1 + om + psi)*S4m1
         dH1m1 = pi1_m1*S1m1 + pi2_m1*S2m1 + pi3_m1*S3m1 + pi4_m1*S4m1 - (nu1 + HCT_m*tauT + psi)*H1m1
         dH2m1 = HCT_m*tauT*H1m1 + gamma*H3m1 - (nu1 + ART_m*tauA1 + psi)*H2m1
         dH3m1 = ART_m*tauA1*H2m1 - (nu2 + gamma + psi)*H3m1
         dY1m1 = nu1*H1m1 - (omega1 + HCT_m*tauT + psi)*Y1m1
         dY2m1 = nu1*H2m1 + HCT_m*tauT*Y1m1 + gamma*Y3m1 - (omega1 + ART_m*tauA2 + psi)*Y2m1
         dY3m1 = nu2*H3m1 + ART_m*tauA2*Y2m1 - (omega2 + gamma + psi)*Y3m1
         dZ1m1 = omega1*Y1m1 - (epi1 + HCT_m*tauT + psi)*Z1m1
         dZ2m1 = omega1*Y2m1 + HCT_m*tauT*Z1m1 + gamma*Z3m1 - (epi1 + ART_m*tauA3 + psi)*Z2m1
         dZ3m1 = omega2*Y3m1 + ART_m*tauA3*Z2m1 - (epi2 + gamma + psi)*Z3m1
         dA1m1 = epi1*Z1m1 - (delta1 + HCT_m*tauT + psi)*A1m1
         dA2m1 = epi1*Z2m1 + HCT_m*tauT*A1m1 + gamma*A3m1 - (delta1 + ART_m*tauA4 + psi)*A2m1
         dA3m1 = epi2*Z3m1 + ART_m*tauA4*A2m1 - (delta2 + gamma + psi)*A3m1
         
         
         dS1m2 = (1-prop_circ)*mu*Nm2 + om*S2m2 - (pi1_m2 + tauC + HCT_m*tauTS*theta + psi)*S1m2
         dS2m2 = HCT_m*tauTS*theta*S1m2 - (pi2_m2 + tauC + om + psi)*S2m2
         dS3m2 = prop_circ*mu*Nm2 + tauC*S1m2 + om*S4m2 - (pi3_m2 + HCT_m*tauTS*theta + psi)*S3m2
         dS4m2 = HCT_m*tauTS*theta*S3m2 + tauC*S2m2 - (pi4_m2 + om + psi)*S4m2
         dH1m2 = pi1_m2*S1m2 + pi2_m2*S2m2 + pi3_m2*S3m2 + pi4_m2*S4m2 - (nu1 + HCT_m*tauT + psi)*H1m2
         dH2m2 = HCT_m*tauT*H1m2 + gamma*H3m2 - (nu1 + ART_m*tauA1 + psi)*H2m2
         dH3m2 = ART_m*tauA1*H2m2 - (nu2 + gamma + psi)*H3m2
         dY1m2 = nu1*H1m2 - (omega1 + HCT_m*tauT + psi)*Y1m2
         dY2m2 = nu1*H2m2 + HCT_m*tauT*Y1m2 + gamma*Y3m2 - (omega1 + ART_m*tauA2 + psi)*Y2m2
         dY3m2 = nu2*H3m2 + ART_m*tauA2*Y2m2 - (omega2 + gamma + psi)*Y3m2
         dZ1m2 = omega1*Y1m2 - (epi1 + HCT_m*tauT + psi)*Z1m2
         dZ2m2 = omega1*Y2m2 + HCT_m*tauT*Z1m2 + gamma*Z3m2 - (epi1 + ART_m*tauA3 + psi)*Z2m2
         dZ3m2 = omega2*Y3m2 + ART_m*tauA3*Z2m2 - (epi2 + gamma + psi)*Z3m2
         dA1m2 = epi1*Z1m2 - (delta1 + HCT_m*tauT + psi)*A1m2
         dA2m2 = epi1*Z2m2 + HCT_m*tauT*A1m2 + gamma*A3m2 - (delta1 + ART_m*tauA4 + psi)*A2m2
         dA3m2 = epi2*Z3m2 + ART_m*tauA4*A2m2 - (delta2 + gamma + psi)*A3m2
         
         dS1m3 = (1-prop_circ)*mu*Nm3 + om*S2m3 - (pi1_m3 + tauC + HCT_m*tauTS*theta + psi)*S1m3
         dS2m3 = HCT_m*tauTS*theta*S1m3 - (pi2_m3 + tauC + om + psi)*S2m3
         dS3m3 = prop_circ*mu*Nm3 + tauC*S1m3 + om*S4m3 - (pi3_m3 + HCT_m*tauTS*theta + psi)*S3m3
         dS4m3 = HCT_m*tauTS*theta*S3m3 + tauC*S2m3 - (pi4_m3 + om + psi)*S4m3
         dH1m3 = pi1_m3*S1m3 + pi2_m3*S2m3 + pi3_m3*S3m3 + pi4_m3*S4m3 - (nu1 + HCT_m*tauT + psi)*H1m3
         dH2m3 = HCT_m*tauT*H1m3 + gamma*H3m3 - (nu1 + ART_m*tauA1 + psi)*H2m3
         dH3m3 = ART_m*tauA1*H2m3 - (nu2 + gamma + psi)*H3m3
         dY1m3 = nu1*H1m3 - (omega1 + HCT_m*tauT + psi)*Y1m3
         dY2m3 = nu1*H2m3 + HCT_m*tauT*Y1m3 + gamma*Y3m3 - (omega1 + ART_m*tauA2 + psi)*Y2m3
         dY3m3 = nu2*H3m3 + ART_m*tauA2*Y2m3 - (omega2 + gamma + psi)*Y3m3
         dZ1m3 = omega1*Y1m3 - (epi1 + HCT_m*tauT + psi)*Z1m3
         dZ2m3 = omega1*Y2m3 + HCT_m*tauT*Z1m3 + gamma*Z3m3 - (epi1 + ART_m*tauA3 + psi)*Z2m3
         dZ3m3 = omega2*Y3m3 + ART_m*tauA3*Z2m3 - (epi2 + gamma + psi)*Z3m3
         dA1m3 = epi1*Z1m3 - (delta1 + HCT_m*tauT + psi)*A1m3
         dA2m3 = epi1*Z2m3 + HCT_m*tauT*A1m3 + gamma*A3m3 - (delta1 + ART_m*tauA4 + psi)*A2m3
         dA3m3 = epi2*Z3m3 + ART_m*tauA4*A2m3 - (delta2 + gamma + psi)*A3m3
         
         
         ## Prevalence
         PRf1 = (H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1)/Nf1
         PRf2 = (H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2)/Nf2
         PRf3 = (H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3)/Nf3
         PRm1 = (H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1)/Nm1
         PRm2 = (H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2)/Nm2
         PRm3 = (H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/Nm3
         
         PRtot = (( H1f1+H2f1+H3f1+Y1f1+Y2f1+Y3f1+Z1f1+Z2f1+Z3f1+A1f1+A2f1+A3f1
                    +H1f2+H2f2+H3f2+Y1f2+Y2f2+Y3f2+Z1f2+Z2f2+Z3f2+A1f2+A2f2+A3f2
                    +H1f3+H2f3+H3f3+Y1f3+Y2f3+Y3f3+Z1f3+Z2f3+Z3f3+A1f3+A2f3+A3f3
                    +H1m1+H2m1+H3m1+Y1m1+Y2m1+Y3m1+Z1m1+Z2m1+Z3m1+A1m1+A2m1+A3m1
                    +H1m2+H2m2+H3m2+Y1m2+Y2m2+Y3m2+Z1m2+Z2m2+Z3m2+A1m2+A2m2+A3m2
                    +H1m3+H2m3+H3m3+Y1m3+Y2m3+Y3m3+Z1m3+Z2m3+Z3m3+A1m3+A2m3+A3m3)/N)
         
         propART = (H3f1+Y3f1+Z3f1+A3f1+H3f2+Y3f2+Z3f2+A3f2+H3f3+Y3f3+Z3f3+A3f3
                    +H3m1+Y3m1+Z3m1+A3m1+H3m2+Y3m2+Z3m2+A3m2+H3m3+Y3m3+Z3m3+A3m3)/(H3f1+Y3f1+Z3f1+A3f1+H3f2+Y3f2+Z3f2+A3f2+H3f3+Y3f3+Z3f3+A3f3
                                                                                   +H3m1+Y3m1+Z3m1+A3m1+H3m2+Y3m2+Z3m2+A3m2+H3m3+Y3m3+Z3m3+A3m3
                                                                                   +H2f1+Y2f1+Z2f1+A2f1+H2f2+Y2f2+Z2f2+A2f2+H2f3+Y2f3+Z2f3+A2f3
                                                                                   +H2m1+Y2m1+Z2m1+A2m1+H2m2+Y2m2+Z2m2+A2m2+H2m3+Y2m3+Z2m3+A2m3
                                                                                   +H1f1+Y1f1+Z1f1+A1f1+H1f2+Y1f2+Z1f2+A1f2+H1f3+Y1f3+Z1f3+A1f3
                                                                                   +H1m1+Y1m1+Z1m1+A1m1+H1m2+Y1m2+Z1m2+A1m2+H1m3+Y1m3+Z1m3+A1m3)
         ART.male = (H3m1+Y3m1+Z3m1+A3m1+H3m2+Y3m2+Z3m2+A3m2+H3m3+Y3m3+Z3m3+A3m3)/(Nm1+Nm2+Nm3-S.male)
         ART.female = (H3f1+Y3f1+Z3f1+A3f1+H3f2+Y3f2+Z3f2+A3f2+H3f3+Y3f3+Z3f3+A3f3)/(Nf1+Nf2+Nf3-S.female)
         
         # Proportion aware of HIV+ status
         prop.aware = (H2f1+H3f1+H2f2+H3f2+H2f3+H3f3+H2m1+H3m1+H2m2+H3m2+H2m3+H3m3+
                         Y2f1+Y3f1+Y2f2+Y3f2+Y2f3+Y3f3+Y2m1+Y3m1+Y2m2+Y3m2+Y2m3+Y3m3+
                         Z2f1+Z3f1+Z2f2+Z3f2+Z2f3+Z3f3+Z2m1+Z3m1+Z2m2+Z3m2+Z2m3+Z3m3+
                         A2f1+A3f1+A2f2+A3f2+A2f3+A3f3+A2m1+A3m1+A2m2+A3m2+A2m3+A3m3)/(Htot+Ytot+Ztot+Atot)
         aware.male = (H2m1+H3m1+H2m2+H3m2+H2m3+H3m3+Y2m1+Y3m1+Y2m2+Y3m2+Y2m3+Y3m3+Z2m1+Z3m1+Z2m2+Z3m2+Z2m3+Z3m3+A2m1+A3m1+A2m2+A3m2+A2m3+A3m3)/(Nm1+Nm2+Nm3-S.male)
         aware.female = (H2f1+H3f1+H2f2+H3f2+H2f3+H3f3+Y2f1+Y3f1+Y2f2+Y3f2+Y2f3+Y3f3+Z2f1+Z3f1+Z2f2+Z3f2+Z2f3+Z3f3+A2f1+A3f1+A2f2+A3f2+A2f3+A3f3)/(Nf1+Nf2+Nf3-S.female)
         
         # Proportion of susceptibles on PrEP
         propPREP = (S2f1+S2f2+S2f3+S2m1+S2m2+S2m3+S4m1+S4m2+S4m3)/Stot
         
         # Proportion of those aware on ART
         ART.aware = (H3f1+Y3f1+Z3f1+A3f1+H3f2+Y3f2+Z3f2+A3f2+H3f3+Y3f3+Z3f3+A3f3
                      +H3m1+Y3m1+Z3m1+A3m1+H3m2+Y3m2+Z3m2+A3m2+H3m3+Y3m3+Z3m3+A3m3)/(
                        H2f1+H3f1+H2f2+H3f2+H2f3+H3f3+H2m1+H3m1+H2m2+H3m2+H2m3+H3m3+
                          Y2f1+Y3f1+Y2f2+Y3f2+Y2f3+Y3f3+Y2m1+Y3m1+Y2m2+Y3m2+Y2m3+Y3m3+
                          Z2f1+Z3f1+Z2f2+Z3f2+Z2f3+Z3f3+Z2m1+Z3m1+Z2m2+Z3m2+Z2m3+Z3m3+
                          A2f1+A3f1+A2f2+A3f2+A2f3+A3f3+A2m1+A3m1+A2m2+A3m2+A2m3+A3m3)
         ART.aware.male = (H3m1+Y3m1+Z3m1+A3m1+H3m2+Y3m2+Z3m2+A3m2+H3m3+Y3m3+Z3m3+A3m3)/(H2m1+H3m1+H2m2+H3m2+H2m3+H3m3+Y2m1+Y3m1+Y2m2+Y3m2+Y2m3+Y3m3+Z2m1+Z3m1+Z2m2+Z3m2+Z2m3+Z3m3+A2m1+A3m1+A2m2+A3m2+A2m3+A3m3)
         ART.aware.female = (H3f1+Y3f1+Z3f1+A3f1+H3f2+Y3f2+Z3f2+A3f2+H3f3+Y3f3+Z3f3+A3f3)/(H2f1+H3f1+H2f2+H3f2+H2f3+H3f3+Y2f1+Y3f1+Y2f2+Y3f2+Y2f3+Y3f3+Z2f1+Z3f1+Z2f2+Z3f2+Z2f3+Z3f3+A2f1+A3f1+A2f2+A3f2+A2f3+A3f3)
         
         # Proportion of men circumcised         
         propCIRC = (S3m1+S3m2+S3m3+S4m1+S4m2+S4m3)/(Nm1+Nm2+Nm3)
         
         HCTrate = tauT
         PREPrateM = HCT_m*tauTS*theta
         PREPrateF = tauTS*theta
         
         ## output
         list(c(dS1f1,dH1f1,dY1f1,dZ1f1,dA1f1,
                dS1f2,dH1f2,dY1f2,dZ1f2,dA1f2,
                dS1f3,dH1f3,dY1f3,dZ1f3,dA1f3,
                dS1m1,dH1m1,dY1m1,dZ1m1,dA1m1,
                dS1m2,dH1m2,dY1m2,dZ1m2,dA1m2,
                dS1m3,dH1m3,dY1m3,dZ1m3,dA1m3,
                dS2f1,dH2f1,dY2f1,dZ2f1,dA2f1,
                dS2f2,dH2f2,dY2f2,dZ2f2,dA2f2,
                dS2f3,dH2f3,dY2f3,dZ2f3,dA2f3,
                dS2m1,dH2m1,dY2m1,dZ2m1,dA2m1,
                dS2m2,dH2m2,dY2m2,dZ2m2,dA2m2,
                dS2m3,dH2m3,dY2m3,dZ2m3,dA2m3,
                dS3m1,dS3m2,dS3m3,dS4m1,dS4m2,dS4m3,
                dH3f1,dY3f1,dZ3f1,dA3f1,
                dH3f2,dY3f2,dZ3f2,dA3f2,
                dH3f3,dY3f3,dZ3f3,dA3f3,
                dH3m1,dY3m1,dZ3m1,dA3m1,
                dH3m2,dY3m2,dZ3m2,dA3m2,
                dH3m3,dY3m3,dZ3m3,dA3m3),
              IRf1=IRf1,IRf2=IRf2,IRf3=IRf3,IRm1=IRm1,IRm2=IRm2,IRm3=IRm3,
              PRf1=PRf1,PRf2=PRf2,PRf3=PRf3,PRm1=PRm1,PRm2=PRm2,PRm3=PRm3,
              IRtot=IRtot,PRtot=PRtot,AMtot=AMtot,prop.aware=propaware,
              Nf1=Nf1, Nf2=Nf2, Nf3=Nf3, Nm1=Nm1, Nm2=Nm2, Nm3=Nm3, N=N, propART=propART, propCIRC=propCIRC,ARTinc=ARTinc,
              pi1_f1=pi1_f1, pi1_m1=pi1_m1, pi1_f2=pi1_f2, pi1_m2=pi1_m2, pi1_f3=pi1_f3, pi1_m3=pi1_m3,
              rho_f1_1_main=rho_f1_1_main,
              AMf1=AMf1, AMf2=AMf2, AMf3=AMf3, AMm1=AMm1, AMm2=AMm2, AMm3=AMm3,
              Stot=Stot,Htot=Htot,Ytot=Ytot,Ztot=Ztot,Atot=Atot,ART.aware=ART.aware,ART.drop=ART.drop,
              cond.main=cond_main,cond.reg=cond_reg,cond.casual=cond_casual,
              HCTrate = HCTrate, tauC = tauC,
              aware.male = aware.male, aware.female = aware.female,
              ART.male = ART.male, ART.female = ART.female, ART.aware.male = ART.aware.male, ART.aware.female = ART.aware.female,
              propPREP = propPREP, PREPrateM = PREPrateM, PREPrateF = PREPrateF, PrEP.drop = PrEP.drop, propPREP.male=propPREP.male, propPREP.female=propPREP.female,
              NItot=NItot, NIf1=NIf1, NIf2=NIf2, NIf3=NIf3, NIm1=NIm1, NIm2=NIm2, NIm3=NIm3)
       }
  )
}

########################
### SOLVE ODE SYSTEM ###
########################

state=c(S1f1=S1f1,H1f1=H1f1,Y1f1=Y1f1,Z1f1=Z1f1,A1f1=A1f1,
        S1f2=S1f2,H1f2=H1f2,Y1f2=Y1f2,Z1f2=Z1f2,A1f2=A1f2,
        S1f3=S1f3,H1f3=H1f3,Y1f3=Y1f3,Z1f3=Z1f3,A1f3=A1f3,
        S1m1=S1m1,H1m1=H1m1,Y1m1=Y1m1,Z1m1=Z1m1,A1m1=A1m1,
        S1m2=S1m2,H1m2=H1m2,Y1m2=Y1m2,Z1m2=Z1m2,A1m2=A1m2,
        S1m3=S1m3,H1m3=H1m3,Y1m3=Y1m3,Z1m3=Z1m3,A1m3=A1m3,
        S2f1=S2f1,H2f1=H2f1,Y2f1=Y2f1,Z2f1=Z2f1,A2f1=A2f1,
        S2f2=S2f2,H2f2=H2f2,Y2f2=Y2f2,Z2f2=Z2f2,A2f2=A2f2,
        S2f3=S2f3,H2f3=H2f3,Y2f3=Y2f3,Z2f3=Z2f3,A2f3=A2f3,
        S2m1=S2m1,H2m1=H2m1,Y2m1=Y2m1,Z2m1=Z2m1,A2m1=A2m1,
        S2m2=S2m2,H2m2=H2m2,Y2m2=Y2m2,Z2m2=Z2m2,A2m2=A2m2,
        S2m3=S2m3,H2m3=H2m3,Y2m3=Y2m3,Z2m3=Z2m3,A2m3=A2m3,
        S3m1=S3m1,S3m2=S3m2,S3m3=S3m3,S4m1=S4m1,S4m2=S4m2,S4m3=S4m3,
        H3f1=H3f1,Y3f1=Y3f1,Z3f1=Z3f1,A3f1=A3f1,
        H3f2=H3f2,Y3f2=Y3f2,Z3f2=Z3f2,A3f2=A3f2,
        H3f3=H3f3,Y3f3=Y3f3,Z3f3=Z3f3,A3f3=A3f3,
        H3m1=H3m1,Y3m1=Y3m1,Z3m1=Z3m1,A3m1=A3m1,
        H3m2=H3m2,Y3m2=Y3m2,Z3m2=Z3m2,A3m2=A3m2,
        H3m3=H3m3,Y3m3=Y3m3,Z3m3=Z3m3,A3m3=A3m3)

times=seq(0,720,by=1) #run for 45 years (= 1 generation of sexually active persons)

parameters=c(mu=mu, psi=psi,
             nu1=1/avg.nu1, omega1=1/avg.omega1, epi1=1/avg.epi1, delta1=1/avg.delta1,
             nu2=1/avg.nu2, omega2=1/avg.omega2, epi2=1/avg.epi2, delta2=1/avg.delta2)

require(deSolve) 

trial1_FM=as.data.frame(ode(y=state,times=times,func=SHYZA_FM,parms=parameters))
tail(trial1_FM)


#####################################################################################
## Output to assess intervention impact

# Relative reduction in prevalence
RP_f1 = (trial1_FM_BASE[721,"PRf1"] - trial1_FM[721,"PRf1"]) / trial1_FM_BASE[721,"PRf1"]
RP_f2 = (trial1_FM_BASE[721,"PRf2"] - trial1_FM[721,"PRf2"]) / trial1_FM_BASE[721,"PRf2"]
RP_f3 = (trial1_FM_BASE[721,"PRf3"] - trial1_FM[721,"PRf3"]) / trial1_FM_BASE[721,"PRf3"]
RP_m1 = (trial1_FM_BASE[721,"PRm1"] - trial1_FM[721,"PRm1"]) / trial1_FM_BASE[721,"PRm1"]
RP_m2 = (trial1_FM_BASE[721,"PRm2"] - trial1_FM[721,"PRm2"]) / trial1_FM_BASE[721,"PRm2"]
RP_m3 = (trial1_FM_BASE[721,"PRm3"] - trial1_FM[721,"PRm3"]) / trial1_FM_BASE[721,"PRm3"]
RP_tot = (trial1_FM_BASE[721,"PRtot"] - trial1_FM[721,"PRtot"]) / trial1_FM_BASE[721,"PRtot"]

# Relative reduction in incidence
RI_f1 = (trial1_FM_BASE[721,"IRf1"] - trial1_FM[721,"IRf1"]) / trial1_FM_BASE[721,"IRf1"]
RI_f2 = (trial1_FM_BASE[721,"IRf2"] - trial1_FM[721,"IRf2"]) / trial1_FM_BASE[721,"IRf2"]
RI_f3 = (trial1_FM_BASE[721,"IRf3"] - trial1_FM[721,"IRf3"]) / trial1_FM_BASE[721,"IRf3"]
RI_m1 = (trial1_FM_BASE[721,"IRm1"] - trial1_FM[721,"IRm1"]) / trial1_FM_BASE[721,"IRm1"]
RI_m2 = (trial1_FM_BASE[721,"IRm2"] - trial1_FM[721,"IRm2"]) / trial1_FM_BASE[721,"IRm2"]
RI_m3 = (trial1_FM_BASE[721,"IRm3"] - trial1_FM[721,"IRm3"]) / trial1_FM_BASE[721,"IRm3"]
RI_tot = (trial1_FM_BASE[721,"IRtot"] - trial1_FM[721,"IRtot"]) / trial1_FM_BASE[721,"IRtot"]

# Relative reduction in AIDS-related mortality
RA_f1 = (trial1_FM_BASE[721,"AMf1"] - trial1_FM[721,"AMf1"]) / trial1_FM_BASE[721,"AMf1"]
RA_f2 = (trial1_FM_BASE[721,"AMf2"] - trial1_FM[721,"AMf2"]) / trial1_FM_BASE[721,"AMf2"]
RA_f3 = (trial1_FM_BASE[721,"AMf3"] - trial1_FM[721,"AMf3"]) / trial1_FM_BASE[721,"AMf3"]
RA_m1 = (trial1_FM_BASE[721,"AMm1"] - trial1_FM[721,"AMm1"]) / trial1_FM_BASE[721,"AMm1"]
RA_m2 = (trial1_FM_BASE[721,"AMm2"] - trial1_FM[721,"AMm2"]) / trial1_FM_BASE[721,"AMm2"]
RA_m3 = (trial1_FM_BASE[721,"AMm3"] - trial1_FM[721,"AMm3"]) / trial1_FM_BASE[721,"AMm3"]
RA_tot = (trial1_FM_BASE[721,"AMtot"] - trial1_FM[721,"AMtot"]) / trial1_FM_BASE[721,"AMtot"]

# Reduction in total number of infections
RNI_tot = (sum(trial1_FM_BASE[721, c("Htot","Ytot","Ztot","Atot")]) - sum(trial1_FM[721, c("Htot","Ytot","Ztot","Atot")])) / sum(trial1_FM_BASE[721, c("Htot","Ytot","Ztot","Atot")])

# Reduction in number of NEW infections
RNInew_tot = (sum(trial1_FM_BASE[541:721, "NItot"]) - sum(trial1_FM[541:721, "NItot"])) / sum(trial1_FM_BASE[541:721, "NItot"])
RNInew_f1 = (sum(trial1_FM_BASE[541:721, "NIf1"]) - sum(trial1_FM[541:721, "NIf1"])) / sum(trial1_FM_BASE[541:721, "NIf1"])
RNInew_f2 = (sum(trial1_FM_BASE[541:721, "NIf2"]) - sum(trial1_FM[541:721, "NIf2"])) / sum(trial1_FM_BASE[541:721, "NIf2"])
RNInew_f3 = (sum(trial1_FM_BASE[541:721, "NIf3"]) - sum(trial1_FM[541:721, "NIf3"])) / sum(trial1_FM_BASE[541:721, "NIf3"])
RNInew_m1 = (sum(trial1_FM_BASE[541:721, "NIm1"]) - sum(trial1_FM[541:721, "NIm1"])) / sum(trial1_FM_BASE[541:721, "NIm1"])
RNInew_m2 = (sum(trial1_FM_BASE[541:721, "NIm2"]) - sum(trial1_FM[541:721, "NIm2"])) / sum(trial1_FM_BASE[541:721, "NIm2"])
RNInew_m3 = (sum(trial1_FM_BASE[541:721, "NIm3"]) - sum(trial1_FM[541:721, "NIm3"])) / sum(trial1_FM_BASE[541:721, "NIm3"])

## SAVE OUTPUT
out = data.frame(sim=sim, RNI_tot=RNI_tot, RNInew_tot=RNInew_tot,
                 RP_f1=RP_f1, RP_f2=RP_f2, RP_f3=RP_f3, RP_m1=RP_m1, RP_m2=RP_m2, RP_m3=RP_m3, RP_tot=RP_tot,
                 RI_f1=RI_f1, RI_f2=RI_f2, RI_f3=RI_f3, RI_m1=RI_m1, RI_m2=RI_m2, RI_m3=RI_m3, RI_tot=RI_tot,
                 RA_f1=RA_f1, RA_f2=RA_f2, RA_f3=RA_f3, RA_m1=RA_m1, RA_m2=RA_m2, RA_m3=RA_m3, RA_tot=RA_tot,
                 RNInew_f1=RNInew_f1, RNInew_f2=RNInew_f2, RNInew_f3=RNInew_f3, RNInew_m1=RNInew_m1, RNInew_m2=RNInew_m2, RNInew_m3=RNInew_m3)

write.table(out,file=paste0("uncertainty_sim_", sim, ".csv"),append=F,col.names=FALSE,row.names=FALSE,sep=",")


