args <- commandArgs(trailingOnly =  TRUE)
time <- as.numeric(args[1])
lambda <- as.numeric(args[2])
psi <- as.numeric(args[3])
C_f <- as.numeric(args[4])
C_m <- as.numeric(args[5])
high_v1 <- as.numeric(args[6])
high_v2 <- as.numeric(args[7])
beta_f <- as.numeric(args[8])
beta_m <- as.numeric(args[9])
cond_m1 <- as.numeric(args[10])
cond_m2 <- as.numeric(args[11])
cond_m3 <- as.numeric(args[12])
cond_f2 <- as.numeric(args[13])
cond_f3 <- as.numeric(args[14])
cond_eff <- as.numeric(args[15])
acts_m3_casual <- as.numeric(args[16])
acts_f3_casual <- as.numeric(args[17])
acts_m2_reg <- as.numeric(args[18])
acts_m3_reg <- as.numeric(args[19])
acts_f2_reg <- as.numeric(args[20])
acts_f3_reg <- as.numeric(args[21])
acts_m1_main <- as.numeric(args[22])
acts_m2_main <- as.numeric(args[23])
acts_m3_main <- as.numeric(args[24])
acts_f1_main <- as.numeric(args[25])
acts_f2_main <- as.numeric(args[26])
acts_f3_main <- as.numeric(args[27])
nu <- as.numeric(args[28])
omega <- as.numeric(args[29])
epi <- as.numeric(args[30])
delta <- as.numeric(args[31])
sim <- as.integer(args[32])
dm1 <- as.numeric(args[33])
dm2 <- as.numeric(args[34])
dm3 <- as.numeric(args[35])
df1 <- as.numeric(args[36])
df2 <- as.numeric(args[37])
df3 <- as.numeric(args[38])
N <- as.integer(args[39])
A_main <- as.numeric(args[40])
A_regular <- as.numeric(args[41])
A_casual <- as.numeric(args[42])
int_cond <- as.numeric(args[43])
circ <- as.numeric(args[44])
prep_m <- as.numeric(args[45])
prep_f <- as.numeric(args[46])
ART_m <- as.numeric(args[47])
ART_f <- as.numeric(args[48])
prop_hct_m <- as.numeric(args[49])
prop_hct_f <- as.numeric(args[50])


library(deSolve)

### ODES

SHYZA_FM=function(t, state, parameters)
{

  with(as.list(c(state,parameters)),
       {

         ## ODE
         ### Females       
         dSf1 = lambda*Nf1 - (pi_f1 + psi)*Sf1
         dHf1 = pi_f1*Sf1 - (nu+psi)*Hf1
         dYf1 = nu*Hf1 - (omega+psi)*Yf1
         dZf1 = omega*Yf1-(epi_f1+psi)*Zf1
         dAf1 = epi_f1*Zf1-(delta+psi)*Af1
         
         dSf2 = lambda*Nf2 - (pi_f2 + psi)*Sf2
         dHf2 = pi_f2*Sf2 - (nu+psi)*Hf2
         dYf2 = nu*Hf2 - (omega+psi)*Yf2
         dZf2 = omega*Yf2-(epi_f2+psi)*Zf2
         dAf2 = epi_f2*Zf2-(delta+psi)*Af2
         
         dSf3 = lambda*Nf3 - (pi_f3 + psi)*Sf3
         dHf3 = pi_f3*Sf3 - (nu+psi)*Hf3
         dYf3 = nu*Hf3 - (omega+psi)*Yf3
         dZf3 = omega*Yf3-(epi_f3+psi)*Zf3
         dAf3 = epi_f3*Zf3-(delta+psi)*Af3
         
         ### Males
         
         dSm1 = lambda*Nm1 - (pi_m1 + psi)*Sm1
         dHm1 = pi_m1*Sm1 - (nu+psi)*Hm1
         dYm1 = nu*Hm1 - (omega+psi)*Ym1
         dZm1 = omega*Ym1-(epi_m1+psi)*Zm1
         dAm1 = epi_m1*Zm1-(delta+psi)*Am1
         
         dSm2 = lambda*Nm2 - (pi_m2 + psi)*Sm2
         dHm2 = pi_m2*Sm2 - (nu+psi)*Hm2
         dYm2 = nu*Hm2 - (omega+psi)*Ym2
         dZm2 = omega*Ym2-(epi_m2+psi)*Zm2
         dAm2 = epi_m2*Zm2-(delta+psi)*Am2
         
         dSm3 = lambda*Nm3 - (pi_m3 + psi)*Sm3
         dHm3 = pi_m2*Sm3 - (nu+psi)*Hm3
         dYm3 = nu*Hm3 - (omega+psi)*Ym3
         dZm3 = omega*Ym3-(epi_m3+psi)*Zm3
         dAm3 = epi_m3*Zm3-(delta+psi)*Am3
         
         list(c(dSf1,dHf1,dYf1,dZf1,dAf1,dSf2,dHf2,dYf2,dZf2,dAf2,dSf3,dHf3,dYf3,dZf3,dAf3,
                dSm1,dHm1,dYm1,dZm1,dAm1,dSm2,dHm2,dYm2,dZm2,dAm2,dSm3,dHm3,dYm3,dZm3,dAm3))
       }
  )
}

dm1=dm1; dm2=dm2; df1=df1; df2=df2
dm3=dm3 ; df3=df3

A_main=A_main  
A_regular=A_regular
A_casual=A_casual

Nf1=816.02
Nf2=29.62448
Nf3=68.16089
Nm1=788.9862
Nm2=145.3718
Nm3=120.6077


Hm1=2.363762
Hm2=0.642904
Hf1=4.356965
Hf2=0.2421061
Hf3=0.6300652
Hm3=0.5334462

Ym1=60.79127
Ym2=16.6315
Yf1=112.8477
Yf2=6.257212
Yf3=16.21821
Ym3=13.77717

Zm1=14.67272
Zm2=4.021324
Zf1=27.29696
Zf2=1.513043
Zf3=3.917305
Zm3=3.330755

Am1=15.36668
Am2=4.220504
Af1=28.66579
Af2=1.588855
Af3=4.108595
Am3=3.496842

N = Nf1+Nf2+Nf3+Nm1+Nm2+Nm3

Sf1=Nf1-Hf1-Yf1-Zf1-Af1
Sf2=Nf2-Hf2-Yf2-Zf2-Af2 
Sf3=Nf3-Hf3-Yf3-Zf3-Af3 

Sm1=Nm1-Hm1-Ym1-Zm1-Am1
Sm2=Nm2-Hm2-Ym2-Zm2-Am2
Sm3=Nm3-Hm3-Ym3-Zm3-Am3

#---------------------------------------------------------------------------------------------#
# Intervention parameters

circ = circ
circ_eff = 0.6 # relative susceptibility of circumcised men = 0.4 = reduction of 60%

int_cond = int_cond

# ART
ART_eff = 0.96 # reduced infectiousness on ART
ART_m = ART_m
ART_f = ART_f

# ART_f=ART_m=0 # for model without ART

ART_m_lv = 0.8*ART_m
ART_m_hv = 0.15*ART_m
ART_m_pa = 0.05*ART_m
ART_f_lv = 0.8*ART_f
ART_f_hv = 0.15*ART_f
ART_f_pa = 0.05*ART_f

ex_ART = 95 # prolonged duration of pre-AIDS stage for those on ART
ex_m1 = (ex_ART*ART_m)
ex_m2 = (ex_ART*ART_m)
ex_m3 = (ex_ART*ART_m)
ex_f1 = (ex_ART*ART_f)
ex_f2 = (ex_ART*ART_f)
ex_f3 = (ex_ART*ART_f)

epi_m1 = epi + ex_m1
epi_m2 = epi + ex_m2
epi_m3 = epi + ex_m3
epi_f1 = epi + ex_f1
epi_f2 = epi + ex_f2
epi_f3 = epi + ex_f3

# PrEP
prep_eff = 0.67 # relative susceptibility on PrEP
prep_m = prep_m
prep_f = prep_f

# Counseling and testing
prop_hct_m = prop_hct_m
prop_hct_f = prop_hct_f

cond_hct = 0.05

circ_hct = 0.05

circ = circ + (circ*circ_hct*prop_hct_m)
circ = ifelse(circ>=1, 1, circ)

eta = 0.1 # reduction in number of partners due to HCT


#------------------------------------------------------------------------------#

## MAIN PARTNERSHIPS
# FEMALES
rho_f1_1_main = (1 - A_main)* (dm1*(Nm1-Am1)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) + A_main 
rho_f1_2_main = (1 - A_main)* (dm2*(Nm2-Am2)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) 
rho_f1_3_main = (1 - A_main)* (dm3*(Nm3-Am3)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) 

rho_f2_1_main = (1 - A_main)* (dm1*(Nm1-Am1)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) 
rho_f2_2_main = (1 - A_main)* (dm2*(Nm2-Am2)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) + A_main 
rho_f2_3_main = (1 - A_main)* (dm3*(Nm3-Am3)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) 

rho_f3_1_main = (1 - A_main)* (dm1*(Nm1-Am1)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) 
rho_f3_2_main = (1 - A_main)* (dm2*(Nm2-Am2)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3)))
rho_f3_3_main = (1 - A_main)* (dm3*(Nm3-Am3)/(dm1*(Nm1-Am1)+dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) + A_main 

# MALES
rho_m1_1_main = (1 - A_main)* (df1*(Nf1-Af1)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) + A_main 
rho_m1_2_main = (1 - A_main)* (df2*(Nf2-Af2)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) 
rho_m1_3_main = (1 - A_main)* (df3*(Nf3-Af3)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) 

rho_m2_1_main = (1 - A_main)* (df1*(Nf1-Af1)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) 
rho_m2_2_main = (1 - A_main)* (df2*(Nf2-Af2)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) + A_main 
rho_m2_3_main = (1 - A_main)* (df3*(Nf3-Af3)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) 

rho_m3_1_main = (1 - A_main)* (df1*(Nf1-Af1)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) 
rho_m3_2_main = (1 - A_main)* (df2*(Nf2-Af2)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) 
rho_m3_3_main = (1 - A_main)* (df3*(Nf3-Af3)/(df1*(Nf1-Af1)+df2*(Nf2-Af2)+df3*(Nf3-Af3))) + A_main 

## REGULAR PARTNERSHIPS
# FEMALES
rho_f2_2_reg = (1 - A_regular)* (dm2*(Nm2-Am2)/(dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) + A_regular 
rho_f2_3_reg = (1 - A_regular)* (dm3*(Nm3-Am3)/(dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) 

rho_f3_2_reg = (1 - A_regular)* (dm2*(Nm2-Am2)/(dm2*(Nm2-Am2)+dm3*(Nm3-Am3)))
rho_f3_3_reg = (1 - A_regular)* (dm3*(Nm3-Am3)/(dm2*(Nm2-Am2)+dm3*(Nm3-Am3))) + A_regular

# MALES
rho_m2_2_reg = (1 - A_regular)* (df2*(Nf2-Af2)/(df2*(Nf2-Af2)+df3*(Nf3-Af3))) + A_regular 
rho_m2_3_reg = (1 - A_regular)* (df3*(Nf3-Af3)/(df2*(Nf2-Af2)+df3*(Nf3-Af3))) 

rho_m3_2_reg = (1 - A_regular)* (df2*(Nf2-Af2)/(df2*(Nf2-Af2)+df3*(Nf3-Af3))) 
rho_m3_3_reg = (1 - A_regular)* (df3*(Nf3-Af3)/(df2*(Nf2-Af2)+df3*(Nf3-Af3))) + A_regular

## CASUAL PARTNERSHIPS (always 1 because A_casual=0)
# FEMALES
rho_f3_3_casual = (1 - A_casual)* (dm3*(Nm3-Am3)/(dm3*(Nm3-Am3))) + A_casual

# MALES
rho_m3_3_casual = (1 - A_casual)* (df3*(Nf3-Af3)/(df3*(Nf3-Af3))) + A_casual

## COMBINED
rho_f1_1 = rho_f1_1_main
rho_f1_2 = rho_f1_2_main
rho_f1_3 = rho_f1_3_main

rho_f2_1 = rho_f2_1_main 
rho_f2_2 = rho_f2_2_main * rho_f2_2_reg
rho_f2_3 = rho_f2_3_main * rho_f2_3_reg

rho_f3_1 = rho_f3_1_main
rho_f3_2 = rho_f3_2_main * rho_f3_2_reg
rho_f3_3 = rho_f3_3_main * rho_f3_3_reg * rho_f3_3_casual

rho_m1_1 = rho_m1_1_main
rho_m1_2 = rho_m1_2_main
rho_m1_3 = rho_m1_3_main

rho_m2_1 = rho_m2_1_main
rho_m2_2 = rho_m2_2_main * rho_m2_2_reg
rho_m2_3 = rho_m2_3_main * rho_m2_3_reg

rho_m3_1 = rho_m3_1_main
rho_m3_2 = rho_m3_2_main * rho_m3_2_reg
rho_m3_3 = rho_m3_3_main * rho_m3_3_reg * rho_m3_3_casual


C_f1=C_f2=C_f3 = C_f  
C_m1=C_m2=C_m3= C_m

acts_f1_main=acts_f1_main - (acts_f1_main*eta*prop_hct_f)
acts_m1_main=acts_m1_main - (acts_m1_main*eta*prop_hct_m)

high_v1 = high_v1
beta_f = beta_f 
beta_m = beta_m

cond_f1 = cond_m1
cond_f1 = cond_f1 + (cond_f1*int_cond)
cond_f1 = cond_f1 + (cond_f1*cond_hct*prop_hct_f)
cond_f1 = ifelse(cond_f1 >= 1, 1, cond_f1)
cond_m1 = cond_m1 + (cond_m1*int_cond)
cond_m1 = cond_m1 + (cond_m1*cond_hct*prop_hct_m)
cond_m1 = ifelse(cond_m1 >= 1, 1, cond_m1)
cond_eff = cond_eff



Gf1 = (1 - high_v1*beta_f*(1-ART_eff*ART_m_hv)*(1-prep_eff*prep_f)*(1-cond_f1*cond_eff))^acts_f1_main 
Gm1 = (1 - high_v1*beta_m*(1-ART_eff*ART_f_hv)*(1-prep_eff*prep_m)*(1-cond_m1*cond_eff)*(1-circ*circ_eff))^acts_m1_main

acts_f2_main=acts_f2_main - (acts_f2_main*eta*prop_hct_f)
acts_m2_main=acts_m2_main - (acts_m2_main*eta*prop_hct_m)
acts_f2_reg=acts_f2_reg - (acts_f2_reg*eta*prop_hct_f)
acts_m2_reg=acts_m2_reg - (acts_m2_reg*eta*prop_hct_m)

cond_f2 = cond_f2 + (cond_f2*int_cond) 
cond_f2 = cond_f2 + (cond_f2*cond_hct*prop_hct_f) # HCT impact
cond_f2 = ifelse(cond_f2 >= 1, 1, cond_f2)
cond_m2 = cond_m2 + (cond_m2*int_cond)
cond_m2 = cond_m2 + (cond_m2*cond_hct*prop_hct_m)
cond_m2 = ifelse(cond_m2 >= 1, 1, cond_m2)

Gf2 = (1 - high_v1*beta_f*(1-ART_eff*ART_m_hv)*(1-prep_eff*prep_f)*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg)  
Gm2 = (1 - high_v1*beta_m*(1-ART_eff*ART_f_hv)*(1-prep_eff*prep_m)*(1-cond_m2*cond_eff)*(1-circ*circ_eff))^(acts_m2_main+acts_m2_reg)
acts_f3_main=acts_f3_main - (acts_f3_main*eta*prop_hct_f)
acts_m3_main=acts_m3_main - (acts_m3_main*eta*prop_hct_m)
acts_f3_casual=acts_f3_casual - (acts_f3_casual*eta*prop_hct_f)
acts_m3_casual=acts_m3_casual - (acts_m3_casual*eta*prop_hct_m)
acts_f3_reg=acts_f3_reg - (acts_f3_reg*eta*prop_hct_f)
acts_m3_reg=acts_m3_reg - (acts_m3_reg*eta*prop_hct_m)

cond_f3 = cond_f3 + (cond_f3*int_cond)
cond_f3 = cond_f3 + (cond_f3*cond_hct*prop_hct_f)
cond_f3 = ifelse(cond_f3 >= 1, 1, cond_f3)
cond_m3 = cond_m3 + (cond_m3*int_cond)
cond_m3 = cond_m3 + (cond_m3*cond_hct*prop_hct_m)
cond_m3 = ifelse(cond_m3 >= 1, 1, cond_m3)

Gf3 = (1 - high_v1*beta_f*(1-ART_eff*ART_m_hv)*(1-prep_eff*prep_f)*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg)  
Gm3 = (1 - high_v1*beta_m*(1-ART_eff*ART_f_hv)*(1-prep_eff*prep_m)*(1-cond_m3*cond_eff)*(1-circ*circ_eff))^(acts_m3_main+acts_m3_casual+acts_m3_reg)



Pf2=(1-beta_f*(1-ART_eff*ART_m_lv)*(1-prep_eff*prep_f)*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg) 
Pm2=(1-beta_m*(1-ART_eff*ART_f_lv)*(1-prep_eff*prep_m)*(1-cond_m2*cond_eff)*(1-circ*circ_eff))^(acts_m2_main+acts_m2_reg) 

Pf1=(1-beta_f*(1-ART_eff*ART_m_lv)*(1-prep_eff*prep_f)*(1-cond_f1*cond_eff))^acts_f1_main 
Pm1=(1-beta_m*(1-ART_eff*ART_f_lv)*(1-prep_eff*prep_m)*(1-cond_m1*cond_eff)*(1-circ*circ_eff))^acts_m1_main 

Pf3=(1-beta_f*(1-ART_eff*ART_m_lv)*(1-prep_eff*prep_f)*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg) 
Pm3=(1-beta_m*(1-ART_eff*ART_f_lv)*(1-prep_eff*prep_m)*(1-cond_m3*cond_eff)*(1-circ*circ_eff))^(acts_m3_main+acts_m3_casual+acts_m3_reg) 

high_v2 = high_v2


Xf2=(1-high_v2*beta_f*(1-ART_eff*ART_m_pa)*(1-prep_eff*prep_f)*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg) 
Xm2=(1-high_v2*beta_m*(1-ART_eff*ART_f_pa)*(1-prep_eff*prep_m)*(1-cond_m2*cond_eff)*(1-circ*circ_eff))^(acts_m2_main+acts_m2_reg) 

Xf1=(1-high_v2*beta_f*(1-ART_eff*ART_m_pa)*(1-prep_eff*prep_f)*(1-cond_f1*cond_eff))^acts_f1_main 
Xm1=(1-high_v2*beta_m*(1-ART_eff*ART_f_pa)*(1-prep_eff*prep_m)*(1-cond_m1*cond_eff)*(1-circ*circ_eff))^acts_m1_main 

Xf3=(1-high_v2*beta_f*(1-ART_eff*ART_m_pa)*(1-prep_eff*prep_f)*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg) 
Xm3=(1-high_v2*beta_m*(1-ART_eff*ART_f_pa)*(1-prep_eff*prep_m)*(1-cond_m3*cond_eff)*(1-circ*circ_eff))^(acts_m3_main+acts_f3_casual+acts_m3_reg) 



D_f1_1 = (1 - (Gf1*Hm1 + Pf1*Ym1 + Xf1*Zm1)/(Ym1 + Hm1 + Zm1))
D_f2_1 = (1 - (Gf2*Hm1 + Pf2*Ym1 + Xf2*Zm1)/(Ym1 + Hm1 + Zm1)) 
D_f3_1 = (1 - (Gf3*Hm1 + Pf3*Ym1 + Xf3*Zm1)/(Ym1 + Hm1 + Zm1)) 

D_f1_2 = (1 - (Gf1*Hm2 + Pf1*Ym2 + Xf1*Zm2)/(Ym2 + Hm2 + Zm2)) 
D_f2_2 = (1 - (Gf2*Hm2 + Pf2*Ym2 + Xf2*Zm2)/(Ym2 + Hm2 + Zm2)) 
D_f3_2 = (1 - (Gf3*Hm2 + Pf3*Ym2 + Xf3*Zm2)/(Ym2 + Hm2 + Zm2)) 

D_f1_3 = (1 - (Gf1*Hm3 + Pf1*Ym3 + Xf1*Zm3)/(Ym3 + Hm3 + Zm3)) 
D_f2_3 = (1 - (Gf2*Hm3 + Pf2*Ym3 + Xf2*Zm3)/(Ym3 + Hm3 + Zm3)) 
D_f3_3 = (1 - (Gf3*Hm3 + Pf3*Ym3 + Xf3*Zm3)/(Ym3 + Hm3 + Zm3)) 

# male 
D_m1_1 = (1 - (Gm1*Hf1 + Pm1*Yf1 + Xm1*Zf1)/(Yf1 + Hf1 + Zf1)) 
D_m2_1 = (1 - (Gm2*Hf1 + Pm2*Yf1 + Xm2*Zf1)/(Yf1 + Hf1 + Zf1)) 
D_m3_1 = (1 - (Gm3*Hf1 + Pm3*Yf1 + Xm3*Zf1)/(Yf1 + Hf1 + Zf1)) 

D_m1_2 = (1 - (Gm1*Hf2 + Pm1*Yf2 + Xm1*Zf2)/(Yf2 + Hf2 + Zf2)) 
D_m2_2 = (1 - (Gm2*Hf2 + Pm2*Yf2 + Xm2*Zf2)/(Yf2 + Hf2 + Zf2)) 
D_m3_2 = (1 - (Gm3*Hf2 + Pm3*Yf2 + Xm3*Zf2)/(Yf2 + Hf2 + Zf2)) 

D_m1_3 = (1 - (Gm1*Hf3 + Pm1*Yf3 + Xm1*Zf3)/(Yf3 + Hf3 + Zf3)) 
D_m2_3 = (1 - (Gm2*Hf3 + Pm2*Yf3 + Xm2*Zf3)/(Yf3 + Hf3 + Zf3)) 
D_m3_3 = (1 - (Gm3*Hf3 + Pm3*Yf3 + Xm3*Zf3)/(Yf3 + Hf3 + Zf3)) 


phi_f1_1 = C_f1*D_f1_1 
phi_f1_2 = C_f1*D_f1_2 
phi_f1_3 = C_f1*D_f1_3 

phi_f2_1 = C_f2*D_f2_1 
phi_f2_2 = C_f2*D_f2_2 
phi_f2_3 = C_f2*D_f2_3 

phi_f3_1 = C_f3*D_f3_1 
phi_f3_2 = C_f3*D_f3_2 
phi_f3_3 = C_f3*D_f3_3 

# male
phi_m1_1 = C_m1*D_m1_1 
phi_m1_2 = C_m1*D_m1_2 
phi_m1_3 = C_m1*D_m1_3 

phi_m2_1 = C_m2*D_m2_1 
phi_m2_2 = C_m2*D_m2_2 
phi_m2_3 = C_m2*D_m2_3 

phi_m3_1 = C_m3*D_m3_1 
phi_m3_2 = C_m3*D_m3_2 
phi_m3_3 = C_m3*D_m3_3 


pi_f1= 1 - ((1-phi_f1_1)^(df1*rho_f1_1)*(1-phi_f1_2)^(df1*rho_f1_2)*(1-phi_f1_3)^(df1*rho_f1_3)) 
pi_f2= 1 - ((1-phi_f2_1)^(df2*rho_f2_1)*(1-phi_f2_2)^(df2*rho_f2_2)*(1-phi_f2_3)^(df2*rho_f2_3)) 
pi_f3= 1 - ((1-phi_f3_1)^(df3*rho_f3_1)*(1-phi_f3_2)^(df3*rho_f3_2)*(1-phi_f3_3)^(df3*rho_f3_3)) 

pi_m1= 1 - ((1-phi_m1_1)^(dm1*rho_m1_1)*(1-phi_m1_2)^(dm1*rho_m1_2)*(1-phi_m1_3)^(dm1*rho_m1_3)) 
pi_m2= 1 - ((1-phi_m2_1)^(dm2*rho_m2_1)*(1-phi_m2_2)^(dm2*rho_m2_2)*(1-phi_m2_3)^(dm2*rho_m2_3)) 
pi_m3= 1 - ((1-phi_m3_1)^(dm3*rho_m3_1)*(1-phi_m3_2)^(dm3*rho_m3_2)*(1-phi_m3_3)^(dm3*rho_m3_3)) 

state=c(Sf1=Sf1,Hf1=Hf1,Yf1=Yf1,Zf1=Zf1,Af1=Af1,
        Sf2=Sf2,Hf2=Hf2,Yf2=Yf2,Zf2=Zf2,Af2=Af2,
        Sf3=Sf3,Hf3=Hf3,Yf3=Yf3,Zf3=Zf3,Af3=Af3,
        Sm1=Sm1,Hm1=Hm1,Ym1=Ym1,Zm1=Zm1,Am1=Am1,
        Sm2=Sm2,Hm2=Hm2,Ym2=Ym2,Zm2=Zm2,Am2=Am2,
        Sm3=Sm3,Hm3=Hm3,Ym3=Ym3,Zm3=Zm3,Am3=Am3)

times=seq(0,180,by=0.01) 

parameters=c(nu=1/nu, lambda=lambda, psi=psi,omega=1/omega, epi_m1=1/epi_m1, epi_m2=1/epi_m2, epi_m3=1/epi_m3, epi_f1=1/epi_f1, epi_f2=1/epi_f2, epi_f3=1/epi_f3,
             delta=1/delta,pi_f1=pi_f1,pi_m1=pi_m1,pi_f2=pi_f2,pi_m2=pi_m2,pi_f3=pi_f3,pi_m3=pi_m3)

require(deSolve)

trial1_FM=as.data.frame(ode(y=state,times=times,func=SHYZA_FM,parms=parameters))

N_end_f1=sum(trial1_FM[18001,2:6])
N_end_f2=sum(trial1_FM[18001,7:11])
N_end_f3=sum(trial1_FM[18001,12:16])
N_end_m1=sum(trial1_FM[18001,17:21])
N_end_m2=sum(trial1_FM[18001,22:26])
N_end_m3=sum(trial1_FM[18001,27:31])

# Number of cases
HIV_f1=(trial1_FM[18001,6]+trial1_FM[18001,3]+trial1_FM[18001,4]+trial1_FM[18001,5]) 
HIV_f2=(trial1_FM[18001,11]+trial1_FM[18001,8]+trial1_FM[18001,9]+trial1_FM[18001,10]) 
HIV_f3=(trial1_FM[18001,16]+trial1_FM[18001,13]+trial1_FM[18001,14]+trial1_FM[18001,15]) 
HIV_m1=(trial1_FM[18001,21]+trial1_FM[18001,18]+trial1_FM[18001,19]+trial1_FM[18001,20]) 
HIV_m2=(trial1_FM[18001,26]+trial1_FM[18001,23]+trial1_FM[18001,24]+trial1_FM[18001,25]) 
HIV_m3=(trial1_FM[18001,31]+trial1_FM[18001,28]+trial1_FM[18001,29]+trial1_FM[18001,30])

HIV_all=(sum(trial1_FM[18001,3:6])+sum(trial1_FM[18001,8:11])+sum(trial1_FM[18001,13:16])+sum(trial1_FM[18001,18:21])
         +sum(trial1_FM[18001,23:26])
         +sum(trial1_FM[18001,28:31]))

#Proportion of HIV positives in each class:
prev_f1=(trial1_FM[18001,6]+trial1_FM[18001,3]+trial1_FM[18001,4]+trial1_FM[18001,5])/N_end_f1 # 50.33%
prev_f2=(trial1_FM[18001,11]+trial1_FM[18001,8]+trial1_FM[18001,9]+trial1_FM[18001,10])/N_end_f2 # 50.60%
prev_f3=(trial1_FM[18001,16]+trial1_FM[18001,13]+trial1_FM[18001,14]+trial1_FM[18001,15])/N_end_f3 # 55.41%
prev_m1=(trial1_FM[18001,21]+trial1_FM[18001,18]+trial1_FM[18001,19]+trial1_FM[18001,20])/N_end_m1 # 74.36%
prev_m2=(trial1_FM[18001,26]+trial1_FM[18001,23]+trial1_FM[18001,24]+trial1_FM[18001,25])/N_end_m2 # 79.92%
prev_m3=(trial1_FM[18001,31]+trial1_FM[18001,28]+trial1_FM[18001,29]+trial1_FM[18001,30])/N_end_m3 # 51.73%

# HIV prevalence in the total population
N_total_end=sum(trial1_FM[18001,2:31])
HIV_prev=(sum(trial1_FM[18001,3:6])+sum(trial1_FM[18001,8:11])+sum(trial1_FM[18001,13:16])+sum(trial1_FM[18001,18:21])
          +sum(trial1_FM[18001,23:26])
          +sum(trial1_FM[18001,28:31]))/N_total_end
prev_all=round(HIV_prev*100,2)

# Relative reduction in prevalence
RR_f1 = (0.2122097 - prev_f1) / 0.2122097
RR_f2 = (0.3240974 - prev_f2) / 0.3240974
RR_f3 = (0.3649333 - prev_f3) / 0.3649333
RR_m1 = (0.1181192 - prev_m1) / 0.1181192
RR_m2 = (0.175524 - prev_m2) / 0.175524
RR_m3 = (0.1752642 - prev_m3) / 0.1752642
RR_all = (0.1765018 - HIV_prev) / 0.1765018

# Difference number of infections
cum_f1 = 173.1674 - HIV_f1
cum_f2 = 9.601216 - HIV_f2
cum_f3 = 24.87418 - HIV_f3
cum_m1 = 93.19444 - HIV_m1
cum_m2 = 25.51623 - HIV_m2
cum_m3 = 21.13822 - HIV_m3
cum_all = 347.4917 - HIV_all

### OUTPUT ###

prev=data.frame(sim=sim, N_total=round(N_total_end,0), prev_all=round(prev_all,2), prev_f1=round(prev_f1*100,2), prev_f2=round(prev_f2*100,2), 
                prev_f3=round(prev_f3*100,2), prev_m1=round(prev_m1*100,2), prev_m2=round(prev_m2*100,2), prev_m3=round(prev_m3*100,2),
                RR_f1=RR_f1, RR_f2=RR_f2, RR_f3=RR_f3, RR_m1=RR_m1, RR_m2=RR_m2, RR_m3=RR_m3, RR_all=RR_all,
                cum_f1=cum_f1, cum_f2=cum_f2, cum_f3=cum_f3, cum_m1=cum_m1, cum_m2=cum_m2, cum_m3=cum_m3, cum_all=cum_all)

write.table(prev,file=paste0("output_sim_", sim, ".csv"),append=F,col.names=FALSE,row.names=FALSE,sep=",")




