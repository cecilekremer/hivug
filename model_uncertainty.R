args <- commandArgs(trailingOnly =  TRUE)
lambda <- as.numeric(args[1])
psi <- as.numeric(args[2])
C_f <- as.numeric(args[3])
C_m <- as.numeric(args[4])
high_v1 <- as.numeric(args[5])
high_v2 <- as.numeric(args[6])
beta_f <- as.numeric(args[7])
beta_m <- as.numeric(args[8])
cond_m1 <- as.numeric(args[9])
cond_m2 <- as.numeric(args[10])
cond_m3 <- as.numeric(args[11])
cond_f2 <- as.numeric(args[12])
cond_f3 <- as.numeric(args[13])
cond_eff <- as.numeric(args[14])
acts_m3_casual <- as.numeric(args[15])
acts_f3_casual <- as.numeric(args[16])
acts_m2_reg <- as.numeric(args[17])
acts_m3_reg <- as.numeric(args[18])
acts_f2_reg <- as.numeric(args[19])
acts_f3_reg <- as.numeric(args[20])
acts_m1_main <- as.numeric(args[21])
acts_m2_main <- as.numeric(args[22])
acts_m3_main <- as.numeric(args[23])
acts_f1_main <- as.numeric(args[24])
acts_f2_main <- as.numeric(args[25])
acts_f3_main <- as.numeric(args[26])
nu <- as.numeric(args[27])
omega <- as.numeric(args[28])
epi <- as.numeric(args[29])
delta <- as.numeric(args[30])
dm1 <- as.numeric(args[31])
dm2 <- as.numeric(args[32])
dm3 <- as.numeric(args[33])
df1 <- as.numeric(args[34])
df2 <- as.numeric(args[35])
df3 <- as.numeric(args[36])
N <- as.integer(args[37])
A_main <- as.numeric(args[38])
A_regular <- as.numeric(args[39])
A_casual <- as.numeric(args[40])
int_cond <- as.numeric(args[41])
circ <- as.numeric(args[42])
prep_m <- as.numeric(args[43])
prep_f <- as.numeric(args[44])
ART_m <- as.numeric(args[45])
ART_f <- as.numeric(args[46])
prop_hct_m <- as.numeric(args[47])
prop_hct_f <- as.numeric(args[48])
sim <- as.numeric(args[49])


library(deSolve)

#######################
### CREATE BASELINE ###
#######################

### ODES
SHYZA_FM=function(t, state, parameters)
{
  with(as.list(c(state,parameters)),
       {
         ### Females       
         dSf1 = lambda*Nf1 - (pi_f1 + psi)*Sf1
         dHf1 = pi_f1*Sf1 - (nu+psi)*Hf1
         dYf1 = nu*Hf1 - (omega+psi)*Yf1
         dZf1 = omega*Yf1-(epi+psi)*Zf1
         dAf1 = epi*Zf1-(delta+psi)*Af1
         
         dSf2 = lambda*Nf2 - (pi_f2 + psi)*Sf2
         dHf2 = pi_f2*Sf2 - (nu+psi)*Hf2
         dYf2 = nu*Hf2 - (omega+psi)*Yf2
         dZf2 = omega*Yf2-(epi+psi)*Zf2
         dAf2 = epi*Zf2-(delta+psi)*Af2
         
         dSf3 = lambda*Nf3 - (pi_f3 + psi)*Sf3
         dHf3 = pi_f3*Sf3 - (nu+psi)*Hf3
         dYf3 = nu*Hf3 - (omega+psi)*Yf3
         dZf3 = omega*Yf3-(epi+psi)*Zf3
         dAf3 = epi*Zf3-(delta+psi)*Af3
         
         ### Males
         dSm1 = lambda*Nm1 - (pi_m1 + psi)*Sm1
         dHm1 = pi_m1*Sm1 - (nu+psi)*Hm1
         dYm1 = nu*Hm1 - (omega+psi)*Ym1
         dZm1 = omega*Ym1-(epi+psi)*Zm1
         dAm1 = epi*Zm1-(delta+psi)*Am1
         
         dSm2 = lambda*Nm2 - (pi_m2 + psi)*Sm2
         dHm2 = pi_m2*Sm2 - (nu+psi)*Hm2
         dYm2 = nu*Hm2 - (omega+psi)*Ym2
         dZm2 = omega*Ym2-(epi+psi)*Zm2
         dAm2 = epi*Zm2-(delta+psi)*Am2
         
         dSm3 = lambda*Nm3 - (pi_m3 + psi)*Sm3
         dHm3 = pi_m2*Sm3 - (nu+psi)*Hm3
         dYm3 = nu*Hm3 - (omega+psi)*Ym3
         dZm3 = omega*Ym3-(epi+psi)*Zm3
         dAm3 = epi*Zm3-(delta+psi)*Am3
         
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
N=N

Nf1=0.43*N; Nm1=0.33*N
Nf2=0.02*N; Nm2=0.07*N
Nf3=0.05*N ; Nm3=0.1*N 

Hm1=floor(0.4*0.17*Nm1)
Hm2=floor(0.4*0.17*Nm2)
Hf1=floor(0.4*0.17*Nf1)
Hf2=floor(0.4*0.17*Nf2)
Hf3=floor(0.4*0.17*Nf3)
Hm3=floor(0.4*0.17*Nm3)

Ym1=floor(0.35*0.17*Nm1)
Ym2=floor(0.35*0.17*Nm2)
Yf1=floor(0.35*0.17*Nf1)
Yf2=floor(0.35*0.17*Nf2)
Yf3=floor(0.35*0.17*Nf3)
Ym3=floor(0.35*0.17*Nm3)

Zm1=floor(0.2*0.17*Nm1)
Zm2=floor(0.2*0.17*Nm2)
Zf1=floor(0.2*0.17*Nf1)
Zf2=floor(0.2*0.17*Nf2)
Zf3=floor(0.2*0.17*Nf3)
Zm3=floor(0.2*0.17*Nm3)

Am1=floor(0.05*0.17*Nm1)
Am2=floor(0.05*0.17*Nm2)
Af1=floor(0.05*0.17*Nf1)
Af2=floor(0.05*0.17*Nf2)
Af3=floor(0.05*0.17*Nf3)
Am3=floor(0.05*0.17*Nm3)

Sf1=Nf1-Hf1-Yf1-Zf1-Af1
Sf2=Nf2-Hf2-Yf2-Zf2-Af2 
Sf3=Nf3-Hf3-Yf3-Zf3-Af3 

Sm1=Nm1-Hm1-Ym1-Zm1-Am1
Sm2=Nm2-Hm2-Ym2-Zm2-Am2
Sm3=Nm3-Hm3-Ym3-Zm3-Am3


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

acts_f1_main=acts_f1_main
acts_m1_main=acts_m1_main
high_v1 = high_v1
beta_f = beta_f 
beta_m = beta_m
cond_f1 = cond_m1 = cond_m1
cond_eff = cond_eff

Gf1 = (1 - high_v1*beta_f*(1-cond_f1*cond_eff))^acts_f1_main 
Gm1 = (1 - high_v1*beta_m*(1-cond_m1*cond_eff))^acts_m1_main 

acts_f2_main=acts_f2_main
acts_m2_main=acts_m2_main
acts_f2_reg=acts_f2_reg
acts_m2_reg=acts_m2_reg
cond_f2 = cond_f2
cond_m2 = cond_m2

Gf2 = (1 - high_v1*beta_f*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg)  
Gm2 = (1 - high_v1*beta_m*(1-cond_m2*cond_eff))^(acts_m2_main+acts_m2_reg) 

acts_f3_main=acts_f3_main
acts_m3_main=acts_m3_main
cond_f3 = cond_f3
cond_m3 = cond_m3
acts_f3_casual=acts_f3_casual
acts_m3_casual=acts_m3_casual
acts_f3_reg=acts_f3_reg
acts_m3_reg=acts_m3_reg

Gf3 = (1 - high_v1*beta_f*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg)  
Gm3 = (1 - high_v1*beta_m*(1-cond_m3*cond_eff))^(acts_m3_main+acts_m3_casual+acts_m3_reg)


Pf2=(1-beta_f*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg) 
Pm2=(1-beta_m*(1-cond_m2*cond_eff))^(acts_m2_main+acts_m2_reg) 

Pf1=(1-beta_f*(1-cond_f1*cond_eff))^acts_f1_main 
Pm1=(1-beta_m*(1-cond_m1*cond_eff))^acts_m1_main 

Pf3=(1-beta_f*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg) 
Pm3=(1-beta_m*(1-cond_m3*cond_eff))^(acts_m3_main+acts_m3_casual+acts_m3_reg) 

high_v2 = high_v2

Xf2=(1-high_v2*beta_f*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg) 
Xm2=(1-high_v2*beta_m*(1-cond_m2*cond_eff))^(acts_m2_main+acts_m2_reg) 

Xf1=(1-high_v2*beta_f*(1-cond_f1*cond_eff))^acts_f1_main 
Xm1=(1-high_v2*beta_m*(1-cond_m1*cond_eff))^acts_m1_main 

Xf3=(1-high_v2*beta_f*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg) 
Xm3=(1-high_v2*beta_m*(1-cond_m3*cond_eff))^(acts_m3_main+acts_f3_casual+acts_m3_reg) 




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

times=seq(0,780,by=0.01) 

parameters=c(nu=1/nu, lambda=lambda, psi=psi,omega=1/omega, epi=1/epi, delta=1/delta,pi_f1=pi_f1,pi_m1=pi_m1,pi_f2=pi_f2,pi_m2=pi_m2,pi_f3=pi_f3,pi_m3=pi_m3)

require(deSolve)

trial1_FM=as.data.frame(ode(y=state,times=times,func=SHYZA_FM,parms=parameters))

N_end_f1=sum(trial1_FM[60001,2:6])
N_end_f2=sum(trial1_FM[60001,7:11])
N_end_f3=sum(trial1_FM[60001,12:16])
N_end_m1=sum(trial1_FM[60001,17:21])
N_end_m2=sum(trial1_FM[60001,22:26])
N_end_m3=sum(trial1_FM[60001,27:31])

H_end_f1=sum(trial1_FM[60001,3]); Y_end_f1=sum(trial1_FM[60001,4]); Z_end_f1=sum(trial1_FM[60001,5]); A_end_f1=sum(trial1_FM[60001,6])
H_end_f2=sum(trial1_FM[60001,8]); Y_end_f2=sum(trial1_FM[60001,9]); Z_end_f2=sum(trial1_FM[60001,10]); A_end_f2=sum(trial1_FM[60001,11])
H_end_f3=sum(trial1_FM[60001,13]); Y_end_f3=sum(trial1_FM[60001,14]); Z_end_f3=sum(trial1_FM[60001,15]); A_end_f3=sum(trial1_FM[60001,16])
H_end_m1=sum(trial1_FM[60001,18]); Y_end_m1=sum(trial1_FM[60001,19]); Z_end_m1=sum(trial1_FM[60001,20]); A_end_m1=sum(trial1_FM[60001,21])
H_end_m2=sum(trial1_FM[60001,23]); Y_end_m2=sum(trial1_FM[60001,24]); Z_end_m2=sum(trial1_FM[60001,25]); A_end_m2=sum(trial1_FM[60001,26])
H_end_m3=sum(trial1_FM[60001,28]); Y_end_m3=sum(trial1_FM[60001,29]); Z_end_m3=sum(trial1_FM[60001,30]); A_end_m3=sum(trial1_FM[60001,31])

# Number of cases
HIV_f1=(trial1_FM[60001,6]+trial1_FM[60001,3]+trial1_FM[60001,4]+trial1_FM[60001,5]) 
HIV_f2=(trial1_FM[60001,11]+trial1_FM[60001,8]+trial1_FM[60001,9]+trial1_FM[60001,10]) 
HIV_f3=(trial1_FM[60001,16]+trial1_FM[60001,13]+trial1_FM[60001,14]+trial1_FM[60001,15]) 
HIV_m1=(trial1_FM[60001,21]+trial1_FM[60001,18]+trial1_FM[60001,19]+trial1_FM[60001,20]) 
HIV_m2=(trial1_FM[60001,26]+trial1_FM[60001,23]+trial1_FM[60001,24]+trial1_FM[60001,25]) 
HIV_m3=(trial1_FM[60001,31]+trial1_FM[60001,28]+trial1_FM[60001,29]+trial1_FM[60001,30])

HIV_all=(sum(trial1_FM[60001,3:6])+sum(trial1_FM[60001,8:11])+sum(trial1_FM[60001,13:16])+sum(trial1_FM[60001,18:21])
         +sum(trial1_FM[60001,23:26])
         +sum(trial1_FM[60001,28:31]))

#Proportion of HIV positives in each class:
prev_f1=(trial1_FM[60001,6]+trial1_FM[60001,3]+trial1_FM[60001,4]+trial1_FM[60001,5])/N_end_f1 # 50.33%
prev_f2=(trial1_FM[60001,11]+trial1_FM[60001,8]+trial1_FM[60001,9]+trial1_FM[60001,10])/N_end_f2 # 50.60%
prev_f3=(trial1_FM[60001,16]+trial1_FM[60001,13]+trial1_FM[60001,14]+trial1_FM[60001,15])/N_end_f3 # 55.41%
prev_m1=(trial1_FM[60001,21]+trial1_FM[60001,18]+trial1_FM[60001,19]+trial1_FM[60001,20])/N_end_m1 # 74.36%
prev_m2=(trial1_FM[60001,26]+trial1_FM[60001,23]+trial1_FM[60001,24]+trial1_FM[60001,25])/N_end_m2 # 79.92%
prev_m3=(trial1_FM[60001,31]+trial1_FM[60001,28]+trial1_FM[60001,29]+trial1_FM[60001,30])/N_end_m3 # 51.73%

# HIV prevalence in the total population
N_total_end=sum(trial1_FM[60001,2:31])
HIV_prev=(sum(trial1_FM[60001,3:6])+sum(trial1_FM[60001,8:11])+sum(trial1_FM[60001,13:16])+sum(trial1_FM[60001,18:21])
          +sum(trial1_FM[60001,23:26])
          +sum(trial1_FM[60001,28:31]))/N_total_end
prev_all=round(HIV_prev*100,2)

############################################
## Baseline (15 years extra no intervention)
N_end_f1_base=sum(trial1_FM[78001,2:6])
N_end_f2_base=sum(trial1_FM[78001,7:11])
N_end_f3_base=sum(trial1_FM[78001,12:16])
N_end_m1_base=sum(trial1_FM[78001,17:21])
N_end_m2_base=sum(trial1_FM[78001,22:26])
N_end_m3_base=sum(trial1_FM[78001,27:31])

H_end_f1_base=sum(trial1_FM[78001,3]); Y_end_f1_base=sum(trial1_FM[78001,4]); Z_end_f1_base=sum(trial1_FM[78001,5]); A_end_f1_base=sum(trial1_FM[78001,6])
H_end_f2_base=sum(trial1_FM[78001,8]); Y_end_f2_base=sum(trial1_FM[78001,9]); Z_end_f2_base=sum(trial1_FM[78001,10]); A_end_f2_base=sum(trial1_FM[78001,11])
H_end_f3_base=sum(trial1_FM[78001,13]); Y_end_f3_base=sum(trial1_FM[78001,14]); Z_end_f3_base=sum(trial1_FM[78001,15]); A_end_f3_base=sum(trial1_FM[78001,16])
H_end_m1_base=sum(trial1_FM[78001,18]); Y_end_m1_base=sum(trial1_FM[78001,19]); Z_end_m1_base=sum(trial1_FM[78001,20]); A_end_m1_base=sum(trial1_FM[78001,21])
H_end_m2_base=sum(trial1_FM[78001,23]); Y_end_m2_base=sum(trial1_FM[78001,24]); Z_end_m2_base=sum(trial1_FM[78001,25]); A_end_m2_base=sum(trial1_FM[78001,26])
H_end_m3_base=sum(trial1_FM[78001,28]); Y_end_m3_base=sum(trial1_FM[78001,29]); Z_end_m3_base=sum(trial1_FM[78001,30]); A_end_m3_base=sum(trial1_FM[78001,31])

# Number of cases
HIV_f1_base=(trial1_FM[78001,6]+trial1_FM[78001,3]+trial1_FM[78001,4]+trial1_FM[78001,5]) 
HIV_f2_base=(trial1_FM[78001,11]+trial1_FM[78001,8]+trial1_FM[78001,9]+trial1_FM[78001,10]) 
HIV_f3_base=(trial1_FM[78001,16]+trial1_FM[78001,13]+trial1_FM[78001,14]+trial1_FM[78001,15]) 
HIV_m1_base=(trial1_FM[78001,21]+trial1_FM[78001,18]+trial1_FM[78001,19]+trial1_FM[78001,20]) 
HIV_m2_base=(trial1_FM[78001,26]+trial1_FM[78001,23]+trial1_FM[78001,24]+trial1_FM[78001,25]) 
HIV_m3_base=(trial1_FM[78001,31]+trial1_FM[78001,28]+trial1_FM[78001,29]+trial1_FM[78001,30])

HIV_all_base=(sum(trial1_FM[78001,3:6])+sum(trial1_FM[78001,8:11])+sum(trial1_FM[78001,13:16])+sum(trial1_FM[78001,18:21])
              +sum(trial1_FM[78001,23:26])
              +sum(trial1_FM[78001,28:31]))

#Proportion of HIV positives in each class:
prev_f1_base=(trial1_FM[78001,6]+trial1_FM[78001,3]+trial1_FM[78001,4]+trial1_FM[78001,5])/N_end_f1_base # 50.33%
prev_f2_base=(trial1_FM[78001,11]+trial1_FM[78001,8]+trial1_FM[78001,9]+trial1_FM[78001,10])/N_end_f2_base # 50.60%
prev_f3_base=(trial1_FM[78001,16]+trial1_FM[78001,13]+trial1_FM[78001,14]+trial1_FM[78001,15])/N_end_f3_base # 55.41%
prev_m1_base=(trial1_FM[78001,21]+trial1_FM[78001,18]+trial1_FM[78001,19]+trial1_FM[78001,20])/N_end_m1_base # 74.36%
prev_m2_base=(trial1_FM[78001,26]+trial1_FM[78001,23]+trial1_FM[78001,24]+trial1_FM[78001,25])/N_end_m2_base # 79.92%
prev_m3_base=(trial1_FM[78001,31]+trial1_FM[78001,28]+trial1_FM[78001,29]+trial1_FM[78001,30])/N_end_m3_base # 51.73%

# HIV prevalence in the total population
N_total_end_base=sum(trial1_FM[78001,2:31])
HIV_prev_base=(sum(trial1_FM[78001,3:6])+sum(trial1_FM[78001,8:11])+sum(trial1_FM[78001,13:16])+sum(trial1_FM[78001,18:21])
               +sum(trial1_FM[78001,23:26])
               +sum(trial1_FM[78001,28:31]))/N_total_end_base
prev_all_base=round(HIV_prev_base*100,2)


##################################
### RUN INTERVENTION SCENARIOS ###
##################################


### ODES

SHYZA_FM_INT=function(t, state, parameters)
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


### PARAMETERS


# rate of partner change
dm1=dm1; dm2=dm2; df1=df1; df2=df2
dm3=dm3 ; df3=df3

# assortativity parameters
A_main=A_main  
A_regular=A_regular
A_casual=A_casual

# Take states at end of 50 years (calibration)
Nf1=N_end_f1
Nf2=N_end_f2
Nf3=N_end_f3
Nm1=N_end_m1
Nm2=N_end_m2
Nm3=N_end_m3

Hm1=H_end_m1
Hm2=H_end_m2
Hf1=H_end_f1
Hf2=H_end_f2
Hf3=H_end_f3
Hm3=H_end_m3

Ym1=Y_end_m1
Ym2=Y_end_m2
Yf1=Y_end_f1
Yf2=Y_end_f2
Yf3=Y_end_f3
Ym3=Y_end_m3

Zm1=Z_end_m1
Zm2=Z_end_m2
Zf1=Z_end_f1
Zf2=Z_end_f2
Zf3=Z_end_f3
Zm3=Z_end_m3

Am1=A_end_m1
Am2=A_end_m2
Af1=A_end_f1
Af2=A_end_f2
Af3=A_end_f3
Am3=A_end_m3

N = Nf1+Nf2+Nf3+Nm1+Nm2+Nm3

Sf1=Nf1-Hf1-Yf1-Zf1-Af1
Sf2=Nf2-Hf2-Yf2-Zf2-Af2 
Sf3=Nf3-Hf3-Yf3-Zf3-Af3 

Sm1=Nm1-Hm1-Ym1-Zm1-Am1
Sm2=Nm2-Hm2-Ym2-Zm2-Am2
Sm3=Nm3-Hm3-Ym3-Zm3-Am3

#---------------------------------------------------------------------------------------------#
# Intervention parameters

circ_eff = 0.6 # relative susceptibility of circumcised men = 0.4 = reduction of 60%

ART_eff = 0.96 # reduced infectiousness on ART

# proportion on ART in each disease stage
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

#duration in pre-AIDS stage lifted upwards by the proportion on ART
epi_m1 = epi + ex_m1
epi_m2 = epi + ex_m2
epi_m3 = epi + ex_m3
epi_f1 = epi + ex_f1
epi_f2 = epi + ex_f2
epi_f3 = epi + ex_f3

prep_eff = 0.67 # relative susceptibility on PrEP

# Counseling and testing
cond_hct = 0.05 # increase in condom use
circ_hct = 0.05 # increase in circumcision
eta = 0.1 # reduction in number of partners due to HCT

#total proportion circumcised men
circ = circ + (circ*circ_hct*prop_hct_m)
circ = ifelse(circ>=1, 1, circ)

#------------------------------------------------------------------------------#

## MIXING MATRIX
# MAIN PARTNERSHIPS
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

# REGULAR PARTNERSHIPS
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

# CASUAL PARTNERSHIPS 
# FEMALES
rho_f3_3_casual = (1 - A_casual)* (dm3*(Nm3-Am3)/(dm3*(Nm3-Am3))) + A_casual

# MALES
rho_m3_3_casual = (1 - A_casual)* (df3*(Nf3-Af3)/(df3*(Nf3-Af3))) + A_casual

## COMBINED MIXING MATRIX
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


## Probability of acquiring infection upon exposure
C_f1=C_f2=C_f3 = C_f  
C_m1=C_m2=C_m3= C_m

## Number of sex acts per unit time
acts_f1_main=acts_f1_main - (acts_f1_main*eta*prop_hct_f)
acts_m1_main=acts_m1_main - (acts_m1_main*eta*prop_hct_m)

acts_f2_main=acts_f2_main - (acts_f2_main*eta*prop_hct_f)
acts_m2_main=acts_m2_main - (acts_m2_main*eta*prop_hct_m)
acts_f2_reg=acts_f2_reg - (acts_f2_reg*eta*prop_hct_f)
acts_m2_reg=acts_m2_reg - (acts_m2_reg*eta*prop_hct_m)

acts_f3_main=acts_f3_main - (acts_f3_main*eta*prop_hct_f)
acts_m3_main=acts_m3_main - (acts_m3_main*eta*prop_hct_m)
acts_f3_casual=acts_f3_casual - (acts_f3_casual*eta*prop_hct_f)
acts_m3_casual=acts_m3_casual - (acts_m3_casual*eta*prop_hct_m)
acts_f3_reg=acts_f3_reg - (acts_f3_reg*eta*prop_hct_f)
acts_m3_reg=acts_m3_reg - (acts_m3_reg*eta*prop_hct_m)

## Condom use
cond_eff = cond_eff

cond_f1 = cond_m1
cond_f1 = cond_f1 + (cond_f1*int_cond)
cond_f1 = cond_f1 + (cond_f1*cond_hct*prop_hct_f)
cond_f1 = ifelse(cond_f1 >= 1, 1, cond_f1)
cond_m1 = cond_m1 + (cond_m1*int_cond)
cond_m1 = cond_m1 + (cond_m1*cond_hct*prop_hct_m)
cond_m1 = ifelse(cond_m1 >= 1, 1, cond_m1)

cond_f2 = cond_f2 + (cond_f2*int_cond) 
cond_f2 = cond_f2 + (cond_f2*cond_hct*prop_hct_f) 
cond_f2 = ifelse(cond_f2 >= 1, 1, cond_f2)
cond_m2 = cond_m2 + (cond_m2*int_cond)
cond_m2 = cond_m2 + (cond_m2*cond_hct*prop_hct_m)
cond_m2 = ifelse(cond_m2 >= 1, 1, cond_m2)

cond_f3 = cond_f3 + (cond_f3*int_cond)
cond_f3 = cond_f3 + (cond_f3*cond_hct*prop_hct_f)
cond_f3 = ifelse(cond_f3 >= 1, 1, cond_f3)
cond_m3 = cond_m3 + (cond_m3*int_cond)
cond_m3 = cond_m3 + (cond_m3*cond_hct*prop_hct_m)
cond_m3 = ifelse(cond_m3 >= 1, 1, cond_m3)

## Transmission parameters
high_v1 = high_v1 
high_v2 = high_v2
beta_f = beta_f 
beta_m = beta_m

## Probabilities of not being exposed to HIV when having sex with an infected partner
#partner in initial high viremia stage
Gf1 = (1 - high_v1*beta_f*(1-ART_eff*ART_m_hv)*(1-prep_eff*prep_f)*(1-cond_f1*cond_eff))^acts_f1_main 
Gm1 = (1 - high_v1*beta_m*(1-ART_eff*ART_f_hv)*(1-prep_eff*prep_m)*(1-cond_m1*cond_eff)*(1-circ*circ_eff))^acts_m1_main
Gf2 = (1 - high_v1*beta_f*(1-ART_eff*ART_m_hv)*(1-prep_eff*prep_f)*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg)  
Gm2 = (1 - high_v1*beta_m*(1-ART_eff*ART_f_hv)*(1-prep_eff*prep_m)*(1-cond_m2*cond_eff)*(1-circ*circ_eff))^(acts_m2_main+acts_m2_reg)
Gf3 = (1 - high_v1*beta_f*(1-ART_eff*ART_m_hv)*(1-prep_eff*prep_f)*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg)  
Gm3 = (1 - high_v1*beta_m*(1-ART_eff*ART_f_hv)*(1-prep_eff*prep_m)*(1-cond_m3*cond_eff)*(1-circ*circ_eff))^(acts_m3_main+acts_m3_casual+acts_m3_reg)
#partner in low viremia stage
Pf2=(1-beta_f*(1-ART_eff*ART_m_lv)*(1-prep_eff*prep_f)*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg) 
Pm2=(1-beta_m*(1-ART_eff*ART_f_lv)*(1-prep_eff*prep_m)*(1-cond_m2*cond_eff)*(1-circ*circ_eff))^(acts_m2_main+acts_m2_reg) 
Pf1=(1-beta_f*(1-ART_eff*ART_m_lv)*(1-prep_eff*prep_f)*(1-cond_f1*cond_eff))^acts_f1_main 
Pm1=(1-beta_m*(1-ART_eff*ART_f_lv)*(1-prep_eff*prep_m)*(1-cond_m1*cond_eff)*(1-circ*circ_eff))^acts_m1_main 
Pf3=(1-beta_f*(1-ART_eff*ART_m_lv)*(1-prep_eff*prep_f)*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg) 
Pm3=(1-beta_m*(1-ART_eff*ART_f_lv)*(1-prep_eff*prep_m)*(1-cond_m3*cond_eff)*(1-circ*circ_eff))^(acts_m3_main+acts_m3_casual+acts_m3_reg) 
#partner in pre-aids stage
Xf2=(1-high_v2*beta_f*(1-ART_eff*ART_m_pa)*(1-prep_eff*prep_f)*(1-cond_f2*cond_eff))^(acts_f2_main+acts_f2_reg) 
Xm2=(1-high_v2*beta_m*(1-ART_eff*ART_f_pa)*(1-prep_eff*prep_m)*(1-cond_m2*cond_eff)*(1-circ*circ_eff))^(acts_m2_main+acts_m2_reg) 
Xf1=(1-high_v2*beta_f*(1-ART_eff*ART_m_pa)*(1-prep_eff*prep_f)*(1-cond_f1*cond_eff))^acts_f1_main 
Xm1=(1-high_v2*beta_m*(1-ART_eff*ART_f_pa)*(1-prep_eff*prep_m)*(1-cond_m1*cond_eff)*(1-circ*circ_eff))^acts_m1_main 
Xf3=(1-high_v2*beta_f*(1-ART_eff*ART_m_pa)*(1-prep_eff*prep_f)*(1-cond_f3*cond_eff))^(acts_f3_main+acts_f3_casual+acts_f3_reg) 
Xm3=(1-high_v2*beta_m*(1-ART_eff*ART_f_pa)*(1-prep_eff*prep_m)*(1-cond_m3*cond_eff)*(1-circ*circ_eff))^(acts_m3_main+acts_f3_casual+acts_m3_reg) 

## Probability that infected partner transmits HIV to the susceptible partner
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

## Probability of infection = prob of acquiring infection upon exposure * prob that infected partner transmits infection
# female
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

## FORCE OF INFECTION
pi_f1= 1 - ((1-phi_f1_1)^(df1*rho_f1_1)*(1-phi_f1_2)^(df1*rho_f1_2)*(1-phi_f1_3)^(df1*rho_f1_3)) 
pi_f2= 1 - ((1-phi_f2_1)^(df2*rho_f2_1)*(1-phi_f2_2)^(df2*rho_f2_2)*(1-phi_f2_3)^(df2*rho_f2_3)) 
pi_f3= 1 - ((1-phi_f3_1)^(df3*rho_f3_1)*(1-phi_f3_2)^(df3*rho_f3_2)*(1-phi_f3_3)^(df3*rho_f3_3)) 

pi_m1= 1 - ((1-phi_m1_1)^(dm1*rho_m1_1)*(1-phi_m1_2)^(dm1*rho_m1_2)*(1-phi_m1_3)^(dm1*rho_m1_3)) 
pi_m2= 1 - ((1-phi_m2_1)^(dm2*rho_m2_1)*(1-phi_m2_2)^(dm2*rho_m2_2)*(1-phi_m2_3)^(dm2*rho_m2_3)) 
pi_m3= 1 - ((1-phi_m3_1)^(dm3*rho_m3_1)*(1-phi_m3_2)^(dm3*rho_m3_2)*(1-phi_m3_3)^(dm3*rho_m3_3)) 


### ODE PARAMETERS
state=c(Sf1=Sf1,Hf1=Hf1,Yf1=Yf1,Zf1=Zf1,Af1=Af1,
        Sf2=Sf2,Hf2=Hf2,Yf2=Yf2,Zf2=Zf2,Af2=Af2,
        Sf3=Sf3,Hf3=Hf3,Yf3=Yf3,Zf3=Zf3,Af3=Af3,
        Sm1=Sm1,Hm1=Hm1,Ym1=Ym1,Zm1=Zm1,Am1=Am1,
        Sm2=Sm2,Hm2=Hm2,Ym2=Ym2,Zm2=Zm2,Am2=Am2,
        Sm3=Sm3,Hm3=Hm3,Ym3=Ym3,Zm3=Zm3,Am3=Am3)

times=seq(0.01,180,by=0.01) # Interventions run for 15 years

parameters=c(nu=1/nu, lambda=lambda, psi=psi,omega=1/omega, epi_m1=1/epi_m1, epi_m2=1/epi_m2, epi_m3=1/epi_m3, epi_f1=1/epi_f1, epi_f2=1/epi_f2, epi_f3=1/epi_f3,
             delta=1/delta,pi_f1=pi_f1,pi_m1=pi_m1,pi_f2=pi_f2,pi_m2=pi_m2,pi_f3=pi_f3,pi_m3=pi_m3)

require(deSolve)
trial1_FM_int=as.data.frame(ode(y=state,times=times,func=SHYZA_FM_INT,parms=parameters))

### MODEL OUTPUT ###

N_end_f1_int=sum(trial1_FM_int[18000,2:6])
N_end_f2_int=sum(trial1_FM_int[18000,7:11])
N_end_f3_int=sum(trial1_FM_int[18000,12:16])
N_end_m1_int=sum(trial1_FM_int[18000,17:21])
N_end_m2_int=sum(trial1_FM_int[18000,22:26])
N_end_m3_int=sum(trial1_FM_int[18000,27:31])

# Number of cases
HIV_f1_int=(trial1_FM_int[18000,6]+trial1_FM_int[18000,3]+trial1_FM_int[18000,4]+trial1_FM_int[18000,5]) 
HIV_f2_int=(trial1_FM_int[18000,11]+trial1_FM_int[18000,8]+trial1_FM_int[18000,9]+trial1_FM_int[18000,10]) 
HIV_f3_int=(trial1_FM_int[18000,16]+trial1_FM_int[18000,13]+trial1_FM_int[18000,14]+trial1_FM_int[18000,15]) 
HIV_m1_int=(trial1_FM_int[18000,21]+trial1_FM_int[18000,18]+trial1_FM_int[18000,19]+trial1_FM_int[18000,20]) 
HIV_m2_int=(trial1_FM_int[18000,26]+trial1_FM_int[18000,23]+trial1_FM_int[18000,24]+trial1_FM_int[18000,25]) 
HIV_m3_int=(trial1_FM_int[18000,31]+trial1_FM_int[18000,28]+trial1_FM_int[18000,29]+trial1_FM_int[18000,30])

HIV_all_int=(sum(trial1_FM_int[18000,3:6])+sum(trial1_FM_int[18000,8:11])+sum(trial1_FM_int[18000,13:16])+sum(trial1_FM_int[18000,18:21])
             +sum(trial1_FM_int[18000,23:26])
             +sum(trial1_FM_int[18000,28:31]))

# Prevalence in each class
prev_f1_int=(trial1_FM_int[18000,6]+trial1_FM_int[18000,3]+trial1_FM_int[18000,4]+trial1_FM_int[18000,5])/N_end_f1_int 
prev_f2_int=(trial1_FM_int[18000,11]+trial1_FM_int[18000,8]+trial1_FM_int[18000,9]+trial1_FM_int[18000,10])/N_end_f2_int 
prev_f3_int=(trial1_FM_int[18000,16]+trial1_FM_int[18000,13]+trial1_FM_int[18000,14]+trial1_FM_int[18000,15])/N_end_f3_int 
prev_m1_int=(trial1_FM_int[18000,21]+trial1_FM_int[18000,18]+trial1_FM_int[18000,19]+trial1_FM_int[18000,20])/N_end_m1_int 
prev_m2_int=(trial1_FM_int[18000,26]+trial1_FM_int[18000,23]+trial1_FM_int[18000,24]+trial1_FM_int[18000,25])/N_end_m2_int 
prev_m3_int=(trial1_FM_int[18000,31]+trial1_FM_int[18000,28]+trial1_FM_int[18000,29]+trial1_FM_int[18000,30])/N_end_m3_int 

# HIV prevalence in the total population
N_total_end_int=sum(trial1_FM_int[18000,2:31])
HIV_prev_int=(sum(trial1_FM_int[18000,3:6])+sum(trial1_FM_int[18000,8:11])+sum(trial1_FM_int[18000,13:16])+sum(trial1_FM_int[18000,18:21])
              +sum(trial1_FM_int[18000,23:26])
              +sum(trial1_FM_int[18000,28:31]))/N_total_end_int
prev_all_int=round(HIV_prev_int*100,2)

# Relative reduction in prevalence compared to baseline
RR_f1 = (prev_f1_base - prev_f1_int) / prev_f1_base
RR_f2 = (prev_f2_base - prev_f2_int) / prev_f2_base
RR_f3 = (prev_f3_base - prev_f3_int) / prev_f3_base
RR_m1 = (prev_m1_base - prev_m1_int) / prev_m1_base
RR_m2 = (prev_m2_base - prev_m2_int) / prev_m2_base
RR_m3 = (prev_m3_base - prev_m3_int) / prev_m3_base
RR_all = (HIV_prev_base - HIV_prev_int) / HIV_prev_base

# Relative reduction in prevalence since start of interventions
RR_f1_start = (prev_f1 - prev_f1_int) / prev_f1
RR_f2_start = (prev_f2 - prev_f2_int) / prev_f2
RR_f3_start = (prev_f3 - prev_f3_int) / prev_f3
RR_m1_start = (prev_m1 - prev_m1_int) / prev_m1
RR_m2_start = (prev_m2 - prev_m2_int) / prev_m2
RR_m3_start = (prev_m3 - prev_m3_int) / prev_m3
RR_all_start = (HIV_prev - HIV_prev_int) / HIV_prev

# Difference in number of infections
cum_f1 = HIV_f1_base - HIV_f1_int
cum_f2 = HIV_f2_base - HIV_f2_int
cum_f3 = HIV_f3_base - HIV_f3_int
cum_m1 = HIV_m1_base - HIV_m1_int
cum_m2 = HIV_m2_base - HIV_m2_int
cum_m3 = HIV_m3_base - HIV_m3_int
cum_all = HIV_all_base - HIV_all_int


### WRITE OUTPUT FOR EACH SCENARIO

out=data.frame(sim=sim,N_total=round(N_total_end,0), prev_all=round(prev_all,2), prev_f1=round(prev_f1*100,2), prev_f2=round(prev_f2*100,2), 
               prev_f3=round(prev_f3*100,2), prev_m1=round(prev_m1*100,2), prev_m2=round(prev_m2*100,2), prev_m3=round(prev_m3*100,2),
               N_total_int=round(N_total_end_int,0), prev_all_int=round(prev_all_int,2), prev_f1_int=round(prev_f1_int*100,2), prev_f2_int=round(prev_f2_int*100,2), 
               prev_f3_int=round(prev_f3_int*100,2), prev_m1_int=round(prev_m1_int*100,2), prev_m2_int=round(prev_m2_int*100,2), prev_m3_int=round(prev_m3_int*100,2),
               N_total_base=round(N_total_end_base,0), prev_all_base=round(prev_all_base,2), prev_f1_base=round(prev_f1_base*100,2), prev_f2_base=round(prev_f2_base*100,2), 
               prev_f3_base=round(prev_f3_base*100,2), prev_m1_base=round(prev_m1_base*100,2), prev_m2_base=round(prev_m2_base*100,2), prev_m3_base=round(prev_m3_base*100,2),
               RR_f1=RR_f1, RR_f2=RR_f2, RR_f3=RR_f3, RR_m1=RR_m1, RR_m2=RR_m2, RR_m3=RR_m3, RR_all=RR_all,
               RR_f1_start=RR_f1_start, RR_f2_start=RR_f2_start, RR_f3_start=RR_f3_start, RR_m1_start=RR_m1_start, RR_m2_start=RR_m2_start, RR_m3_start=RR_m3_start, RR_all_start=RR_all_start,
               cum_f1=cum_f1, cum_f2=cum_f2, cum_f3=cum_f3, cum_m1=cum_m1, cum_m2=cum_m2, cum_m3=cum_m3, cum_all=cum_all)

write.table(out,file=paste0("output_sim_", sim, ".csv"),append=F,col.names=FALSE,row.names=FALSE,sep=",")
