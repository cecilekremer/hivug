setwd("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/corrected calibration/FINAL SIM/Interventions/uncertainty")

##################################
### 1. MERGE INPUT-OUTPUT DATA ###
##################################

out=read.csv("output_unct_final.csv", sep=",", header=F)
dim(out)
out2=read.csv("output_extra.csv", sep=",",header=F)
dim(out2)
out=rbind(out,out2)
colnames(out)<-c("sim","N_total","prev_all","prev_f1","prev_f2","prev_f3","prev_m1","prev_m2","prev_m3",
                 "N_total_int","prev_all_int","prev_f1_int","prev_f2_int","prev_f3_int","prev_m1_int","prev_m2_int","prev_m3_int",
                 "N_total_base","prev_all_base","prev_f1_base","prev_f2_base","prev_f3_base","prev_m1_base","prev_m2_base","prev_m3_base",
                 "RR_f1","RR_f2","RR_f3","RR_m1","RR_m2","RR_m3","RR_all",
                 "RR_f1_start","RR_f2_start","RR_f3_start","RR_m1_start","RR_m2_start","RR_m3_start","RR_all_start",
                 "cum_f1","cum_f2","cum_f3","cum_m1","cum_m2","cum_m3","cum_all")


input=read.csv("data_unct_hpc.csv", sep=",", header=T)
dim(input)
#merge input and output datasets
result = merge(out,input,by="sim"); dim(result); #View(result)
result_final = result

#############################
# load("WS_DATA.RData")
dim(result_final)
names(result_final)
attach(result_final)


################################################
### BEST FIT BASELINE                        ###
################################################

summary(result_final$RR_all); summary(result_final$cum_all)
quantile(result_final$RR_all,c(0.025,0.5,0.975))
quantile(result_final$cum_all,c(0.025,0.5,0.975))

summary(((prev_all_base/100*N_total_base) - (prev_all_int/100*N_total_int))/(prev_all_base/100*N_total_base))



# 95% CrI relative reduction prevalence, best fit 
quantile(RR_all,c(0.025,0.5,0.975))
quantile(cum_all,c(0.025,0.5,0.975))


## Intervention parameters: PrEP has highest correlation with outcome
library(ggplot2)
df=data.frame(x=circ,y=RR_all)
p1=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("VMMC") + ylab("Relative reduction in prevalence") +
  scale_fill_distiller(palette="Greys", direction=1)
df=data.frame(x=int_cond,y=RR_all)
p2=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("Increase in condom use") + ylab(" ") +
  scale_fill_distiller(palette="Greys", direction=1)
df=data.frame(x=prep_m,y=RR_all)
p3=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("HIV- men on PrEP") + ylab(" ") + 
  scale_fill_distiller(palette="Greys", direction=1)
df=data.frame(x=prep_f,y=RR_all)
p4=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("HIV- women on PrEP") + ylab(" ") + 
  scale_fill_distiller(palette="Greys", direction=1)
df=data.frame(x=ART_m,y=RR_all)
p5=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("HIV+ men on ART") + ylab("Relative reduction in prevalence") + 
  scale_fill_distiller(palette="Greys", direction=1)
df=data.frame(x=ART_f,y=RR_all)
p6=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("HIV+ women on ART") + ylab(" ") + 
  scale_fill_distiller(palette="Greys", direction=1)
df=data.frame(x=prop_hct_m,y=RR_all)
p7=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("Men accepting HCT") + ylab(" ") + 
  scale_fill_distiller(palette="Greys", direction=1)
df=data.frame(x=prop_hct_f,y=RR_all)
p8=ggplot(df,aes(x=x,y=y)) + geom_bin2d(show.legend=F) + xlab("Women accepting HCT") + ylab(" ") + 
  scale_fill_distiller(palette="Greys", direction=1)
library(gridExtra)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4)





