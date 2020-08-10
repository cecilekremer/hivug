

#Read in input-output dataset
data = read.csv("data_ARF_v28.csv", sep=",", header=T)
attach(data); dim(data); summary(GOF)

#################################
### 2. ACTIVITY REGION FINDER ###
#################################

library(ARF)

## Create binary response variable: 1 if parameter set belongs to top 1%, 0 otherwise

data$top1[rank(data$GOF)<=1000]=1
data$top1[rank(data$GOF)>1000]=0
attach(data); top1=as.factor(top1)
table(top1)

#check GOF of top 1%
test=data[which(data$top1==1),]
summary(test$GOF)

## Generate an ARF object that contains the info for constructing a tree

attach(data)
# remove fixed parameters !
best.arf2=f.arf(top1 ~ lambda + psi + C_f + C_m + high_v2 + beta_f
			+ beta_m + cond_m1 + cond_f2 + cond_f3
			+ cond_eff + acts_m3_casual + acts_m2_reg
			+ acts_m3_reg + acts_f1_main
			+ acts_f2_main + acts_f3_main + acts_m2_main
			+ acts_m3_main + df1 + dm2 + dm3 + N
			+ A_main + A_regular + A_casual + nu, 
			data=data, const=5, minn=20)         


## Generate a report:
f.report(best.arf2, file="ARF_v28.pdf")







