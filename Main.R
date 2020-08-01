source("helper_functions.R")
source("MRMI_function.R")

dat=read.csv("example_data.csv") # Replace with your own data


my_correct_reg=ry~rx1+rx2+rx3+rx4+rx5
my_bad_reg=ry~rx1*rx3+rx1*rx2+rx6+rx7+rx8

my_correct_prop=r~x1+x2+x3+x4+x5
my_bad_prop=r~x1*x3+x1*x2+x6+x7+x8

###### Modeling function for Mu_1 #####

fMRI1<-function(dat,correct_reg=my_correct_reg,bad_reg=my_bad_reg,
                correct_prop=my_correct_prop,
                bad_prop=my_bad_prop){
  
  y0<-dat[,"y0"]
  y1<-dat[,"y1"]
  r<-dat[,"r"]
  ry<-y1[r==1]
  
  n<-length(y1)
  nr<-length(ry)
  nm<-n-nr
  
  ################ Begin section to modify #########################
  x1<-dat[,"x1"]
  x2<-dat[,"x2"]
  x3<-dat[,"x3"]
  x4<-dat[,"x4"]
  x5<-dat[,"x5"]
  x6<-dat[,"x6"]
  x7<-dat[,"x7"]
  x8<-dat[,"x8"]
  # If you have more than 8 covariates, add more here, if you have fewer, remove
  # as many as needed
  
  rx1<-x1[r==1]
  rx2<-x2[r==1]
  rx3<-x3[r==1]
  rx4<-x4[r==1]
  rx5=x5[r==1]
  rx6<-x6[r==1]
  rx7<-x7[r==1]
  rx8=x8[r==1]
  # Add or remove from this list accordingly
  
  # Add or remove from these accordingly
  df_reg=data.frame(ry,rx1,rx2,rx3,rx4,rx5,rx6,rx7,rx8)
  df_prop=data.frame(r,x1,x2,x3,x4,x5,x6,x7,x8)
  df1=data.frame(rx1=dat$x1,rx2=dat$x2,rx3=dat$x3,rx4=dat$x4,rx5=dat$x5,rx6=dat$x6,
                        rx7=dat$x7,rx8=dat$x8)
  ##################### End section to modify #####################
  
  ######0111:
  l_O1<-lm(bad_reg,data=df_reg)
  em1=predict(l_O1,df1)
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  l_M2<-glm(bad_prop,family=binomial(logit),data=df_prop) #zmodel
  ep2<-l_M2$fitted.values
  
  lm_p<-lm(r~0+ep1+ep2)
  ebeta_p<-lm_p$coefficient
  ep_0111<-ebeta_p[1]*ep1+ebeta_p[2]*ep2
  em_0111<-em1
  
  ep_0111<-(ep_0111-mean(ep_0111))/sd(ep_0111)
  em_0111<-(em_0111-mean(em_0111))/sd(em_0111)
  
  ######1011:
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  l_M2<-glm(bad_prop,family=binomial(logit),data=df_prop) #zmodel
  ep2<-l_M2$fitted.values
  
  lm_p<-lm(r~0+ep1+ep2)
  ebeta_p<-lm_p$coefficient
  ep_1011<-ebeta_p[1]*ep1+ebeta_p[2]*ep2
  em_1011<-em1
  
  ep_1011<-(ep_1011-mean(ep_1011))/sd(ep_1011)
  em_1011<-(em_1011-mean(em_1011))/sd(em_1011)
  
  ######1101:
  
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  rem1<-em1[r==1]
  
  l_O2<-lm(bad_reg,data=df_reg) #zmodel
  em2=predict(l_O2,df1)
  rem2<-em2[r==1]

  l_M1<-glm(bad_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  
  lm_m<-lm(ry~0+rem1+rem2)
  ebeta_m<-lm_m$coefficient
  em_1101<-ebeta_m[1]*em1+ebeta_m[2]*em2
  ep_1101<-ep1
  
  ep_1101<-(ep_1101-mean(ep_1101))/sd(ep_1101)
  em_1101<-(em_1101-mean(em_1101))/sd(em_1101)
  
  ######1110:     # Left off here with replacements
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  rem1<-em1[r==1]
  
  l_O2<-lm(bad_reg,data=df_reg)
  em2=predict(l_O2,df1)
  rem2<-em2[r==1]
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  
  lm_m<-lm(ry~0+rem1+rem2)
  ebeta_m<-lm_m$coefficient
  em_1110<-ebeta_m[1]*em1+ebeta_m[2]*em2
  ep_1110<-ep1
  
  ep_1110<-(ep_1110-mean(ep_1110))/sd(ep_1110)
  em_1110<-(em_1110-mean(em_1110))/sd(em_1110)
  
  ######1111:
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  
  rem1<-em1[r==1]
  l_O2<-lm(bad_reg,data=df_reg)
  em2=predict(l_O2,df1)
  
  rem2<-em2[r==1]
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  l_M2<-glm(bad_prop,family=binomial(logit),data=df_prop)
  ep2<-l_M2$fitted.values
  
  lm_m<-lm(ry~0+rem1+rem2)
  ebeta_m<-lm_m$coefficient
  em_1111<-ebeta_m[1]*em1+ebeta_m[2]*em2
  lm_p<-lm(r~0+ep1+ep2)
  ebeta_p<-lm_p$coefficient
  ep_1111<-ebeta_p[1]*ep1+ebeta_p[2]*ep2
  
  ep_1111<-(ep_1111-mean(ep_1111))/sd(ep_1111)
  em_1111<-(em_1111-mean(em_1111))/sd(em_1111)
  
  
  #########1010
  
  l_O1<-lm(correct_reg,data=df_reg) #reg
  em_1010=predict(l_O1,df1)
 
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop) #prop
  ep_1010<-l_M1$fitted.values
  
  ep_1010<-(ep_1010-mean(ep_1010))/sd(ep_1010)
  em_1010<-(em_1010-mean(em_1010))/sd(em_1010)
  
  #########0101
  
  l_O1<-lm(bad_reg,data=df_reg) #reg
  em_0101=predict(l_O1,df1)
  
  l_M1<-glm(bad_prop,family=binomial(logit),data=df_prop) #prop
  ep_0101<-l_M1$fitted.values
  
  ep_0101<-(ep_0101-mean(ep_0101))/sd(ep_0101)
  em_0101<-(em_0101-mean(em_0101))/sd(em_0101)
  
  #########1001
  
  l_O1<-lm(correct_reg,data=df_reg) #reg
  em_1001=predict(l_O1,df1)
  
  l_M1<-glm(bad_prop,family=binomial(logit),data=df_prop) #prop
  ep_1001<-l_M1$fitted.values
  
  ep_1001<-(ep_1001-mean(ep_1001))/sd(ep_1001)
  em_1001<-(em_1001-mean(em_1001))/sd(em_1001)
  
  #########0110

  l_O1<-lm(bad_reg,data=df_reg) #reg
  em_0110=predict(l_O1,df1)
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop) #prop
  ep_0110<-l_M1$fitted.values
  
  ep_0110<-(ep_0110-mean(ep_0110))/sd(ep_0110)
  em_0110<-(em_0110-mean(em_0110))/sd(em_0110)
  
  
  etheta<-cbind(em_0111,ep_0111,em_1011,ep_1011,em_1101,ep_1101,em_1110,ep_1110,em_1111,ep_1111,em_1010,ep_1010,em_0101,ep_0101,
                em_1001,ep_1001,em_0110,ep_0110)
  return(etheta)
}

############ Modeling function for Mu_0 ######################

fMRI0<-function(dat,correct_reg=my_correct_reg,bad_reg=my_bad_reg,
                correct_prop=my_correct_prop,
                bad_prop=my_bad_prop){
  
  y0<-dat[,"y0"]
  y1<-dat[,"y1"]
  r<-dat[,"r"]
  ry<-y0[r==0]
  n<-length(y0)
  nr<-length(ry)
  nm<-n-nr
  
  ################## Begin section to modify ####################
  
  x1<-dat[,"x1"]
  x2<-dat[,"x2"]
  x3<-dat[,"x3"]
  x4<-dat[,"x4"]
  x5<-dat[,"x5"]
  x6<-dat[,"x6"]
  x7<-dat[,"x7"]
  x8<-dat[,"x8"]

 
  rx1<-x1[r==0]
  rx2<-x2[r==0]
  rx3<-x3[r==0]
  rx4<-x4[r==0]
  rx5=x5[r==0]
  rx6<-x6[r==0]
  rx7<-x7[r==0]
  rx8=x8[r==0]
  

  df_reg=data.frame(ry,rx1,rx2,rx3,rx4,rx5,rx6,rx7,rx8)
  df_prop=data.frame(r,x1,x2,x3,x4,x5,x6,x7,x8)

  df1=data.frame(rx1=dat$x1,rx2=dat$x2,rx3=dat$x3,rx4=dat$x4,rx5=dat$x5,rx6=dat$x6,
                 rx7=dat$x7,rx8=dat$x8)
  
  ###############End section to modify ########################
  
  ######0111:
  l_O1<-lm(bad_reg,data=df_reg)
  em1=predict(l_O1,df1)
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  l_M2<-glm(bad_prop,family=binomial(logit),data=df_prop)
  ep2<-l_M2$fitted.values
  
  lm_p<-lm(r~0+ep1+ep2)
  ebeta_p<-lm_p$coefficient
  ep_0111<-1-(ebeta_p[1]*ep1+ebeta_p[2]*ep2)
  em_0111<-em1
  
  ep_0111<-(ep_0111-mean(ep_0111))/sd(ep_0111)
  em_0111<-(em_0111-mean(em_0111))/sd(em_0111)
  
  ######1011:
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  l_M2<-glm(bad_prop,family=binomial(logit),data=df_prop)
  ep2<-l_M2$fitted.values
  
  lm_p<-lm(r~0+ep1+ep2)
  ebeta_p<-lm_p$coefficient
  ep_1011<-1-(ebeta_p[1]*ep1+ebeta_p[2]*ep2)
  em_1011<-em1
  
  ep_1011<-(ep_1011-mean(ep_1011))/sd(ep_1011)
  em_1011<-(em_1011-mean(em_1011))/sd(em_1011)
  
  ######1101:
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  rem1<-em1[r==0]
  
  l_O2<-lm(bad_reg,data=df_reg)
  em2=predict(l_O2,df1)
  rem2<-em2[r==0]
  l_M1<-glm(bad_prop,family=binomial(logit),data=df_prop)
  ep1<-1-l_M1$fitted.values # ADDED A 1-here 
  
  lm_m<-lm(ry~0+rem1+rem2)
  
  ebeta_m<-lm_m$coefficient
  em_1101<-ebeta_m[1]*em1+ebeta_m[2]*em2
  ep_1101<-ep1
  
  ep_1101<-(ep_1101-mean(ep_1101))/sd(ep_1101)
  em_1101<-(em1-mean(em1))/sd(em1)
  
  ######1110:
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  rem1<-em1[r==0]
  
  l_O2<-lm(bad_reg,data=df_reg)
  em2=predict(l_O2,df1)
  rem2<-em2[r==0]
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-1-l_M1$fitted.values # ADDED A 1- HERE
  
  lm_m<-lm(ry~0+rem1+rem2)
  ebeta_m<-lm_m$coefficient
  em_1110<-ebeta_m[1]*em1+ebeta_m[2]*em2
  ep_1110<-ep1
  
  ep_1110<-(ep_1110-mean(ep_1110))/sd(ep_1110)
  em_1110<-(em_1110-mean(em_1110))/sd(em_1110)
  
  ######1111:
  l_O1<-lm(correct_reg,data=df_reg)
  em1=predict(l_O1,df1)
  rem1<-em1[r==0]
  
  l_O2<-lm(bad_reg,data=df_reg)
  em2=predict(l_O2,df1)
  rem2<-em2[r==0]
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop)
  ep1<-l_M1$fitted.values
  l_M2<-glm(bad_prop,family=binomial(logit),data=df_prop)
  ep2<-l_M2$fitted.values
  
  lm_m<-lm(ry~0+rem1+rem2)
  ebeta_m<-lm_m$coefficient
  em_1111<-ebeta_m[1]*em1+ebeta_m[2]*em2
  lm_p<-lm(r~0+ep1+ep2)
  ebeta_p<-lm_p$coefficient
  ep_1111<-1-(ebeta_p[1]*ep1+ebeta_p[2]*ep2)
  
  ep_1111<-(ep_1111-mean(ep_1111))/sd(ep_1111)
  em_1111<-(em_1111-mean(em_1111))/sd(em_1111)
  
  
  #########1010
 
  l_O1<-lm(correct_reg,data=df_reg)
  em_1010=predict(l_O1,df1)
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop) #prop
  ep_1010<-1-l_M1$fitted.values
  
  ep_1010<-(ep_1010-mean(ep_1010))/sd(ep_1010)
  em_1010<-(em_1010-mean(em_1010))/sd(em_1010)
  
  #########0101

  l_O1<-lm(bad_reg,data=df_reg) #reg
  em_0101=predict(l_O1,df1)
  
  l_M1<-glm(bad_prop,family=binomial(logit),data=df_prop) #prop
  ep_0101<-1-l_M1$fitted.values
  
  ep_0101<-(ep_0101-mean(ep_0101))/sd(ep_0101)
  em_0101<-(em_0101-mean(em_0101))/sd(em_0101)
  
  #########1001

  l_O1<-lm(correct_reg,data=df_reg) #reg
  em_1001=predict(l_O1,df1)
  
  l_M1<-glm(bad_prop,family=binomial(logit),data=df_prop) #prop
  ep_1001<-1-l_M1$fitted.values
  
  ep_1001<-(ep_1001-mean(ep_1001))/sd(ep_1001)
  em_1001<-(em_1001-mean(em_1001))/sd(em_1001)
  
  #########0110

  l_O1<-lm(bad_reg,data=df_reg) #reg
  em_0110=predict(l_O1,df1)
  
  l_M1<-glm(correct_prop,family=binomial(logit),data=df_prop) #prop
  ep_0110<-1-l_M1$fitted.values
  
  ep_0110<-(ep_0110-mean(ep_0110))/sd(ep_0110)
  em_0110<-(em_0110-mean(em_0110))/sd(em_0110)
  
  etheta<-cbind(em_0111,ep_0111,em_1011,ep_1011,em_1101,ep_1101,em_1110,ep_1110,em_1111,ep_1111,em_1010,ep_1010,em_0101,ep_0101,
                em_1001,ep_1001,em_0110,ep_0110)
  return(etheta)
}

# Run the function to get your results

# Change these as desired:
# w=lambda
# K=number of units in nearest neighborhood
# L= number of imputations

Res_MRMI<-fMRMI_final(dat,w=.5,K=3,L=3)
Res_MRMI





