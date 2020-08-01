##########################Important functions:

#####1.One dimensional Newton method: For solving root

newton=function(fun,derf,x0,eps) # Traditional Newton's method, just provide
  # a fun, its derivative, an initial guess x0 and a tolerance eps (to tell us when
  # to terminate)
{
           iter=0
           repeat{
                  iter=iter+1
                  x1=x0-fun(x0)/derf(x0)
                  if(abs(x0-x1)<eps||abs(fun(x1))<.1e-10)
                      break
                  x0=x1
                  cat("****** Iter.No:",iter," Current Iterate=",x1,fill=T)
              }
            return(x1)
}



#####2.Two dimensional Newton method:

newton2=function(fun,derf,x0,eps)
{
           iter=0
           repeat{
                  iter=iter+1
                  x1=x0-solve(derf(x0))%*%fun(x0)
                  if(sqrt(sum((x0-x1)^{2}))<eps||sqrt(sum(fun(x1)^{2}))<.1e-10)
                      break
                  x0=x1
                  cat("****** Iter.No:",iter," Current Iterate=",x1,fill=T)
              }
            return(x1)
}


#####3.Lambda function for the scalar:(Which caculate the Langrange multiplier in EL procedure)

Lfs<-function(u){
lambda0<-0
k<-0
gamma0<-1
epsilon<-10^{-8}
nu<-length(u)

repeat{
       D1<-sum(u/(1+lambda0*u))/nu
       D2<-D1/(-sum(u^{2}/(1+lambda0*u)^{2})/nu)
       if(abs(D2)<epsilon){M<-lambda0
                           break
                          }
       repeat{ 
              deta0<-gamma0*D2
              if(sum(as.numeric(1+(lambda0-deta0)*u<=0))>0){ gamma0<-gamma0/2
                                                           }
              else{break}
             }
       lambda0<-lambda0-deta0
       k<-k+1
       gamma0<-(k+1)^{-1/2}
       if(k>=100){ M<-0
                   break
                 }
       }
return(M)
}

#####4.Lambda function for vector:

###u is a r*n matrix 

Lfv<-function(u){

lambda0<-rep(0,dim(u)[1])
k<-0
gamma0<-1
epsilon<-10^{-8}
nu<-dim(u)[2]

repeat{
       a<-t(1+t(lambda0)%*%u)
       a<-as.numeric(a)
       B<-NULL
       iter<-0
       repeat{
             iter<-iter+1
             B<-rbind(B,a)
             if(iter==dim(u)[1]){break}
             }
       unew<-u/B
       D1<-apply(unew,1,mean)
       D2<-solve(-nu^{-1}*unew%*%t(unew))%*%D1
       if(sqrt(sum(D2^{2}))<epsilon){M<-lambda0
                                     break
                                    }
       repeat{ 
              deta0<-gamma0*D2
              if(sum(as.numeric(t(1+t(lambda0-deta0)%*%u)<=0))>0){ gamma0<-gamma0/2
                                                                 }
              else{break}
             }
       lambda0<-lambda0-deta0
       k<-k+1
###    gamma0<-(k+1)^{-1/2}
       gamma0<-1
       if(k>=200){ M<-rep(0,dim(u)[1])
                   break
                 }
       }
M<-as.numeric(M)
return(M)
}

fsamp<-function(tt){
  ress<-sample(tt,1)
  return(ress)
}

######################################################
#  R function for solving g_1(lambda)=0 in Wu (2005) #
######################################################

#########################################
# R Function for EL Lagrange Multiplier #
#     x is of dimension 2 or higher     #
# Input: u=(x_1,x_2,...,x_n)            #
#        ds=(d_1,d_2,...,d_n)           #
#        (Design Wights: d_i=1/pi_i)    #
#        ds=(1,1,...1) for iid data     #
#        mu: benchmark means for x      #
# Output: lambda=M                      #
#                                       #
# Written by Changbao Wu, March, 2000   #
# Modified by Randy Sitter, June, 2006  #
#########################################

Lag2=function(u,ds,mu){ 
  n=length(ds)
  u=u-rep(1,n)%*%t(mu)
  M=0*mu
  dif=1
  tol=1e-8
  k=0
  while(dif>tol & k<=50){
  D1=t(u)%*%((ds/(1+u%*%M))*rep(1,n))
  DD=-t(u)%*%(c((ds/(1+u%*%M)^2))*u)
  D2=solve(DD,D1,tol=1e-40)
  dif=max(abs(D2))
  rule=1
  while(rule>0){
    rule=0
    if(min(1+t(M-D2)%*%t(u))<=0) rule=rule+1
    if(rule>0) D2=D2/2
               }
    M=M-D2
    k=k+1
                   }
  if(k>=50) M=0*mu
  return(as.vector(M))
          }

# Main fMRMI functino


























