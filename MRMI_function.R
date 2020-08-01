fMRMI<-function(dat,w,K,L){
  y0<-dat[,"y0"]
  y1<-dat[,"y1"]
  r<-dat[,"r"]
  ry1<-y1[r==1]
  ry0=y0[r==0]
  
  n<-dim(dat)[1]
  nry1<-length(ry1)
  nry0=length(ry0)
  nm_y1<-n-nry1
  nm_y0=n-nry0
  id<-1:n
  s0_y1<-fMRI1(dat)
  s0_y0=fMRI0(dat)
  THETA<-NULL
  iter<-0
  repeat{
    iter<-iter+1
    ids<-sample(id,n,replace=T)
    dats<-dat[ids,]
    s_y1<-fMRI1(dats)
    s_y0<-fMRI0(dats)
    rs<-dats[,"r"]
    y1s<-dats[,"y1"]
    y0s=dats[,"y0"]
    ry1s<-y1s[rs==1]
    ry0s=y0s[rs==0]
    
    Etheta<-NULL
    EVAR<-NULL
    for(iter2 in 1:9){
      em_0_y1<-s0_y1[,(2*iter2-1)]
      ep_0_y1<-s0_y1[,(2*iter2)]
      m_em_y1<-em_0_y1[r==0]
      m_ep_y1<-ep_0_y1[r==0]
      em_s_y1<-s_y1[,(2*iter2-1)]
      ep_s_y1<-s_y1[,(2*iter2)]
      r_em_y1<-em_s_y1[rs==1]
      r_ep_y1<-ep_s_y1[rs==1]
      nrs_y1<-length(r_ep_y1)
      
      OM_a1_y1<-kronecker(m_em_y1,(-r_em_y1),FUN = "+")
      OM_a2_y1<-matrix(OM_a1_y1,nm_y1,nrs_y1,byrow=T)
      PM_a1_y1<-kronecker(m_ep_y1,(-r_ep_y1),FUN = "+")
      PM_a2_y1<-matrix(PM_a1_y1,nm_y1,nrs_y1,byrow=T)
      Dist_y1<-sqrt(OM_a2_y1^{2}*w+PM_a2_y1^{2}*(1-w))
      sidDist_y1<-t(apply(Dist_y1,1,order)) 
      sidDist_K_y1<-sidDist_y1[,1:K] # Error here
      idmy_y1<-apply(sidDist_K_y1,1,fsamp)
      
      imp_y1<-ry1s[idmy_y1]
      
      em_0_y0<-s0_y0[,(2*iter2-1)]
      ep_0_y0<-s0_y0[,(2*iter2)]
      m_em_y0<-em_0_y0[r==1]
      m_ep_y0<-ep_0_y0[r==1]
      em_s_y0<-s_y0[,(2*iter2-1)]
      ep_s_y0<-s_y0[,(2*iter2)]
      r_em_y0<-em_s_y0[rs==0]
      r_ep_y0<-ep_s_y0[rs==0]
      nrs_y0<-length(r_ep_y0)
      
      OM_a1_y0<-kronecker(m_em_y0,(-r_em_y0),FUN = "+")
      OM_a2_y0<-matrix(OM_a1_y0,nm_y0,nrs_y0,byrow=T)
      PM_a1_y0<-kronecker(m_ep_y0,(-r_ep_y0),FUN = "+")
      PM_a2_y0<-matrix(PM_a1_y0,nm_y0,nrs_y0,byrow=T)
      Dist_y0<-sqrt(OM_a2_y0^{2}*w+PM_a2_y0^{2}*(1-w))
      sidDist_y0<-t(apply(Dist_y0,1,order))
      sidDist_K_y0<-sidDist_y0[,1:K]     # Error coming from here, figure out
      idmy_y0<-apply(sidDist_K_y0,1,fsamp)
      
      imp_y0<-ry0s[idmy_y0]
      
      part_one=ry1-imp_y0
      part_two=imp_y1-ry0
      final=append(part_one,part_two)
      
      etheta<-(sum(final))/n
      Etheta<-c(Etheta,etheta)
      
      eVAR<-var(final)/n
      EVAR<-c(EVAR,eVAR)
    }
    
    THETAs<-c(Etheta,EVAR)
    THETA<-rbind(THETA,THETAs)
    if(iter==L){break}
  }
  THETA1<-THETA[,1:9]
  THETA2<-THETA[,10:18]
  ETHETA<-apply(THETA1,2,mean)
  Bm<-apply(THETA1,2,var)
  Um<-apply(THETA2,2,mean)
  EV<-(1+1/L)*Bm+Um
  StdError=sqrt(EV)
  LB<-ETHETA-1.96*sqrt(EV)
  UB<-ETHETA+1.96*sqrt(EV)
  ind_CI<-as.logical((LB<0 & UB<0)|(LB>0 & UB >0))
  len_CI<-UB-LB
  RESF<-list(means=ETHETA,SEs=StdError,Low=LB,High=UB,length=len_CI,sig=ind_CI)
  return(RESF)
}

fMRMI_final=function(dat,w=.5,K=3,L=3){return(fMRMI(dat,w,K,L))}