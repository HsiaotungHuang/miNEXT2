#' Mixture diversity of 2 communities
#' \code{miNEXT} computes the mixture diversity of 3 communities.
#' @import dplyr 
#' @import iNEXT
#' @importFrom chaoUtility Boot_p
#' @param data a matrix/data.frame.
#' @param knots an integer specifying the number of points of m1. Default to be 15.
#' @param m1 a vector specifying the values of m1 where the mixture diversity will be computed.
#' @param nboot an integer specifying the number of bootstrap times to build confidence interval. 
#' Use 0 to skip bootstrap.
#' @return a list of 2 components: \code{$Each} a table of diversity of single community; \code{$Mixture} a 
#' table of mixture diversity.
#' @export
#' 
#' 



Abundance<-function(data=NULL,knots=10,size=NULL,nboots=0,Sdecompose=T ){
  data1mn<-Create_data1mn(data, knots, size)
  data1<-data1mn$data1
  m<-data1mn$m
  n<-data1mn$n
  
  result_abun<-NULL
  result_abun_CI<-NULL
  result_abun<-Abun(data=data1,knots=knots,size=size,Sdecompose=Sdecompose)
  esti<-result_abun
  if(nboots>1){
    boot_sample<-Abun_CreatBootstrapSample(data=data1,nboots =nboots)
    out_boot <- sapply(1:nboots,function(k){
      #  print(paste("boots",k))
      bdata<-boot_sample[[k]]
      bdata <- bdata[rowSums(bdata)>0,]
      out<-NULL
      out =  Abun(data=bdata,knots=knots,size=size,Sdecompose=Sdecompose)
      out012 <- matrix(c(out$q0[,5], out$q1[,5], out$q2[,5]),
                       ncol = 3, dimnames = list(NULL,c("q0","q1","q2")))
      out_p <- matrix(c(out$q0_ana[,3], out$q0_ana[,4], out$q0_ana[,5]),
                      ncol = 3, dimnames = list(NULL,colnames(out$q0_ana)[3:5])) %>% rbind(., .[1:(nrow(out$q0)-nrow(out$q0_ana)),])
      cbind(out012, out_p)
      
    }, simplify = "array")
    
    sd_boot <- apply(out_boot,MARGIN = c(1,2), sd) 
    esti$q0 <- esti$q0 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot[,1]) , UCL = (.[,5] + 1.96*sd_boot[,1]), s.e. = sd_boot[,1])
    esti$q1 <- esti$q1 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot[,2]) , UCL = (.[,5] + 1.96*sd_boot[,2]), s.e. = sd_boot[,2])
    esti$q2 <- esti$q2 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot[,3]) , UCL = (.[,5] + 1.96*sd_boot[,3]), s.e. = sd_boot[,3])
    tmp <- nrow(esti$q0_ana)
    esti$q0_ana <- esti$q0_ana %>% cbind(., LCL = (.[1:tmp,3:5] - 1.96*sd_boot[1:tmp,4:6]) , UCL = (.[1:tmp,3:5] + 1.96*sd_boot[1:tmp,4:6]),
                                         s.e. = sd_boot[1:tmp,4:6])
    esti$q0_ana[ esti$q0_ana < 0 ] = 0
    
  
    sd012<-c(unique(sd_boot[,1]),unique(sd_boot[,2]),unique(sd_boot[,3]))
    esti$q012 <- esti$q012 %>% mutate(LCL = Mixture-qnorm(0.975)*sd012,UCL = Mixture+qnorm(0.975)*sd012,
                              se=sd012)
    
  }
  return(esti)
}
    
    

#' Mixture diversity of 2 communities
#' @param data1 a matrix/data.frame.
#' @param m a matrix specifying the values of m1,m2 where the mixture diversity will be computed.
#' @param n a matrix specifying the total observed values : n1,n2.
#' @return a table: \code{$Mix2} a table of mixture diversity for q=0,1,2.

Abun_Mix <- function(data, knots = 10, size = NULL ){
  data1mn<-Create_data1mn(data, knots, size)
  data1<-data1mn$data1
  m<-data1mn$m
  n<-data1mn$n
  
  x1<-data1[,1]
  x2<-data1[,2]
  
  #compute the mixture diversity
  #make sure if needed extrapolation?
  mid_inext<- apply(m, 1, function(z) sum(z<=n))
  datap <- NULL
  mid_in<-NULL
  mid_ext<-NULL
  
  if( sum( mid_inext==2) < nrow(m)){
    datap <- Boot_p(x = data1,Bootype = "JADE",datatype = "abundance") %>% 
      sapply(., function(x){x[1:nrow(data1)]})
    #zzdata1<<-data1
    #zzdatap<<-datap
    mid_in<-which(mid_inext %in% 2)    
    mid_ext<-which(mid_inext %in% 1) 
    m_in<-m[mid_in,]
    m_ext<-m[mid_ext,]
    #q=1
    q1_in<-Dq01_in(x1,x2,m_in,q=1)
    q1_ext<-Dq1_ext(x1,x2,m_ext,n,datap)
    q1_mix<-c(q1_in,q1_ext)
    #q=0
    q0_in<-Dq01_in(x1,x2,m_in,q=0)
    q0_ext<-Dq0_ext(x1,x2,m_ext,n,datap)
    q0_mix<-c(q0_in,q0_ext)
    
  }else{
    datap <- NULL
    q1_mix<-Dq01_in(x1,x2,mid_inext,q=1)
    q0_mix<-Dq01_in(x1,x2,mid_inext,q=0)
  }
  
  #q=2 in and ext are the same formula
  q2_mix<-Dq2(x1,x2,m,n1,n2)
  
  m1<-m[,1]
  m2<-m[,2]
  Mix2<-tibble(m1=rep(m1,3),m2=rep(m2,3),q=rep(c(0,1,2),each=nrow(m)),Mixture =c(q0_mix,q1_mix,q2_mix))
  return(Mix2)
}


Abun <- function(data, knots = 10, size = NULL,Sdecompose=T){
  data1mn<-Create_data1mn(data, knots, size)
  data1<-data1mn$data1
  m<-data1mn$m
  n<-data1mn$n
  
  x1 = data1[, 1]
  x2 = data1[, 2]
  m1<-m[,1]
  m2<-m[,2]
  n1<-n[1]
  n2<-n[2]
  
  #compute diveristy of each single community
  Each<- lapply(1:ncol(data1),function(i){
    #mtmp <-  seq(0,max(n[i],n[1]),length.out = knots) %>% round(0)
    mtmp <-  m[,i]
    cname<-names(n)[i]
    output <- iNEXT(x = data1[,i],q = c(0,1,2),datatype = "abundance", size = mtmp,se = F)$iNextEst
    output$qD[is.na(output$qD)] <- 0
    output %>% as_tibble() %>% select(m, q = order, !!cname := qD) 
  }) 
  
  
  
  #compute the mixture diversity
  Mix2<-Abun_Mix(data1, knots,size)
  
  #merge each and mixture into one dataframe
  each1<-Each[[1]] %>% mutate(m1=m)
  each2<-Each[[2]] %>% mutate(m2=m)
  Result_q012<-Mix2 %>% left_join(each1,by=c("m1","q")) %>% select(-m)
  Result_q012<-Result_q012 %>% left_join(each2,by=c("m2","q")) %>% select(-m)
  Result_q012<- Result_q012 %>% select(m1,m2,q,5:6,Mixture) 
  
  Result_q012.1<-Result_q012
  tmp<-Result_q012.1 %>% filter(m2==n2)
  Result_q012.1<-rbind(Result_q012.1,tmp)
  Result_q012.1<-Result_q012.1 %>% arrange(q,-m1,m2)
  
  output0<-Result_q012.1 %>% filter(q==0) %>% select(-q) %>% data.matrix()
  output1<-Result_q012.1 %>% filter(q==1) %>% select(-q) %>% data.matrix()
  output2<-Result_q012.1 %>% filter(q==2) %>% select(-q) %>% data.matrix()
  
 
  
  #Species decompose
  if(Sdecompose==T){
    q0_ana<-SDecompose(data,knots, size,q=0)
  }
  
  
  
  list(q0 = output0, q1 = output1, q2 = output2, q0_ana = q0_ana,q012 =Result_q012)
  
}
