
#' Incidence(data, allpts = FALSE, size = NULL, knots = 20 ,nboots = 0) for incidence data, comupute composite diversity of any sample and species composition (shared and unique species)
#' along with their confidence interval.
#' @param data a Sx2 dataframe, the intact assemblage (main) assemblage should be the first column.
#' @param allpts specifying whether to compute all combinations of sampling units of two assemblages. Default is FALSE.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param knots the number of points that the mixture diveristy will be computed. Default is 20.
#' @param nboots the number of replication bootstrap times. Use 0 to skip bootstrap which might take more time. Default is 0.                 
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. All estimators
#' are presented along with their confidence interval.
#' @export
Incidence <- function(data, allpts = FALSE, size = NULL, knots = 20, nboots = 0,Sdecompose=T){
  
  data_extr<-Change4_c2tlc1t(data)
  extr<-data_extr$extr
  data<-data_extr$data
  
  t<-Create_t(data, allpts , size , knots )
  if(is.null(colnames(data))){colnames(data) <- paste0("site",1:ncol(data))}
  nT <- data[1,]
  T1<-as.numeric(nT[1])
  T2<-as.numeric(nT[2])
  
  
  result_inci<-NULL
  result_inci_CI<-NULL
  result_inci<-Inci(data=data,knots=knots,size=size,Sdecompose=Sdecompose)
  esti<-result_inci
  if(nboots>1){
    boot_sample<-Inci_CreatBootstrapSample(data=data,nboots =nboots)
    sepe <- nrow(boot_sample)
    
    out_boot <- sapply(1:nboots, function(j){
      
     
      d1 <- c(T1,boot_sample[1:(sepe/2),j])
      d2 <- c(T2,boot_sample[(sepe/2+1):sepe,j])
      data_b <- matrix( c(d1,d2), ncol = 2 ,dimnames = list(NULL, names(datap) ))
      data_b <- data_b[rowSums(data_b)>0,] 
      data_b<-as.data.frame(data_b)
      
      out = Inci(data=data_b, allpts, size, knots )
      
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
    esti
   
    
  }
  return(esti)
  
}

#' Mixture diversity of 2 communities
#' @param data a  Sx2 dataframe, the intact assemblage (main) assemblage should be the first column.
#' @param t a matrix specifying the values of t1,t2 where the mixture diversity will be computed.
#' @param nT a matrix specifying the values of T1,T2
#' @return a table: \code{$Mix2} a table of mixture diversity for q=0,1,2.
Inci_Mix <- function(data, t,nT){
  
  #compute the mixture diversity
  #make sure if needed extrapolation?
  tid_inext<- apply(t, 1, function(z) sum(z<=nT))
  datap <- NULL
  tid_in<-NULL
  tid_ext<-NULL
  
  x1<-data[,1]
  x2<-data[,2]
  
  data1<-data[-1,]
  y1<-data1[,1]
  y2<-data1[,2]
 
  if( sum( tid_inext==2) < nrow(t)){
   # datap <- Boot_p(x = data,Bootype = "JADE",datatype = "incidence_freq",zero=T)
    p1_est = Boot_p(x = x1,Bootype = "JADE",datatype = "incidence_freq",zero=T)
    p2_est = Boot_p(x = x2,Bootype = "JADE",datatype = "incidence_freq",zero=T)
    datap = data1
    datap[,1]<-p1_est[1:nrow(datap)]
    datap[,2]<-p2_est[1:nrow(datap)]
  
    tid_in<-which(tid_inext %in% 2)    
    tid_ext<-which(tid_inext %in% 1) 
    t_in<-t[tid_in,]
    t_ext<-t[tid_ext,]
    
    #q=0,1 in
    q01_in<-Dq01_in_inc(x1,x2,t_in)
    q0_in<-q01_in[,1]
    q1_in<-q01_in[,2]
    
    #q=0,1 ext
    q01_ext<-Dq01_ext_inc(x1,x2,t_ext,datap)
    q0_ext<-q01_ext[,1]
    q1_ext<-q01_ext[,2]
    
    q0_mix<-c(q0_in,q0_ext)
    q1_mix<-c(q1_in,q1_ext)
    
    
  }else{
    datap <- NULL
    q01_mix<-Dq01_in_inc(x1,x2,t)
    q0_mix<-q01_mix[,1]
    q1_mix<-q01_mix[,2]
  }
  
  #q=2 in and ext are the same formula
  q2_mix<-Dq2_inci(x1,x2,t)
  
  
  t1<-t[,1]
  t2<-t[,2]
  Mix2<-tibble(t1=rep(t1,3),t2=rep(t2,3),q=rep(c(0,1,2),each=nrow(t)),Diversity =c(q0_mix,q1_mix,q2_mix))
  return(Mix2)
  
}


#' Inci(data, allpts = FALSE, size = NULL, knots = 20 ) for incidence data, comupute composite diversity of any sample and species composition (shared and unique species)
#' @param data a Sx2 dataframe, the intact assemblage (main) assemblage should be the first column.
#' @param allpts specifying whether to compute all combinations of sampling units of two assemblages. Default is FALSE.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param knots the number of points that the mixture diveristy will be computed. Default is 20.
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. 
Inci <- function(data, allpts = FALSE, size = NULL, knots = 20,Sdecompose=T ){
  
  
  data_extr<-Change4_c2tlc1t(data)
  extr<-data_extr$extr
  data<-data_extr$data
  
  t<-Create_t(data, allpts , size , knots )
  
  if(is.null(colnames(data))){colnames(data) <- paste0("site",1:ncol(data))}
  nT <- data[1,]
  T1<-as.numeric(nT[1])
  T2<-as.numeric(nT[2])

  
  #compute diveristy of each single community
  Each_inc<- lapply(1:ncol(data),function(i){
    #mtmp <-  seq(0,max(n[i],n[1]),length.out = knots) %>% round(0)
    ttmp <-  t[,i]
    cname<-names(nT)[i]
    output <- iNEXT(x = data[,i],q = c(0,1,2),datatype = "incidence_freq", size = ttmp,se = F)$iNextEst
    output$qD[is.na(output$qD)] <- 0
    output %>% as_tibble() %>% select(t, q = order, !!cname := qD) 
  }) 
  
 
  #compute the mixture diversity
  Mix2<-Inci_Mix(data, t,nT)
  
  
  
  #merge each and mixture into one dataframe
  each1<-Each_inc[[1]] %>% mutate(t1=t)
  each2<-Each_inc[[2]] %>% mutate(t2=t)
  Result_q012<-Mix2 %>% left_join(each1,by=c("t1","q")) %>% select(-t)
  Result_q012<-Result_q012 %>% left_join(each2,by=c("t2","q")) %>% select(-t)
  Result_q012<- Result_q012 %>% select(t1,t2,q,5:6,Mixture=Diversity) 
  
  Result_q012.1<-Result_q012
  tmp<-Result_q012.1 %>% filter(t2==T2)
  Result_q012.1<-rbind(Result_q012.1,tmp)
  Result_q012.1<-Result_q012.1 %>% arrange(q,-t1,t2)
  
  output0<-Result_q012.1 %>% filter(q==0) %>% select(-q) %>% data.matrix()
  output1<-Result_q012.1 %>% filter(q==1) %>% select(-q) %>% data.matrix()
  output2<-Result_q012.1 %>% filter(q==2) %>% select(-q) %>% data.matrix()
  
  #Species decompose
  if(Sdecompose==T){
    q0_ana<-SDecompose_Inc(data,allpts,size,knots,q=0)
   
  }
  
  colnames(output0)[c(1,2)] = c("m1","m2")
  
  colnames(output1)[c(1,2)] = c("m1","m2")
  colnames(output2)[c(1,2)] = c("m1","m2")
  colnames(q0_ana)[c(1,2)] = c("m1","m2")
 
  return(list(q0 = output0, q1 = output1, q2 = output2, q0_ana = q0_ana,q012 =Result_q012))
}