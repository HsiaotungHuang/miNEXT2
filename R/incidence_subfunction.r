

#######################################################################################################################
#' R code for q=0,1 2-habitat mixture of rarefactioncurves
#' @param x1,x2 a vector/matrix/list of species incidence frequency and first row is T
#' @param t an int matrix for number of individuals for x1,x2 
#' @return  a vector/matrix/list with species diversity for q=0 and 1
#' @examples
#' data(spider)
#' data1<-spider


Dq01_in_inc<-function(x1,x2,t){
  out = apply(t,1,function(tt) {
    t1 = tt[1]
    t2 = tt[2]
    D0_rare(x1, x2, t1,t2)
  })
  t(out)
}


#######################################################################################################################
#' R code for q=0 2-habitat mixture of  extrapolation curves
#' @param x1 @param x2 a vector/matrix/list of species abundance frequency
#' @param m_ext an int matrix for number of individuals for x1,x2(extra parts)
#' @param datap probability for x1,x2
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution




Dq01_ext_inc<-function(x1,x2,t_ext,datap){
  p1_hat = datap[, 1]
  p2_hat = datap[, 2]
  T1<-x1[1]
  T2<-x2[1]
  t1<-t_ext[,1]
  h0hat<-apply(t_ext,1,function(z) {
    t1_ex = z[1]
    t2_ex = z[2]
    h0_hat_inci(p1_hat, p2_hat, t1_ex, t2_ex,T1,T2)
  })
  
  y1<-x1[-1]
  y2<-x2[-1]
  h1hat<-apply(t_ext,1,function(z) {
    t1_ex = z[1]
    t2_ex = z[2]
    h1_hat_inci(p1_hat, p2_hat,y1,y2, t1_ex, t2_ex,T1,T2)
  })
  
  
  t1T2<-cbind(t1,T2)
  Dq01_t1T2<-Dq01_in_inc(x1,x2,t1T2)
 
  Dq0_ext<-Dq01_t1T2[,1]+h0hat
  Dq1_ext<-Dq01_t1T2[,2]*exp(h1hat)
  
  return(cbind(Dq0_ext,Dq1_ext))
}


#######################################################################################################################
#' R code for incidence data q=2 2-habitat mixture of rarefaction and extrapolation curves
#' @param x1 @param x2 a vector/matrix/list of species incidence
#' @param t1 @param t2 plots for x1,x2
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution


Dq2_inci<-function(x1,x2,t){
  T1<-x1[1]
  T2<-x2[1]
  y1<-x1[-1]
  y2<-x2[-1]
  out = apply(t,1,function(mm) {
    t1 = mm[1]
    t2 = mm[2]
   
    if(t1 !=0 || t2!=0){
      u1 = sum(y1)
      u2 = sum(y2)
      u1_hat = t1*u1/T1;
      u2_hat = t2*u2/T2;
      D2 = 1/(u1_hat+u2_hat)+sum(y1*(y1-1)/(T1*(T1-1)))*t1*(t1-1)/(u1_hat+u2_hat)^2+sum(y2*(y2-1)/(T2*(T2-1)))*t2*(t2-1)/(u1_hat+u2_hat)^2+2*sum(y1*y2)*t1*t2/(T1*T2)/(u1_hat+u2_hat)^2
      D2<-1/D2
      
    }
    else if(t1==0 && t2==0)
      D2 = 0
  })
  out
  
}



Change4_c2tlc1t<-function(data){
  data <- data[rowSums(data)>0,]
  if(sum(data[,1]>0) < sum(data[,2]>0)){
    data = data[,c(2,1)]
  }
  
  ts = as.numeric(data[1,])
  T1 = ts[1]
  T2 = ts[2]
  
  d0ob = c(sum(data[,1]>0),sum(data[,2]>0))
  if(T1 == T2){
    extr = 0
    if(d0ob[1]<d0ob[2]){
      data = data[,c(2,1)]
    }
  }else if(((T1-T2)*(d0ob[1]-d0ob[2])) > 0){
    extr = 1
    if(T1<T2){
      data = data[,c(2,1)]
    }
  }else if (((T1-T2) * (d0ob[1]-d0ob[2])) < 0 ){
    extr = 0
    if(T1>T2){
      data = data[,c(2,1)]
    }
  }
  return(list(extr=extr,data=data))
}


Create_t<-function(data, allpts = FALSE, size = NULL, knots = 20  ){
  ts = as.numeric(data[1,])
  T1 = ts[1]
  T2 = ts[2]
  if(allpts == TRUE){
    t1 = seq(T1,0)
    t2 = T1-t1
  }else if (allpts == FALSE & !is.null(size)){
    if(max(size)>T1)
      stop(paste0("Because the community ",colnames(data)[1]," has higher divrsity, its sampling units should be viewed main."))
    
    t1 = sort(c(size,(T1-T2),T1),decreasing = T)
    t1 = t1[!duplicated(t1[t1>=0])]
    t2 = T1-t1
  }else if (allpts == FALSE & is.null(size)){
    t1 = round(seq(T1,0,length.out = knots),0)
    t1 = sort(c(t1,(T1-T2)),decreasing = T)
    t1 = t1[t1>=0] %>% unique
    t2 = T1-t1
  }
  t = cbind(t1,t2)
  return(t)
}


#######################################################################################################################
#' R code for q=0 Species composition information in a mixed habitat
#' @param data a matrix/data.frame for species abundance frequency of 2 sites
#' @param knots an integer specifying the number of points of m1. Default to be 10.
#' @param size an vector
#' @param q an integer equal to 0
#' @export

SDecompose_Inc<-function(data, allpts = FALSE, size = NULL, knots = 20,q=0){
  data_extr<-Change4_c2tlc1t(data)
  extr<-data_extr$extr
  data<-data_extr$data
  
  t<-Create_t(data, allpts , size , knots )
  if(is.null(colnames(data))){colnames(data) <- paste0("site",1:ncol(data))}
  nT <- data[1,]
  T1<-as.numeric(nT[1])
  T2<-as.numeric(nT[2])
  t1<-t[,1]
  t2<-t[,2]
  
  t1_in = sort(t1[t1>=(T1-T2)],decreasing = T) 
  t2_in = t1_in[1] - t1_in
  
  data1 = data[-1,]
  q0_ana<-NULL
  if(q==0){
    datash = data1[(data1[,1]>0 & data1[,2]>0), , drop = F]
    dataun1 = data1[(data1[,1]>0 & data1[,2]==0), , drop = F]
    dataun2 = data1[(data1[,1]==0 & data1[,2]>0), , drop = F]
    
    un1 = sapply(t1_in,function(i){
      sum(un_inci(yi = dataun1[,1], T1, i))
    })
    un2 = sapply(t2_in,function(i){
      sum(un_inci(yi = dataun2[,2], T = T2, t = i))
    })
    sh12 = sapply(1:length(t1_in),function(i){
      sum(sh_inci(yi1 = datash[,1],yi2 = datash[,2], T1 = T1, t1 =  t1_in[i],
                  T2 = T2, t2 = t2_in[i]))
    })
    q0_ana = data.frame(m1 = t1_in, m2 = t2_in,
                        q0_un1 = un1, q0_un2 = un2,
                        q0_sh = sh12)
    
    
    colnames(q0_ana)[3:5] = c(paste0("Uni.",colnames(data1)[1]),
                              paste0("Uni.",colnames(data1)[2]),
                              "Share")
 
  }
  
  return(q0_ana)
  
}


# Functions for bootstrap
#' @importFrom chaoUtility Boot_p
#' @importFrom stats rmultinom
Inci_CreatBootstrapSample <- function(data, nboots = 0){
  data <- data[rowSums(data)>0,]
  data_boot<-NULL
  if(nboots>1){
    
    nT1 <- data[1,1]
    nT2 <- data[1,2]
    datap = data[-1,]
    
    x1 = data[, 1]
    x2 = data[, 2]
    y1 <- x1[-1]
    y2 <- x2[-2]
    p1_est = Boot_p(x = x1,datatype = "incidence")
    p2_est = Boot_p(x = x2,datatype = "incidence")
    datap[,1]<-p1_est[1:nrow(datap)]
    datap[,2]<-p2_est[1:nrow(datap)]
    undetec1 <- p1_est[-c(1:nrow(datap))]
    undetec2 <- p2_est[-c(1:nrow(datap))]
    
    Q1 = sum(rowSums(data[-1,])==1)
    Q2 = sum(rowSums(data[-1,])==2)
    nT = nT1+nT2
    Q0.hat<-Chat1_f0Fun(Q1,Q2,nT)[2]
    datap = matrix(data = 0,nrow = Q0.hat, ncol = 2,
                   dimnames = list(NULL, names(datap))) %>% rbind(datap,.)
    zero1 <-  which(datap[,1]==0)
    zero2 <-  which(datap[,2]==0)
    data_boot <- sapply(1:nboots,function(k){
      fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
      fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
      datapp <- datap 
      datapp[fill1,1] <- undetec1
      datapp[fill2,2] <- undetec2
      data1 <- sapply(datapp[,1], function(i) rbinom(n = 1, size = nT1, prob = i)) 
      data2 <- sapply(datapp[,2], function(i) rbinom(n = 1, size = nT2, prob = i)) 
      tmp <- sum(data1>0) > sum(data2>0)
      tmp2 <- ((data1==0) & (data2>0)) %>% sum 
      while( (tmp==F) | tmp2==0 ){
        fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
        fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
        datapp <- datap 
        datapp[fill1,1] <- undetec1
        datapp[fill2,2] <- undetec2
        data1 <- sapply(datapp[,1], function(i) rbinom(n = 1, size = nT1, prob = i)) 
        data2 <- sapply(datapp[,2], function(i) rbinom(n = 1, size = nT2, prob = i)) 
        tmp <- sum(data1>0) > sum(data2>0)
        tmp2 <- ((data1==0) & (data2>0)) %>% sum 
      }
      c(data1,data2)
      
    })
    
  }
  return(data_boot)
}


