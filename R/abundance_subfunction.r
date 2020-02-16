#######################################################################################################################
#' R code for q=0,1 2-habitat mixture of rarefactioncurves
#' @param x1 a vector/matrix/list of species abundance frequency
#' @param mm an int matrix for number of individuals for x1,x2,x3 
#' @param q  a numerical value for Dq
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution
#' @examples
#' data(spider)
#' data1<-spider

Dq01_in<-function(x1,x2,mm,q){
  out = apply(mm,1,function(mm) {
    mm1 = mm[1]
    mm2 = mm[2]
    D_share(x1, x2, mm1, mm2, q)
  })
  out
}

#######################################################################################################################
#' R code for abundance data q=2 2-habitat mixture of rarefaction and extrapolation curves
#' @param x1 @param x2 a vector/matrix/list of species abundance frequency
#' @param mm an int matrix for number of individuals for x1,x2 
#' @param n1 @param n2 observations for x1,x2
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution


Dq2<-function(x1,x2,mm,n1,n2){
  out = apply(mm,1,function(mm) {
    m1 = mm[1]
    m2 = mm[2]
    
    if(m1 !=0 || m2!=0){
      
      D2 = 1/(m1+m2)+sum(x1*(x1-1)/(n1*(n1-1)))*m1*(m1-1)/(m1+m2)^2+sum(x2*(x2-1)/(n2*(n2-1)))*m2*(m2-1)/(m1+m2)^2+2*sum(x1*x2/(n1*n2))*m1*m2/(m1+m2)^2
      D2<-1/D2
      
    }
    else if(m1==0 && m2==0)
      D2<-0
  })
  out
}



#######################################################################################################################
#' R code for q=0 2-habitat mixture of  extrapolation curves
#' @param x1 @param x2 a vector/matrix/list of species abundance frequency
#' @param m_ext an int matrix for number of individuals for x1,x2(extra parts)
#' @param datap probability for x1,x2
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution

Dq0_ext<-function(x1,x2,m_ext,n,datap){
  p1_hat = datap[, 1]
  p2_hat = datap[, 2]
  n1<-n[1]
  n2<-n[2]
  mm1<-m_ext[,1]
  mm2s<-m_ext[,2]-n2
  mmext<-cbind(mm1,mm2s)
  h0hat<-apply(mmext,1,function(z) {
      mm1 = z[1]
      mm2s = z[2]
      h0_hat_cpp(p1_hat, p2_hat, mm1, mm2s,n1,n2)
    })
  m1n2<-cbind(mm1,n2)
  Dq0_m1n2<-Dq01_in(x1,x2,m1n2,q=0)
  
  Dq0_ext<-Dq0_m1n2+h0hat
  
  return(Dq0_ext) 
}


#######################################################################################################################
#' R code for q=1 2-habitat mixture of  extrapolation curves
#' @param x1 @param x2 a vector/matrix/list of species abundance frequency
#' @param m_ext an int matrix for number of individuals for x1,x2(extra parts)
#' @param datap probability for x1,x2
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution

Dq1_ext<-function(x1,x2,m_ext,n,datap){
  p1_hat = datap[, 1]
  p2_hat = datap[, 2]
  n1<-n[1]
  n2<-n[2]
  mm1<-m_ext[,1]
  mm2<-m_ext[,2]
  mmext<-cbind(mm1,mm2)
  h1hat<-apply(mmext,1,function(z) {
    print(z)
    mm1 = z[1]
    mm2 = z[2]
    h1_hat_cpp(p1_hat, p2_hat, x1,x2,mm1, mm2,n1,n2)
  })
  m1n2<-cbind(mm1,n2)
  Dq1_m1n2<-Dq01_in(x1,x2,m1n2,q=1)
  zDq1_m1n2<<-Dq1_m1n2
  zh1hat<<-h1hat
  
  Dq1_ext<-Dq1_m1n2*exp(h1hat)
  
  return(Dq1_ext) 
}

#######################################################################################################################
#' R code for q=0 Species composition information in a mixed habitat
#' we can decompose species richness Dq0   of any mixed sample into the sum of three components:
#' the numbers of shared species, species unique to the intact habitat (Assemblage I), and species unique to the modified habitat (Assemblage II). 
#' @param x1 @param x2 a vector/matrix/list of species abundance frequency
#' @param m_ext an int matrix for number of individuals for x1,x2(extra parts)
#' @param datap probability for x1,x2
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution

Dq1_ext<-function(x1,x2,m_ext,n,datap){
  p1_hat = datap[, 1]
  p2_hat = datap[, 2]
  n1<-n[1]
  n2<-n[2]
  mm1<-m_ext[,1]
  mm2<-m_ext[,2]
  mmext<-cbind(mm1,mm2)
  h1hat<-apply(mmext,1,function(z) {
    mm1 = z[1]
    mm2 = z[2]
    h1_hat_cpp(p1_hat, p2_hat, x1,x2,mm1, mm2,n1,n2)
  })
  m1n2<-cbind(mm1,n2)
  Dq1_m1n2<-Dq01_in(x1,x2,m1n2,q=1)
  zDq1_m1n2<<-Dq1_m1n2
  zh1hat<<-h1hat
  
  Dq1_ext<-Dq1_m1n2*exp(h1hat)
  
  return(Dq1_ext) 
}






Chat1_f0Fun <-function(f1, f2, n) {
  if (f2 > 0) {
    f0 <- (n - 1) / n * f1^2 / (2 * f2)
    #C<-1-f0*f1/(n*f0+f1)
    #C <- 1 - f1 / n * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))
    C <- 1 - f1 / n * (n-1)*f1/((n-1)*f1+2*f2)
    
  } else if (f2 == 0 & f1 != 0) {
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
    #C<-1-f0*f1/(n*f0+f1)
    C <- 1 - f1 / n * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
    
  } else {
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
    #f0 <- 0
    C <- 1
  }
  f0 <- ceiling(f0)
  return(c(C, f0))
}


#######################################################################################################################
#' R code for create bootstrap sample: only keep not all 0 species and switch column for n1>n2
#' @importFrom chaoUtility Boot_p
#' @param data a matrix/data.frame for species abundance frequency of 2 sites
#' @param nboots an integer specifying the number of points of m1

#
Abun_CreatBootstrapSample <- function(data, nboots = 0){
  data <- data[rowSums(data)>0,]
  if(sum(data[,1]>0) < sum(data[,2]>0)){
    data = data[,c(2,1)]
  }
  
  if(nboots>1){
    x1 = data[, 1]
    x2 = data[, 2]
    n1 = sum(x1)
    n2 = sum(x2)
    
    
    ##bootstrap p from JADE: two parameters
    # p1 = DetAbu(x1, zero = T)
    # p2 = DetAbu(x2, zero = T)
    # data_p = data
    # data_p[,1]<-p1
    # data_p[,2]<-p2
    # undetec1 <- UndAbu(x1)
    # undetec2 <- UndAbu(x2)
    
    ## bootstrap p from 2013: one parameter
    #p1_est = boot_p_abu(x1)
    #p2_est = boot_p_abu(x2)
    p1_est = Boot_p(x = x1,datatype = "abundance")
    p2_est = Boot_p(x = x2,datatype = "abundance")
    
    data_p = data
    data_p[,1]<-p1_est[1:nrow(data_p)]
    data_p[,2]<-p2_est[1:nrow(data_p)]
    undetec1 <- p1_est[-c(1:nrow(data_p))]
    undetec2 <- p2_est[-c(1:nrow(data_p))]
    
    
    
    f1 = sum(rowSums(data)==1)
    f2 = sum(rowSums(data)==2)
    n=n1+n2
    f0hat<-Chat1_f0Fun(f1,f2,n)[2]
    
    data_p = matrix(data = 0,nrow = f0hat, ncol = 2,
                    dimnames = list(NULL, names(data_p))) %>% rbind(data_p,.)
    zero1 <-  which(data_p[,1]==0)
    zero2 <-  which(data_p[,2]==0)
    
    data_boot <- list() 
    for(k in 1:nboots){
      fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
      fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
      data_pboot <- data_p 
      data_pboot[fill1,1] <- undetec1
      data_pboot[fill2,2] <- undetec2
      nrow(data_pboot)
      # bootx1 <- sapply(data_pboot[,1], function(i) rbinom(n = 1, size = n1, prob = i)) 
      # bootx2 <- sapply(data_pboot[,2], function(i) rbinom(n = 1, size = n2, prob = i)) 
      bootx1<- rmultinom (n = 1, size = n1, prob = data_pboot[,1])
      bootx2<- rmultinom (n = 1, size = n2, prob = data_pboot[,2])
      
      tmp <- sum(bootx1>0) > sum(bootx2>0)
      tmp2 <- ((bootx1==0) & (bootx1>0)) %>% sum 
      while( (tmp==F) | tmp2==0 ){
        fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
        fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
        data_pboot <- data_p 
        data_pboot[fill1,1] <- undetec1
        data_pboot[fill2,2] <- undetec2
        bootx1<- rmultinom (n = 1, size = n1, prob = data_pboot[,1])
        bootx2<- rmultinom (n = 1, size = n2, prob = data_pboot[,2]) 
        tmp <- sum(bootx1>0) > sum(bootx2>0)
        tmp2 <- ((bootx1==0) & (bootx2>0)) %>% sum 
      }
      data_b<-cbind(bootx1,bootx2)
      data_b <- data_b[rowSums(data_b)>0,] 
      colnames(data_b)<-c(colnames(data))
      data_boot[[k]]<-data_b
      
    }
    
  }
  return(data_boot)
}



###for calculate abundance bootstrap confidence interval####

cal_estboot_CI<-function(estboot=NULL,est=NULL){
  out_boot012 <- sapply(1:length(estboot), function(j){
    tmp<-estboot[[j]]
    out012tmp<-matrix(c(tmp$q0[,5], tmp$q1[,5], tmp$q2[,5]), ncol = 3,dimnames = list(NULL,c("q0","q1","q2")))
    
  }, simplify = "array")
  
  
  out_bootana <- sapply(1:length(estboot), function(j){
    tmp<-estboot[[j]]
    outanatmp<-matrix(c(tmp$q0_ana[,3], tmp$q0_ana[,4], tmp$q0_ana[,5]),ncol = 3,dimnames=list(NULL,colnames(tmp$q0_ana)[3:5]))
    
  }, simplify = "array")
  sd_boot_012 <- apply(out_boot012,MARGIN = c(1,2), sd) 
  sd_boot_ana <- apply(out_bootana,MARGIN = c(1,2), sd) 
  
  # mean_boot_012 <- apply(out_boot012,MARGIN = c(1,2), mean) 
  #  mean_boot_ana <- apply(out_bootana,MARGIN = c(1,2), mean) 
  estfinal<-est
  estfinal$q0 <- estfinal$q0 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot_012[,1]) , UCL = (.[,5] + 1.96*sd_boot_012[,1]),s.e.=sd_boot_012[,1])
  estfinal$q1 <- estfinal$q1 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot_012[,2]) , UCL = (.[,5] + 1.96*sd_boot_012[,2]),s.e.=sd_boot_012[,2])
  estfinal$q2 <- estfinal$q2 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot_012[,3]) , UCL = (.[,5] + 1.96*sd_boot_012[,3]),s.e.=sd_boot_012[,3])
  estfinal$q0_ana <- estfinal$q0_ana %>% cbind(., LCL = (.[,3:5] - 1.96*sd_boot_ana[,1:3]) , UCL = (.[,3:5] + 1.96*sd_boot_ana[,1:3]),s.e=sd_boot_ana[,1:3])
  estfinal$q0_ana[ estfinal$q0_ana < 0 ] = 0
  
  estfinal$q012 <- estfinal$q012 %>% mutate(LCL = Mixture-qnorm(0.975)*boot_se,UCL = Mixture+qnorm(0.975)*boot_se,
                            se=boot_se,method = ifelse((m1<=n[1] & m2<=n[2]),"rarefaction","extrapolation"))
  
  return(estfinal)
  
}

Change4_c2slc1s<-function(data){
  data1<-data
  D1 = sum(data1[, 1]>0)
  D2 = sum(data1[, 2]>0)
  if(D1<D2){
    data1 <- data1[,c(2,1)]
  }else{
    data1 <- data1
  }
  return(data1)
}

Create_data1mn<-function(data, knots = 10, size = NULL ){
  data1<-Change4_c2slc1s(data)
  x1 = data1[, 1]
  x2 = data1[, 2]
  
  if(is.null(colnames(data1))){colnames(data1) <- paste0("site",1:ncol(data1))}
  n <- colSums(data1)
  D <- apply(data1,2, function(x) sum(x>0))
  n1<-n[1]
  n2<-n[2]
  
  if(!is.null(size)){
    m1 <- size
    m2 <- n1- m1 
    m = cbind(m1, m2)
  }else{
    if(n1<=n2){
      m2 = round(seq(0, (n1), length.out = knots))
      m1 = n1 - m2
      m = cbind(m1, m2)
     # m2 = round(seq(0, (n2), length.out = knots))
    }else{
      m2 = round(c(seq(0, n2, length.out = knots),
                   seq(n2, (n1), length.out = knots)))
      m1 = n1 - m2
      m = cbind(m1, m2)
      m=unique(m)
    }
  }
  return(list(data1 = data1, m = m,n=n))
}

#######################################################################################################################
#' R code for q=0 Species composition information in a mixed habitat
#' @param data a matrix/data.frame for species abundance frequency of 2 sites
#' @param knots an integer specifying the number of points of m1. Default to be 10.
#' @param size an vector
#' @param q an integer equal to 0
#' @export

SDecompose<-function(data, knots, size,q=0){
  #print(paste("q0_ana  ",Sys.time())) 
  data1mn<-Create_data1mn(data, knots, size)
  data1<-data1mn$data1
  m<-data1mn$m
  n<-data1mn$n
  m1<-m[,1];m2<-m[,2]
  n1<-n[1];n2<-n[2]
  q0_ana<-NULL
  if(q==0){
    m1_i = unique(m1[m1 >= n1-n2])
    m2_i = unique(n1 - m1_i)
    datash = data1[(data1[,1]>0 & data1[,2]>0), , drop=F]
    dataun1 = data1[(data1[,1]>0 & data1[,2]==0), , drop=F]
    dataun2 = data1[(data1[,1]==0 & data1[,2]>0), , drop=F]
    un1 = sapply(m1_i,function(i){
      sum(un_abun(xi = dataun1[,1], n = n1, m = i))
    })
    un2 = sapply(m2_i,function(i){
      sum(un_abun(xi = dataun2[,2], n = n2, m = i))
    })
    sh12 = sapply(1:length(m1_i),function(i){
      sum(sh_abun(xi1 = datash[,1],xi2 = datash[,2], n1 = n1,m1 =  m1_i[i],
                  n2 = n2, m2 = m2_i[i]))
    })
    q0_ana = data.frame(m1 = m1_i, m2 = m2_i,
                        q0_un1 = un1, q0_un2 = un2,
                        q0_sh = sh12)
    colnames(q0_ana)[3:5] = c(paste("Unique to",colnames(data1)[1]),
                              paste("Unique to",colnames(data1)[2]),
                              "Share")  
    
    
  }
  return(q0_ana)
 
}


