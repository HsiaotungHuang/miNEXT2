#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/

double Hypergeometric_inc(int K, int k, int N,  int n) {
  return exp(Rf_lchoose(K,k)+Rf_lchoose(N-K,n-k)-Rf_lchoose(N,n));
  //return Rf_choose(K,k)*Rf_choose(N-K,n-k)/Rf_choose(N,n);
}

// [[Rcpp::export]]
double fk_inc(int k1, int k2,int T1, int T2, int t1, int t2, NumericVector x1, NumericVector y1){
  
  NumericVector x = x1[(x1>=k1) & (y1>=k2)];
  NumericVector y = y1[(x1>=k1) & (y1>=k2)];
  double output = 0;
  if(sum((x1>=k1) & (y1>=k2))==0){
    output = output;
  }else{
    for(int i = 0; i<x.size(); i++){
      output = output + Hypergeometric_inc(x[i], k1, T1, t1)*Hypergeometric_inc(y[i], k2, T2, t2);
    }
  }
  return output;
}


//q=0,1 in a time
// [[Rcpp::export]]
NumericVector D0_rare(NumericVector xi,NumericVector yi,double t1, double t2){
  int T1 = xi(0);
  int T2 = yi(0);
  xi.erase(0);
  yi.erase(0);
  int u1 = sum(xi);
  int u2 = sum(yi);
  double u1_hat = t1*u1/T1;
  double u2_hat = t2*u2/T2;
  
  
  
  double fk1k2_1 = 0;
  double output1 = 0;
  double fk1k2_0 = 0;
  double output0 = 0;
  
  
  for(int k1 = 0; k1<(t1+1) ; k1++){
    for(int k2 = 0; k2<(t2+1) ; k2++){
      if((k1 == 0) & (k2 == 0)){
        fk1k2_0 = fk1k2_0; 
        fk1k2_1 = fk1k2_1; 
      }else{
        double sub1 = fk_inc(k1, k2,T1, T2, t1, t2, xi, yi);
        fk1k2_0 = fk1k2_0 + sub1;
        fk1k2_1 = fk1k2_1 -((k1+k2)/(u1_hat+u2_hat))*log((k1+k2)/(u1_hat+u2_hat))*sub1;
      }
    }
  }
  
  output0 = fk1k2_0;
  output1 = exp(fk1k2_1);
  
  // if((t1<=T1)&(t2<=T2)){
  //   
  // }else{
  //   if((t1<=T1)&(t2>T2)){
  //     
  //   }else{
  //     
  //   }
  // }
  if((t1==0) & (t2==0)){
    output1 = 0; 
  }
  
  
  NumericVector output = NumericVector::create(output0, output1);
  return output;
}



double h0_inci(double pi1, double pi2, int t1, int t2, int T2 ){
  double output = 0;
  output = pow((1-pi1),t1)*pow((1-pi2),T2)*(1-pow((1-pi2),(t2-T2)));
  return output;
}


// [[Rcpp::export]]
double h0_hat_inci(NumericVector pi1, NumericVector pi2, int t1, int t2, int T1, int T2){
  double output_all= 0; 
  //  double output_sh = 0;
  if(t1>0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h0_inci(pi1_tmp[i],pi2_tmp[i],t1,t2,T2)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h0_inci(0,pi2_tmp[i],t1,t2,T2)/(1-pow(1-pi2_tmp[i], T2));
    }
    output_all = sumsh+sumx0;
    //    output_sh = sumsh;
  }
  else if(t1 == 0){
    NumericVector pi2_tmp = pi2[(pi2>0)];
    double sum2 = 0;
    for(int i=0; i < pi2_tmp.size(); i++){
      sum2 = sum2 + h0_inci(0,pi2_tmp[i],0,t2,T2)/(1-pow(1-pi2_tmp[i], T2));
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}


double Qkt(double pi1, double pi2,int t1,int t2,int k1,int k2){
  double result;
  if((pi1 == 0) & (k1==0)){
    result = Rf_dbinom(k2, t2 , pi2, 0);
  }
  else if ((pi2 == 0) & (k2==0)){
    result = Rf_dbinom(k1, t1 , pi1, 0);
  }
  else{
    result = Rf_dbinom(k1, t1 , pi1, 0)*Rf_dbinom(k2, t2 , pi2, 0);
  }
  if((k1 == 0) & (k2==0)){
    result = 0;
  }
  return result;
}

double h1_inci(double pi1, double pi2, double t1, double t2, double UT1, double UT2, int T1, int T2, int s){
  double U1 = t1*UT1/T1;
  double U2 = t2*UT2/T2;
  double tmp1 = 0;
  double tmp2 = 0;
  for(int k2=s; k2 <= t2; k2++){
    for(int k1=s; k1 <= t1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp1 = 0; }
      else{ tmp1 = tmp1 + (k1+k2)/(U1+U2)*log((k1+k2)/(U1+U2))*Qkt(pi1,pi2,t1,t2,k1,k2); }
    }
  }
  tmp1 = -tmp1;
  U2 = UT2;
  for(int k2=s; k2 <= T2; k2++){
    for(int k1=s; k1 <= t1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp2 = 0; }
      else{ tmp2 = tmp2 + (k1+k2)/(U1+U2)*log((k1+k2)/(U1+U2))*Qkt(pi1,pi2,t1,T2,k1,k2); }
    }
  }
  tmp2 = -tmp2;
  double result = tmp1-tmp2;
  return result;
}


double h1_assem2_inci( double pi2,double t2, double UT2, int T2){
  
  double U2 = t2*UT2/T2;
  double tmp1 = 0;
  double tmp2 = 0;
  for(int k2=1; k2 <= t2; k2++){
    tmp1 = tmp1 + (k2)/(U2)*log((k2)/(U2))*Qkt(0,pi2,0,t2,0,k2);
  }
  tmp1 = -tmp1;
  U2 = UT2;
  for(int k2=1; k2 <= T2; k2++){
    tmp2 = tmp2 + (k2)/(U2)*log((k2)/(U2))*Qkt(0,pi2,0,T2,0,k2); 
  }
  tmp2 = -tmp2;
  double result = tmp1-tmp2;
  return result;
}


// [[Rcpp::export]]
double h1_hat_inci(NumericVector pi1, NumericVector pi2, NumericVector yi1, NumericVector yi2, int t1, int t2, int T1, int T2){
  double output_all = 0;
  //  double output_sh = 0;
  if(t1>0){
    NumericVector pi1_tmp = pi1[(yi1>0) & (yi2>0)];
    NumericVector pi2_tmp = pi2[(yi1>0) & (yi2>0)];
    double UT1 = sum(yi1);
    double UT2 = sum(yi2);
    double sumsh = 0;
    //    double sumsh_p = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1_inci(pi1_tmp[i],pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,0)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
      //      sumsh_p = sumsh_p +  h1(pi1_tmp[i],pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,1)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(yi1==0) & (yi2>0)];
    pi2_tmp = pi2[(yi1==0) & (yi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h1_inci(0,pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,0)/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(yi1>0) & (yi2==0)];
    pi2_tmp = pi2[(yi1>0) & (yi2==0)];
    double sumy0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumy0 = sumy0 +  h1_inci(pi1_tmp[i],0,t1,t2,UT1,UT2,T1,T2,0)/(1-pow(1-pi1_tmp[i], T1));
    }
    output_all = sumsh+sumx0+sumy0;
    //    output_sh = sumsh_p;
    
  }
  else if(t1 == 0){
    NumericVector pi2_tmp = pi2[(yi2>0)];
    double UT2 = sum(yi2);
    double sum2 = 0;
    for(int i=0; i < pi2_tmp.size(); i++){
      sum2 = sum2 +  h1_assem2_inci(pi2_tmp[i],t2,UT2,T2)/(1-pow(1-pi2_tmp[i], T2));
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}


// [[Rcpp::export]]
NumericVector un_inci(NumericVector yi, int T, int t){
  int s = yi.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = 1-exp(Rf_lchoose((T-yi[i]),t)-Rf_lchoose(T,t));
  }
  return(out);
}
// [[Rcpp::export]]
NumericVector sh_inci(NumericVector yi1, NumericVector yi2, int T1, int t1, int T2, int t2){
  int s = yi1.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = (1-exp(Rf_lchoose((T1-yi1[i]),t1)-Rf_lchoose(T1,t1)) * exp(Rf_lchoose((T2-yi2[i]),t2)-Rf_lchoose(T2,t2)) );
  }
  return(out);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
