


#################Data generating function: continuous outcome########################

DATAGEN_cont<-function(n){
  
  A<-rbinom(n,1,0.5)
  
  XDATA1<-list()
  for (k in 1:n){
    XDATA1[[k]]<-rnorm(500,mean=(-1*A[k]),sd=1)}
  
  XDATA2<-list()
  for (k in 1:n){
    XDATA2[[k]]<-rnorm(500,mean=(-1*A[k]),sd=1)}
  
  XDATA2M1<-list()
  for (k in 1:n){
    XDATA2M1[[k]]<-XDATA2[[k]]+rnorm(500,mean=0,sd=0.1)}
  
  YDATA<-list()
  t<-seq(from=0.01,to=5,by=0.01)
  
  for (k in 1:n){
    YDATA[[k]]<-rnorm(500,mean=(5-0.5*t+XDATA1[[k]]+XDATA2[[k]]-0.5*XDATA1[[k]]*XDATA2[[k]]-2*A[k]),sd=0.5)}
  
  VDATA<-list()
  for (k in 1:n){
    VDATA[[k]]<-rbinom(500,1,pmin(exp(-3.05-2*t+0.5*XDATA1[[k]]+0.5*XDATA2[[k]]+A[k]+0.5*XDATA1[[k]]*XDATA2[[k]]+0.3*YDATA[[k]]),1))}
  
  
  DATA<-data.frame(ID=rep(1:n,each=500),t_start=rep(seq(from=0,to=4.99,by=0.01),n),t_stop=rep(seq(from=0.01,to=5,by=0.01),n),
                   status=unlist(VDATA),Y=unlist(YDATA),X=rep(A,each=n),Z1=unlist(XDATA1),Z2=unlist(XDATA2))
  
  DATA$Z1Z2=DATA$Z1*DATA$Z2
  
  
  DATA}




#### function to estimate the balancing weights in the proposed IIWEs of Yiu and Su (2022), incorporating selection function
bal_fit_fun_sa<-function(MAT,cons, offs)
{
  library(nleqslv)
  bal_fun<-function(gamma){
    wei<-exp(colSums(gamma*t(MAT)))*offs
    colSums(wei*MAT)-cons}
  
  jac_fun_bal<-function(gamma){
    wei<-exp(colSums(gamma*t(MAT)))*offs
    t(MAT)%*%(MAT*wei)}
  
  lambda<-nleqslv(rep(0,length(cons)),fn=bal_fun,jac=jac_fun_bal,method="Broyden")$x
  exp(colSums(lambda*t(MAT)))*offs
}




