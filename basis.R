scaleData<-function(x){
  u<-max(x)
  l<-min(x)
  (x-l)/(u-l)
}
BsplineNE<-function(t,degree,knots){
  if(t>1) t<-1
  if(t<0) t<-0
  m<-degree+1
  x <- c(rep(0,m),knots,rep(1,m))
  temp <- sapply(1 : (length(x)-1 ),function(i){as.numeric(t< x[i+1]&& t >= x[i])})
  for (k in 2:m){
    for(i in 1 :(length(x) - k)){
      d = 0;e = 0
      if (temp[i]!= 0 && x[i + k - 1]!= x[i])  d = ((t-x[i])*temp[i])/(x[i + k - 1]-x[i])
      if (temp[i + 1] != 0 && x[i + k] != x[i + 1])  e = ((x[i + k]-t)*temp[i + 1])/(x[i + k]-x[i + 1])
      temp[i] = d + e
    }
  }
  if (t == 1)  temp[length(x) - m]= 1
  basis<-temp[1:(length(x) - m)]
  basis
}
trunc<-function(t,knots){
  c(t,t^2,
    t^3,
    sapply(knots,function(k){if(t<k) 0 else{(t-k)^3}}))
}
fourier<-function(t,freq){
  c(sapply(freq,function(x){sin(pi*x*t)}),
    sapply(freq,function(x){cos(pi*x*t)}))
}
#take 30 equall  y distributed knots, degree = 1  
tp<-function(t,num){
  c(trunc(t,1:(num-3)/(num-2)),fourier(t,1:num),BsplineNE(t,1,1:(num-2)/(num-1)))
}
#4*num=q
base<-function(w,x,i,j,num){tp(x[j],num)[i]*w}