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
tp<-function(t){
  c(BsplineNE(t,1,1:9/10),trunc(t,1:9/10),fourier(t,1:10))
}
q<-length(tp(.1))
base<-function(w,x,i,j){tp(x[j])[i]*w}
insert<-function(
  node, #node to insert new operation
  i,    #index of basis function
  j,    #index of x
  w,
  operation="",
  totalNum #total number of nodes of current tree
){
  if(is.null(node)){node<-Node$new(paste(i,",",j,sep=""))}else{
    childi<-node$i
    childj<-node$j
    childw<-node$w
    childL<-node$AddChild(paste(totalNum+1))
    childL$i<-childi
    childL$j<-childj
    childL$w<-childw
    childR<-node$AddChild(paste(totalNum+2))
    childR$i<-i
    childR$j<-j
    childR$w<-w
    node$RemoveAttribute("i")
    node$RemoveAttribute("j")
    node$RemoveAttribute("w")
    node$op<-operation
  }
  node
}
#root: root of tree/subtree
#node
#return path from the root to the node 
#i.e. c(1,2,1,2)
allPath<-function(root,node){
  path<-rep(0,node$level-root$level)
  iters<-length(path)
  current<-node
  while(iters>0){
    path[iters]<-current$position
    current<-current$parent
    iters<-iters-1
  }
  path
}

#Leaves must have i,j,w attributes
#Nodes must have op attributes
#Return a double
evaluateSingle<-function(
  root,#root of the tree/subtree for evaluation 
  x
){
  val<-0
  if(root$totalCount==1){
    val<-base(root$w,x,root$i,root$j)
  }else{
    if(root$op=="+"){
      val<-evaluateSingle(root$children[[1]],x)+evaluateSingle(root$children[[2]],x)
    }else if(root$op=="*"){
      val<-evaluateSingle(root$children[[1]],x)*evaluateSingle(root$children[[2]],x)
    }else{
      print("operation attribute not exist!")
    }
  }
  val
}
#Return a vector of which the length is the same as nrow(x)
evaluateMultiple<-function(
  root,#root of the tree/subtree for evaluation
  x    #matrix each row indicates a single x
){
  val<-rep(0,nrow(x))
  if(root$totalCount==1){
    val<-apply(x,1,
               function(xx){
                 base(root$w,xx,root$i,root$j)
               })
  }else{
    if(root$op=="+"){
      val<-evaluateMultiple(root$children[[1]],x)+
        evaluateMultiple(root$children[[2]],x)
    }else if(root$op=="*"){
      val<-evaluateMultiple(root$children[[1]],x)*
        evaluateMultiple(root$children[[2]],x)
    }else{
      print("operation attribute not exist!")
    }
  }
  val
}

computeRegression<-function(
  tree,#tree structure, 
  #Leaves must have i,j,w attributes
  #Nodes must have op attributes
  x,
  #each row is a single x
  #a vector, the length must equal to nrow(x)
  path
  #path from the root     to the node where to insert new node
){
  opp<-function(i){x<-1;if(i==1){x<-2};x}
  pathNode<-path
  root<-tree$root
  n<-dim(x)[1]
  b<-rep(0,n)
  k<-rep(1,n)
  while(length(pathNode)!=0){
    tempt<-evaluateMultiple(root$children[[opp(pathNode[1])]],x)
    if(root$op=="+"){
      b<-b+tempt*k
    }else if(root$op=="*"){
      k<-k*tempt
    }else{
      print("Operation Attribute Not Exist!")
    }
    root<-root$children[[pathNode[1]]]
    pathNode<-pathNode[-1]
  }
  list(k=k,b=b)  
}
tryInsertNode<-function(tree,#tree structure, 
                        #Leaves must have i,j,w attributes
                        #Nodes must have op attributes
                        x,
                        #each row is a single x
                        y,
                        #a vector, the length must equal to nrow(x)
                        path,
                        constraint){
  RMSE_best<-10^8
  RMSE<-RMSE_best
  i_best<--1
  j_best<--1
  op_best<-""
  w_best<-1000
  intercept_best<-1000 
  sumy2<-crossprod(y)
  if(is.null(tree)){
    for(i in 1:q){
      for(j in 1:dim(x)[2]){
        z<-cbind(apply(x,1,function(xx){
          base(1,xx,i,j)
        }),rep(1,dim(x)[1]))
        # sumz2<-crossprod(z)
        # if(abs(det(sumz2))>1e-10){
        #   sumzy<-crossprod(z,y)
        #   tempt<-solve(sumz2,sumzy)
        #   tempt<-tryCatch(
        #     solve(sumz2,sumzy),
        #     error=function(e) e,
        #     warning=function(w) w)
        #   if(inherits(tempt, "error")){
        #     tempt<-Solve(sumz2,sumzy)
        #   }
        tt<-tryCatch(
          optim(rep(0,2),
                function(c){crossprod(y-z%*%c)}, 
                method = "L-BFGS-B",
                lower=c(-constraint,-Inf),
                upper=c(constraint,Inf))$par,
          error=function(e) e,
          warning=function(w) w)
        if(inherits(tt, "error")){
          tt<-c(0,mean(y))
        }
        w<-tt[1]
        intercept<-tt[2]
        RMSE<-10^10
        try(
          RMSE<-sqrt(crossprod(y-z%*%tt)/length(y)),silent=TRUE
        )
          
          if(RMSE<RMSE_best){
            RMSE_best<-RMSE
            i_best<-i
            j_best<-j
            op_best<-""
            w_best<-w
            intercept_best<-intercept
          }
          }
      }
    }else{
    result<-computeRegression(tree,x,path)
    if(length(path)>0){
      node<-tree$Climb(position=path)
    }else{node<-tree$root}
    tempt<-evaluateMultiple(node,x)
    intcp<- result$b+tempt*result$k
    #
    for(i in 1:q){
      for(j in 1:dim(x)[2]){
        treeTemp<-Clone(tree,attributes = TRUE)
        if(length(path)>0){
          nodeTemp<-treeTemp$Climb(position=path)
        }else{nodeTemp<-treeTemp$root}
        z<-cbind(apply(x,1,function(xx){
          base(1,xx,i,j)
        })*result$k,rep(1,dim(x)[1])) 
        sumz2<-crossprod(z)
        sumzyi<-crossprod(z,y-intcp)
        tt<-tryCatch(
          optim(rep(0,2),
                function(c){crossprod(y-intcp-z%*%c)}, 
                method = "L-BFGS-B",
                lower=c(-constraint,-Inf),
                upper=c(constraint,Inf))$par,
          error=function(e) e,
          warning=function(w) w)
        if(inherits(tt, "error")){
          tt<-c(0,mean(y-intcp))
        }
        w<-tt[1]
        intercept<-tt[2]
        RMSE<-10^10
        try(RMSE<-sqrt(crossprod(y-intcp-z%*%tt)/length(y)),
            silent=TRUE)
        if(RMSE<RMSE_best){  
          RMSE_best<-RMSE
          i_best<-i
          j_best<-j
          op_best<-"+"
          w_best<-w
          intercept_best<-intercept
        }
      }
    }
    #Potential nodes for multiplication
    u<-result$k*tempt
    for(i in 1:q){
      for(j in 1:dim(x)[2]){
        treeTemp<-Clone(tree,attributes = TRUE)
        if(length(path)>0){
          nodeTemp<-treeTemp$Climb(position=path)
        }else{nodeTemp<-treeTemp$root}
        z<-cbind(apply(x,1,function(xx){
          base(1,xx,i,j)
        })*u,rep(1,dim(x)[1]))
        sumz2<-crossprod(z)
        sumzyi<-crossprod(z,y-result$b)
        #if(abs(det(sumz2))>1e-10){
        # tt<-solve(sumz2,sumzyi)
        #}else{
        # tt<-optim(rep(0,2),
        #           function(c){crossprod(y-intcp-z%*%c)}, 
        #           method = "L-BFGS-B",
        #           lower=c(-constraint,-Inf),
        #           upper=c(constraint,Inf))$par
        #}
        tt<-tryCatch(
          optim(rep(0,2),
                function(c){crossprod(y-intcp-z%*%c)}, 
                method = "L-BFGS-B",
                lower=c(-constraint,-Inf),
                upper=c(constraint,Inf))$par,
          error=function(e) e,
          warning=function(w) w)
        if(inherits(tt, "error")){
          tt<-c(0,mean(y-intcp))
        }
        w<-tt[1]
        intercept<-tt[2]
        RMSE<-10^10
        try(
          RMSE<-sqrt(crossprod(y-result$b-z%*%tt)/length(y)),silent=TRUE
        )
        
        if(RMSE<RMSE_best){
          RMSE_best<-RMSE
          i_best<-i
          j_best<-j
          op_best<-"*"
          w_best<-w
          intercept_best<-intercept
        }
      }
    }
  }
  list(i=i_best,j=j_best,RMSE=RMSE_best,op=op_best,w=w_best,intercept=intercept_best)
}

