source("basis.R")
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
  x,
  num
){
  val<-0
  if(root$totalCount==1){
    val<-base(root$w,x,root$i,root$j,num)
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
  x,    #matrix each row indicates a single x
  num
){
  val<-rep(0,nrow(x))
  if(root$totalCount==1){
    val<-apply(x,1,
               function(xx){
                 base(root$w,xx,root$i,root$j,num)
               })
  }else{
    if(root$op=="+"){
      val<-evaluateMultiple(root$children[[1]],x,num)+
        evaluateMultiple(root$children[[2]],x,num)
    }else if(root$op=="*"){
      val<-evaluateMultiple(root$children[[1]],x,num)*
        evaluateMultiple(root$children[[2]],x,num)
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
  path,
  #path from the root     to the node where to insert new node
  num #4num=q
){
  opp<-function(i){x<-1;if(i==1){x<-2};x}
  pathNode<-path
  root<-tree$root
  n<-dim(x)[1]
  b<-rep(0,n)
  k<-rep(1,n)
  while(length(pathNode)!=0){
    tempt<-evaluateMultiple(root$children[[opp(pathNode[1])]],x,num)
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
# update a single weight as a LS estimation
adjustWeights<-function(tree,#tree structure, 
                                         #Leaves must have i,j,w attributes
                                         #Nodes must have op attributes
                                         x,
                                         #each row is a single x
                                         y,
                                         #a vector, the length must equal to nrow(x)
                                         path,
                                         constraint,
                                         num){
  q<-length(tp(.1,num))
  sumy2<-crossprod(y)
    result<-computeRegression(tree,x,path,num)
    node<-tree$Climb(position=path)
    i<-node$i
    j<-node$j
    intcp<- result$b
    z<-cbind(apply(x,1,function(xx){
          base(1,xx,i,j,num)
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
        #update the tree and return the intercept
        diff<-(w-node$w)^2
        node$w<-w
  list(intercept=tt[2],difference=diff)
}
# Update weights by coordinate descent and return intercept
coordianteDescentWeights<-function(tree,#tree structure, 
                                   #Leaves must have i,j,w attributes
                                   #Nodes must have op attributes
                                   x,
                                   #each row is a single x
                                   y,
                                   #a vector, the length must equal to nrow(x)
                                   constraint,
                                   num){
  treeTemp<-Clone(tree,attributes = TRUE)
  diff<-100
  diffCutoff<-10^-3
  iitersNum<-10
  iiters<-1
  intcpt<-0
  while(diff>diffCutoff & iiters<iitersNum){
    for(node in treeTemp$leaves){
      update<-adjustWeights(treeTemp, x,y,allPath(treeTemp$root,node),1,num)
      diff<-diff+update$difference/length(treeTemp$leaves)
      intcpt<-update$intercept
    }
    iiters<-iiters+1
  }

  list(intercept = intcpt, tree = treeTemp) 
}

tryInsertNode<-function(tree,#tree structure, 
                        #Leaves must have i,j,w attributes
                        #Nodes must have op attributes
                        x,
                        #each row is a single x
                        y,
                        #a vector, the length must equal to nrow(x)
                        path,
                        constraint,
                        num){
  q<-length(tp(.1,num))
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
          base(1,xx,i,j,num)
        }),rep(1,dim(x)[1]))
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
    result<-computeRegression(tree,x,path,num)
    if(length(path)>0){
      node<-tree$Climb(position=path)
    }else{node<-tree$root}
    tempt<-evaluateMultiple(node,x,num)
    intcp<- result$b+tempt*result$k
    #
    for(i in 1:q){
      for(j in 1:dim(x)[2]){
        treeTemp<-Clone(tree,attributes = TRUE)
        if(length(path)>0){
          nodeTemp<-treeTemp$Climb(position=path)
        }else{nodeTemp<-treeTemp$root}
        z<-cbind(apply(x,1,function(xx){
          base(1,xx,i,j,num)
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
          base(1,xx,i,j,num)
        })*u,rep(1,dim(x)[1]))
        sumz2<-crossprod(z)
        sumzyi<-crossprod(z,y-result$b)
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
addLeaf<-function(tree,#tree structure, 
                  #Leaves must have i,j,w attributes
                  #Nodes must have op attributes
                  x,
                  #each row is a single x
                  y,
                  test_x,
                  test_y,
                  #a vector, the length must equal to nrow(x)
                  constraint,
                  num){
  r_best<-10^5
  if(is.null(tree)){
    result<-tryInsertNode(tree,as.matrix(x),y,c(),1,num)
    tree<-Node$new("1")
    tree$w<-result$w
    tree$i<-result$i
    tree$j<-result$j
    intercept<-result$intercept
  }else{
    for(node in tree$leaves){
      result<-tryInsertNode(tree,x,y,allPath(tree$root,node),1,num)
      if(result$RMSE<r_best){
        result_best<-result
        r_best<-result$RMSE
        node_best<-node
      }
    }
    t<-insert(
      node_best,
      result_best$i,    #index of basis function
      result_best$j,    #index of x
      result_best$w,
      result_best$op,
      tree$totalCount #total number of nodes of current tree
    )
    updateWeights<-coordianteDescentWeights(tree,x,y,1,num)
    tree<-updateWeights$tree
    intercept<-updateWeights$intercept
  }
  y_hat<-evaluateMultiple(tree,test_x,num)+intercept
  test<-sqrt(crossprod(y_hat-test_y)/length(test_y))
  y_t_hat<-evaluateMultiple(tree,x,num)+intercept
  training<-sqrt(crossprod(y_t_hat-y)/length(y))
  print(paste("intercept",intercept))
  print("training,test")
  print(c(training,test))
  print(tree,"i","j","w","op")
  list(tree=tree,intercept=intercept,rmse_train=training,rmse_test=test)
}
constructTree<-function(x,y,test_x,test_y,constraint,num,k,run_id){
  tree<-NULL
  itersNum<-k
  rMSE_iters<-rep(0,itersNum)
  RMSET_iters<-rep(0,itersNum)
  for(i in 1:itersNum){
    print(i)
    newLeaf<-addLeaf(tree,x,y,test_x,test_y,constraint,num)
    tree<-newLeaf$tree
    rMSE_iters[i]<-newLeaf$rmse_train
    RMSET_iters[i]<-newLeaf$rmse_test
  }
  print("training RMSE after transformation")
  print(rMSE_iters)
  print("test RMSE after transformation")
  print(RMSET_iters)
  result<-rbind(rMSE_iters,RMSET_iters)
  rownames(result)<-paste(c("training","test"),rep(dim(x)[1],2),rep(dim(x)[2],2),rep(num,2),rep(k,2),rep(run_id,2),sep="-")
}

