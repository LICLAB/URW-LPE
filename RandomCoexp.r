
#the randomwalk method
RandomWalk2igraph<-
function(igraphM,VertexWeight,EdgeWeight=TRUE,gamma=0.7){
if(EdgeWeight==TRUE)
{
 adjM<-get.adjacency(igraphM,attr="weight") # convert igraph object to a weight matrix
}
if(EdgeWeight==FALSE)
{
 adjM<-get.adjacency(igraphM) # convert igraph object to a conventional matrix
}

res<-rw(adjM,VertexWeight,gamma)
return(drop(res))
}

rw<-
function(W,p0,gamma) {
  
   p0<-t(p0)
   p0 <- p0/sum(p0)
   PT <- p0
   
   k <- 0
   delta <- 1

  Ng <- dim(W)[2]
  for (i in 1:Ng) {
      sumr<-sum(W[i,])
      if(sumr==0)
      {
      W[i,] <-numeric(length=length(W[i,]))
      }
      if(sumr>0)
      {
      W[i,] <- W[i,]/sum(W[i,])
      }
    }
    
    
   W <- t(W)
   
   while(delta>1e-10) {
      PT1 <- (1-gamma)*W
      PT2 <- PT1 %*% t(PT)
      PT3 <- (gamma*p0)
      PT4 <- t(PT2) + PT3
      delta <- sum(abs(PT4 - PT))
      PT <- PT4
      k <- k + 1
    }
    PT<-t(PT)
    rownames(PT)<-NULL
    return(PT)
}


#Unsupervised Random Walk method
GetESCCscoreCutNet<-
function(Net,FC){

g <- graph.data.frame(Net[,1:2], directed=F, vertices=NULL)


vetextnames<-get.vertex.attribute(g, "name", index=V(g))


logfold<-c()
for(i in(1:length(vetextnames)))
{
  print(i)
  
  TLF<-FC[which(FC[,1]==vetextnames[i]),6]
  
  if(length(TLF)>0)
  {
   
   logfold<-c(logfold,TLF)
  }
if(length(TLF)==0)
  {
   
   logfold<-c(logfold,0)
  }
}
ID<-vetextnames
logfold[(logfold==Inf)]<-max(logfold[!(logfold==Inf)])
logfold[is.na(logfold)]<-0
logfold[(logfold==(-Inf))]<-min(logfold[!(logfold==(-Inf))])


resNodeW<-RandomWalk2igraph(g,log10(1+log10(1+abs(logfold))),EdgeWeight=FALSE)

result<-data.frame(ID=ID,logfold=logfold,ESCCW=resNodeW)
return(result)
}
