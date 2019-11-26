dijkstra<-function(netcostmatrix,s,d){
  
  n<-nrow(netcostmatrix)
  farthestprevhop<-c()
  farthestnexthop<-c()
  parent<-c()
  distance<-c()
  visited<-c()
  for (i in 1:n) {
    farthestprevhop[i]<-i
    farthestnexthop[i]<-i
  }
  
  visited[1:n]<-FALSE
  distance[1:n]<-Inf
  parent[1:n]<-0
  distance[s]<-0
  
  for (i in 1:(n-1)) {
    temp<-c()
    for (h in 1:n) {
      if( !visited[h]) temp<-c(temp,distance[h]) else  temp<-c(temp,Inf)
    }
    t<-min(temp)
    u<-which(temp== min(temp), arr.ind = TRUE)
    if(length(u)>1){
      u<-u[1]
    }
    visited[u]<-TRUE
    for (v in 1:n) {
      if( distance[u]+netcostmatrix[u,v]< distance[v]){
          distance[v]<-(distance[u]+netcostmatrix[u,v])
          parent[v]<-u
        }
    }
  }
  shortestpath<-c()
  if (parent[d] != 0){
    t<-d
    shortestpath<-c(d)
    while(t!=s){
      p<-parent[t]
      shortestpath<-c(p,shortestpath)
      
      if(netcostmatrix[t,farthestprevhop[t]] < netcostmatrix[t,p])
        farthestprevhop[t]<-p
      if (netcostmatrix[p,farthestnexthop[p]] < netcostmatrix[p,t])
        farthestnexthop[p]<-t
      t<-p
    }
  }
  totalcost<-distance[d]
  result<-list(shortestpath=shortestpath,totalcost=totalcost)
  return(result)
}
