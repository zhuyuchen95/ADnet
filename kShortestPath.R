kshortestpath<-function(netcostmatrix,source,destination,k_paths){
  if(source>nrow(netcostmatrix) | destination>nrow(netcostmatrix)){
    warning("The source or destination are more than the size of netcostmatrix")
    shortestpaths<-list()
    totalcosts<-c()
  }else{
    k<-1
    first<-dijkstra(netcostmatrix,source,destination)
    path<-first$shortestpath
    cost<-first$totalcost
    if (is.null(path)){
      shortestpaths<-NULL
      totalcosts<-NULL
    }else{
      path_number<-1
      P<-list()
      P[[path_number]]<-list(path,cost)
      current_P<-path_number
      size_X<-1
      shortestpaths<-list()
      totalcosts<-c()
      X<-list()
      S<-NULL
      X[[size_X]]<-list(path_number,path,cost)
      S[path_number]<-path[1]
      shortestpaths[[k]]<-path
      totalcosts[k]<-cost
      while ( k < k_paths & size_X != 0){
        for (i in 1 : length(X)) {
          if (X[[i]][[1]] == current_P){
            size_X<-size_X-1
            X[[i]]<-NULL
            break
          }
        }
        P_<-P[[current_P]][[1]]
        w<-S[current_P]
        for(i in 1:length(P_)){
          if(w==P_[i])
            w_index_in_path<-i
        }
        for(index_dev_vertex in w_index_in_path:(length(P_)-1)){
          temp_netcostmatrix<-netcostmatrix
          if(index_dev_vertex>1){
            for (i in 1:(index_dev_vertex-1)) {
              v<-P_[i]
              temp_netcostmatrix[v,]<-Inf
              temp_netcostmatrix[,v]<-Inf
            }
          }
          SP_samesubpath<-list()
          index<-1
          SP_samesubpath[[index]]<-P_
          for (i in 1:length(shortestpaths)) {
            if(length(shortestpaths[[i]])>=index_dev_vertex){
              if (all(P_[1:index_dev_vertex]==shortestpaths[[i]][1:index_dev_vertex])){
                index<-index+1
                SP_samesubpath[[index]]<-shortestpaths[[i]]
              }
            }
          }
          v_<-P_[index_dev_vertex]
          for (j in 1:length(SP_samesubpath)){
            next_<-SP_samesubpath[[j]][index_dev_vertex+1]
            temp_netcostmatrix[v_,next_]<-Inf
          }
          sub_P<-P_[1:index_dev_vertex]
          cost_sub_P<-0
          if(length(sub_P)-1>=1){
            for (i in 1:(length(sub_P)-1)) {#i=1
              cost_sub_P<-cost_sub_P+netcostmatrix[sub_P[i],sub_P[i+1]]#应该是一个总权值
              # cat(c(i,cost_sub_P), "\n")
            }
            # !!!invaild(无效的): cost_sub_P<-newcost_sub_P#目前是1-3的权=2，1-3-4的权为4
            
          }
          #已改正
          second<-dijkstra(temp_netcostmatrix,P_[index_dev_vertex],destination)
          dev_p<-second$shortestpath
          c<-second$totalcost
          if(!is.null(dev_p)){
            path_number<-path_number+1
            if(length(sub_P)-1>=1){
              P[[path_number]]<-list(c(sub_P[1:(length(sub_P)-1)],dev_p),cost_sub_P+c)
            }else{
              P[[path_number]]<-list(dev_p,cost_sub_P+c)
            }
            S[path_number]<-P_[index_dev_vertex]
            size_X<-size_X+1
            X[[size_X]]<-list(path_number,P[[path_number]][[1]],P[[path_number]][[2]])
          }
        }
        if(size_X>0){
          shortestxcost<-X[[1]][[3]]
          shortestx<-X[[1]][[1]]
          if(size_X>=2){
            for(i in 2:size_X){
              if (X[[i]][[3]]<shortestxcost){
                shortestx<-X[[i]][[1]]
                shortestxcost<-X[[i]][[3]]
              }
            }
          }
          current_P<-shortestx
          k<-k+1
          shortestpaths[[k]]<-P[[current_P]][[1]]
          totalcosts[k]<-P[[current_P]][[2]]
        }
      }
      final<-list(shortestpaths=shortestpaths,totalcosts=totalcosts)
      return(final)   
    }
  }
}



