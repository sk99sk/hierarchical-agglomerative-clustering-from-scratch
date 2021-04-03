library(ISLR)
nci.labs=NCI60$labs
nci.data=NCI60$data
scaled.data=scale(nci.data)
data.dist = dist(scaled.data)
unscaled.dist = dist(nci.data)


iorder = function(m)
{
  N = nrow(m) + 1
  iorder = rep(0,N)
  iorder[1] = m[N-1,1]
  iorder[2] = m[N-1,2]
  loc = 2
  for(i in seq(N-2,1))
  {
    for(j in seq(1,loc))
    {
      if(iorder[j] == i)
      {
        iorder[j] = m[i,1]
        if(j==loc)
        {
          loc = loc + 1
          iorder[loc] = m[i,2]
        } else
        {
          loc = loc + 1
          for(k in seq(loc, j+2)) iorder[k] = iorder[k-1]
          iorder[j+1] = m[i,2]
        }
      }
    }
  }
  -iorder
}


hc = function(d, method=c("single","complete","average","median"))
{
  if(!is.matrix(d)) d = as.matrix(d)

  method_fn = switch(match.arg(method),
                     single   = min,
                     complete = max,
                     average  = mean,
                     median  = median)
  N = nrow(d)
  diag(d)=Inf
  n = -(1:N)      

  m = matrix(0,nrow=N-1, ncol=2)   
  h = rep(0,N-1)                   
  for(j in seq(1,N-1))
  {

    h[j] = min(d)

    i = which(d - h[j] == 0, arr.ind=TRUE)

    i = i[1,,drop=FALSE]
    p = n[i]

    p = p[order(p)]
    m[j,] = p

    grp = c(i, which(n %in% n[i[1,n[i]>0]]))
    n[grp] = j
    r = apply(d[i,],2,method_fn)

    d[min(i),] = d[,min(i)] = r
    d[min(i),min(i)]        = Inf
    d[max(i),] = d[,max(i)] = Inf
  }

  structure(list(merge = m, height = h, order = iorder(m),
                 labels = nci.labs, method = method, 
                 call = match.call(), dist.method = "euclidean"), 
            class = "hclust")
}




hc_cent = function(d,data)
{
  d = as.matrix(d)
  max_i <- numeric(0)
  N = nrow(d)
  diag(d)=Inf
  n = -(1:N)
  n_new <- t(cbind(n,data))
  
  m = matrix(0,nrow=N-1, ncol=2)   
  h = rep(0,N-1)                   
  for(j in seq(1,N-1))
  {
    
    h[j] = min(d)
    
    i = which(d - h[j] == 0, arr.ind=TRUE)
    i = i[1,,drop=FALSE]
    p = n[i]
    p = p[order(p)]
    m[j,] = p
    
    grp = c(i, which(n %in% n[i[1,n[i]>0]]))
    n[grp] = j
    for(k in grp){
      n_new[1,k] = j
    }
    centroid <- c()
    
    
    for(t in 2:nrow(n_new)){
      sum1 = 0
      c1 = 0
      for (y in 1:ncol(n_new)){
        if(n_new[1,y]==j){
          sum1 = sum1 + n_new[t,y]
          c1 = c1 +1
        }
        
      }
      
      cent = sum1/c1
      centroid <- append(centroid,cent)
      
    }
    
    r <- c()
    
    
    for(z in seq(1:ncol(n_new))){
      if(n_new[1,z]==j){
        for(x in 2:nrow(n_new)){
          n_new[x,z] = centroid[x-1]
        }
      }
    }
    
    euclidean <- function(a, b) sqrt(sum((a - b)^2))
    
    for(l in 1:ncol(n_new)){
      temp = n_new[2:nrow(n_new),l]
      r <- append(r,euclidean(temp,centroid))
      
    }
    
    
    max_i <- c(max_i,max(i))
    d[min(i),] <- d[,min(i)] <- r
    d[min(i),min(i)] = Inf 
    d[max_i,] = d[,max_i] = Inf
    
  }
  
  structure(list(merge = m, height = h, order = iorder(m),
                 labels = nci.labs,method = "centroid", 
                 call = match.call(), dist.method = "euclidean"), 
            class = "hclust")
}



plot(hc_cent(unscaled.dist,nci.data))
#plot(hc(data.dist,method = 'average'))
