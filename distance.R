product<-function(tab){
  r=rowSums(tab)
  c=colSums(tab)
  res=r %*% t(c)
  return(res)
}

distance_abs<-function(tab){
  pr=product(tab)
  diff=tab-pr
  diff=diff * diff
  return(sum(diff)*ncol(tab)*nrow(tab))
}

dda_abs<-function(i,j,k,m,ri,cj){ 
  res=0
  if ((i == k) & (j != m)) res=-cj
  if ((i != k) & (m == j)) res=-ri
  if ((i == k) & (j == m)) res=1 - (ri + cj)
  return(res)
}

ddb_abs<-function(tab, c,r, k, m){
  ttab=tab-product(tab)
  s= 0
  for (i in  1:nrow(tab))
    for (j in 1:ncol(tab)) 
      s = s + 2 * ttab[i, j] * dda_abs(i, j, k,m,r[i], c[j])
  return(s*nrow(tab)*ncol(tab)) 
}
  
derivative_distance_abs<-function(tab){
  r=rowSums(tab)
  c=colSums(tab)
  d=matrix(data=NA,nrow=nrow(tab),ncol=ncol(tab))
  for (i in  1:nrow(tab))
    for (j in 1:ncol(tab)) 
      d[i,j]=ddb_abs(tab,c,r,i,j)
  return(d)
}


distance_rel<-function(tab){
  pr=product(tab)
  diff=tab/pr
  diff=diff-1
  diff=diff * diff
  return(sum(diff))
}
dda_rel<-function(i,j,k,m,ri,cj,tab){ 
  res=0
  if ((i == k) & (j != m)) res=-tab[i,j]/(ri*ri*cj)
  if ((i != k) & (m == j)) res=-tab[i,j]/(ri*cj*cj)
  if ((i == k) & (j == m)) {
    z = ri * cj - (ri + cj) * tab[i,j]
    n=ri*ri*cj*cj
    res=z/n
  }
  return(res)
}
ddb_rel<-function(tab, c,r, k, m){
  ttab=tab/product(tab) 
  s= 0
  for (i in  1:nrow(tab))
    for (j in 1:ncol(tab)) 
      s = s + 2 * (ttab[i, j]-1) * dda_rel(i, j, k,m,r[i], c[j],tab)
  return(s) 
}
derivative_distance_rel<-function(tab){
  r=rowSums(tab)
  c=colSums(tab)
  d=matrix(data=NA,nrow=nrow(tab),ncol=ncol(tab))
  for (i in  1:nrow(tab))
    for (j in 1:ncol(tab)) 
      d[i,j]=ddb_rel(tab,c,r,i,j)
  return(d)
}
