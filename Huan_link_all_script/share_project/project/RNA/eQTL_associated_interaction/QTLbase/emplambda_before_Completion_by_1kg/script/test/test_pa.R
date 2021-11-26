library(doParallel)
library(foreach)

test=function(x){
  while(x>0){
    x=x+1
  }
  return(x)
}


cores=detectCores() #检查可用核数
cl <- makeCluster(10)
registerDoParallel(cl)
result1=foreach(i=1:100) %dopar% test(i)
stopImplicitCluster()