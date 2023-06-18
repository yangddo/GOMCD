
sampling = function(X=X, target=target, result, beta){
  ncD = ncol(X)
  n_target = table(target)
  classP = names(which.min(n_target))
  P_set = subset(X, target == names(which.min(n_target)))
  N_set = subset(X, target != names(which.min(n_target)))
  P_class = rep(names(which.min(n_target)), nrow(P_set))
  N_class = target[target != names(which.min(n_target))]
  sizeP = nrow(P_set)
  sizeN = nrow(N_set)
  
  result$dup_weight=ifelse(is.na(result$dup_weight),0,result$dup_weight)
  result$dup_weight=result$dup_weight/sum(result$dup_weight)
  sample_ratio = (result$dup_weight*beta +
                    as.numeric(prop.table(table(result$clustering)))*(1-beta))*
    sizeP*result$dupSize
  sample_ratio = round(sample_ratio)
  
  # Generate artificial instances for each Mixture
  res<-NULL
  G1 = length(result$mcd)
  for(i in 1:G1){
    synNum <- sample_ratio[i]
    if(is.na(synNum)){synNum=0}
    if(synNum!=0){
      temp_result <- mvtnorm::rmvnorm(synNum,result$mcd[[i]]$center,result$mcd[[i]]$cov,method="svd")
      res<-rbind(res,temp_result)
    }
    
  }
  
  syn_dat<-data.frame(res)
  
  P_set[, ncD + 1] = P_class
  N_set[, ncD + 1] = N_class
  colnames(P_set) = c(colnames(X), "class")
  colnames(N_set) = c(colnames(X), "class")
  syn_dat[, ncD + 1] = rep(names(which.min(n_target)), nrow(syn_dat))
  colnames(syn_dat) = c(colnames(X), "class")
  NewD = rbind(P_set, syn_dat, N_set)
  rownames(NewD) = NULL
  
  
  D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set, 
                  orig_P = P_set, 
                  # dup_size = dupSize, 
                  # alpha = alpha,tail_alpha=tail_alpha,
                  method = "MCDMOTE")
  class(D_result) = "gen_data"
  print(paste(G1," components",sep=""))
  # print(paste("MCDMOTE use ",length(which(dup_weight>0))," components",sep=""))
  return(D_result)
}