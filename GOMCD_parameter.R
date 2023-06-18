GMOTE2 <-function(X,target,dupSize=1,G=5,tail_alpha=0.05,alpha=0.5){
  
  # Default Setting
  ncD = ncol(X)
  n_target = table(target)
  classP = names(which.min(n_target))
  P_set = subset(X, target == names(which.min(n_target)))
  N_set = subset(X, target != names(which.min(n_target)))
  P_class = rep(names(which.min(n_target)), nrow(P_set))
  N_class = target[target != names(which.min(n_target))]
  sizeP = nrow(P_set)
  sizeN = nrow(N_set)
  
  # Pre GMM
  tempx<-P_set
  
  size_temp<-nrow(tempx)
  gmm_t<-mclust::mclustBIC(data=tempx)
  gmm_t<-summary(gmm_t)
  gmm_name<-names(gmm_t)[1]
  
  modelNames<-strsplit(gmm_name,",")[[1]][1]
  G1<-as.numeric(strsplit(gmm_name,",")[[1]][2])
  G1=ifelse(G<G1,G,G1)
  gmm<-mclust::Mclust(data=tempx,G=G1,modelNames = modelNames)
  
  while(min(table(gmm$classification))<=ncD){
    G1=G1-1
    gmm<-mclust::Mclust(data=tempx,G=G1,modelNames = modelNames)
    if(G1==1){break}
  }
  
  # calculate the tail probability in MCD-GMM
  
  mahd.mat<-matrix(NA,ncol=G1,nrow=nrow(tempx))
  
  for(i in 1:G1){
    if(length(which(gmm$classification==i))<=2*ncD){
      mahd = mahalanobis(x=tempx, center = gmm$parameters$mean[,i],
                         cov = gmm$parameters$variance$sigma[,,i], tol=1e-100)
    }else{
      mcd=tryCatch({
        robustbase::covMcd(tempx[gmm$classification==i,],alpha = alpha,
                           tolSolve = 1e-100)
      }, warning=function(e){
        center = gmm$parameters$mean[,i]
        cov = gmm$parameters$variance$sigma[,,i]
        list(center=center,cov = cov)
      })
      noise=matrix(0,ncol=ncD,nrow=ncD)
      for(k in 1:10){
        sss=FALSE
        tryCatch({solve(mcd$cov+noise);sss=TRUE},
                 error=function(e){
                   noise1 = mvtnorm::rmvnorm(ncD, rep(0,ncD), diag(ncD)*1e-10)
                   noise1[lower.tri(noise1)] = t(noise1)[lower.tri(noise1)]
                   cc=mcd$cov+noise
                   diag(cc)=abs(diag(cc))
                   diag(noise1)=diag(cc)-diag(mcd$cov)
                   noise<<-noise1
                 })
        if(sss){
          mcd$cov = mcd$cov+noise
          break}
      }
      
      mahd = mahalanobis(x=tempx, center = mcd$center,
                         cov = mcd$cov, tol=1e-100)
    }
    mahd.mat[,i] <- mahd
  }
  
  cutoff = qchisq(1-tail_alpha, df=ncD)
  # Remove Outliers of Pre GMM
  temp_mat<-apply(matrix(1:size_temp),MARGIN = 1,function(x){
    min(mahd.mat[x,])    })
  tempx<-tempx[temp_mat<cutoff,]
  
  # Ready for GMM removed Outliers
  size_temp<-nrow(tempx)
  # dupSize<-round(sizeP/size_temp)*dupSize
  
  gmm_t<-mclust::mclustBIC(data=tempx)
  gmm_t<-summary(gmm_t)
  gmm_name<-names(gmm_t)[1]
  
  modelNames<-strsplit(gmm_name,",")[[1]][1]
  G1<-as.numeric(strsplit(gmm_name,",")[[1]][2])
  G1=ifelse(G<G1,G,G1)
  # Construct GMM
  gmm<-mclust::Mclust(data=tempx,G=G1,modelNames = modelNames)
  while(min(table(gmm$classification))<=ncD){
    G1=G1-1
    gmm<-mclust::Mclust(data=tempx,G=G1,modelNames = modelNames)
    if(G1==1){break}
  }
  
  mcd=list()
  for(i in 1:G1){
    if(length(which(gmm$classification==i))<=2*ncD){
      mcd[[i]]=list()
      center = gmm$parameters$mean[,i]
      cov = gmm$parameters$variance$sigma[,,i]
      
      mcd[[i]]$center = center
      mcd[[i]]$cov = cov
    }else{
      mcd[[i]]=tryCatch({
        robustbase::covMcd(tempx[gmm$classification==i,],alpha = alpha,
                           tolSolve = 1e-80)},
        warning=function(e){
          center = gmm$parameters$mean[,i]
          cov = gmm$parameters$variance$sigma[,,i]
          list(center=center, cov = cov)
        })
      noise=matrix(0,ncol=ncD,nrow=ncD)
      for(k in 1:10){
        sss=FALSE
        tryCatch({solve(mcd[[i]]$cov+noise);sss=TRUE},
                 error=function(e){
                   noise1 = mvtnorm::rmvnorm(ncD, rep(0,ncD), diag(ncD)*1e-10)
                   noise1[lower.tri(noise1)] = t(noise1)[lower.tri(noise1)]
                   cc=mcd[[i]]$cov+noise
                   diag(cc)=abs(diag(cc))
                   diag(noise1)=diag(cc)-diag(mcd[[i]]$cov)
                   noise<<-noise1
                 })
        if(sss){
          mcd[[i]]$cov = mcd[[i]]$cov+noise
          break}
      }
      
    }
  }
  
  pred_N = matrix(NA,ncol=G1,nrow=sizeN)
  for(i in 1:G1){
    pred_N[,i] = mvtnorm::dmvnorm(N_set,mcd[[i]]$center, mcd[[i]]$cov)*gmm$parameters$pro[i]
  }
  Neg_class = apply(pred_N,1,which.max)
  
  # calculate the tail probability of Negative in MCD-GMM
  # calculate sampling weight for each cluster
  nMD = list()
  pMD = list()
  pF=rep(NA,G1)
  Fweight = rep(NA,G1)
  for(i in 1:G1){
    pMD[[i]]=mahalanobis(x=tempx[which(gmm$classification==i),], center = mcd[[i]]$center,
                         cov = mcd[[i]]$cov, tol=1e-80)
    idx = which(pMD[[i]]<cutoff)
    pMD[[i]] = pMD[[i]][idx]
    
    if(length(which(Neg_class==i))!=0){
      nMD[[i]]=mahalanobis(x=N_set[Neg_class==i,], center = mcd[[i]]$center,
                           cov = mcd[[i]]$cov, tol=1e-80)
      idx = which(nMD[[i]]<cutoff)
      if(length(idx)!=0){
        nMD[[i]]=nMD[[i]][idx]
        pF[i] = pf(mean(pMD[[i]])/mean(nMD[[i]]),
                   length(pMD[[i]])*ncD,
                   length(nMD[[i]])*ncD)
        Fweight[i] = mean(pMD[[i]])/mean(nMD[[i]])
      }else{
        nMD[[i]]=NA; pF[i]=0; Fweight[i]=0
      }
    }else{
      nMD[[i]]=NA; pF[i]=0; Fweight[i]=0
    }
  }
  dup_weight=rep(NA,G1)
  for(i in 1:G1){
    # dup_weight[i] = length(pMD[[i]]) / length(nMD[[i]]) * Fweight[i]
    dup_weight[i] = pF[i]
  }
  print(dup_weight)
  dup_weight=ifelse(is.na(dup_weight),0,dup_weight)
  if(sum(dup_weight)==0){
    for(i in 1:G1){dup_weight[i]=1}
  }
  print(G1)
  return(list(mcd=mcd,dup_weight=dup_weight,
              # Fweight=Fweight,
              dupSize=dupSize, clustering = gmm$classification))}