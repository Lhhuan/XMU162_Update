






file <- read.delim(x,header=T,stringsAsFactors=F)
file <- list(file[which(file$P.value<0.05),])
down_cluster_homer[[i]] <- as.data.frame(down_cluster_homer[[i]])
clusterTF <- down_cluster_homer[[i]]
str_tmp <- strsplit(clusterTF[,1],'()')
for(j in 1:length(str_tmp)){
	clusterTF$TF[j] <- substring(clusterTF[j,1], 1, which(str_tmp[[j]]=='(')[1]-1)
	clusterTF$motif[j] <- substring(clusterTF[j,1], which(str_tmp[[j]]=='(')[1]+1, which(str_tmp[[j]]==')')[1]-1)
	if(is.na(clusterTF$TF[j])==TRUE){clusterTF$TF[j] <- down_cluster_homer[[i]][j,1]}
	if(is.na(clusterTF$motif[j])==TRUE){clusterTF$motif[j] <- down_cluster_homer[[i]][j,1]}
}
down_cluster_homer[[i]] <- clusterTF
