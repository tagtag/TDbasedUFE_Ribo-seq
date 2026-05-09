x<-read.csv("GSE291652_RNA_seq_raw_count.tsv.gz",sep="\t")
x1<-read.csv("GSE291653_Ribo_seq_raw_count.tsv.gz",sep="\t")
require(biomaRt)
require(readxl)
x2 <- read_excel("13059_2025_3800_MOESM4_ESM.xlsx",sheet=4)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_genes <- unlist(x2[,1])
result <- getBM(attributes = c("mgi_symbol", "ensembl_gene_id"), 
                filters = "mgi_symbol", 
                values = mouse_genes, 
                mart = mart)
index <- match(result[match(unlist(x2[,1]),result[,1]),2],x[,1])
Z <- array(NA,c(dim(x)[1],3,6,3))
X<-data.matrix(x[,c(-1,-20)])
dim(X)<-c(18175,3,6)
Z[,,,1] <- X
X<-data.matrix(x1[,c(-1,-20)])
dim(X)<-c(18175,3,6)
Z[,,,2] <- X
 X <- matrix(NA,dim(x2)[1],18)
 X[,c(1:13,15:18)]<-data.matrix(x2[,-1])
X[,14] <- X[,13]
dim(X)<-c(dim(X)[1],3,6)
Z[index[!is.na(index)],,,3] <- X[!is.na(index),,]
Z[is.na(Z)]<-0

Z<-aperm(Z,c(1,3,2,4))
Z[,,,1]<-log10(Z[,,,1]+1)
Z[,,,2]<-log10(Z[,,,2]+1)
Z<-apply(Z,2:4,scale)
Z[is.na(Z)]<-0
require(rTensor)
HOSVD<-hosvd(as.tensor(Z),c(10,6,3,3))

pdf(file="Uk.pdf",width=5,height=15)
par(mfrow=c(3,1))
plot(HOSVD$U[[4]][,1],type="h",ylim=c(-0.7,0.0),cex.lab=1.5,cex.axis=2,ylab="U1k",xlab="k");abline(0,0,col=2,lty=2,lwd=2)
plot(HOSVD$U[[4]][,2],type="h",cex.lab=1.5,cex.axis=2,ylab="U2k",xlab="k");abline(0,0,col=2,lty=2,lwd=2)
plot(HOSVD$U[[4]][,3],type="h",cex.lab=1.5,cex.axis=2,ylab="U3k",xlab="k");abline(0,0,col=2,lty=2,lwd=2)
par(mfrow=c(1,1))
dev.off()

pdf(file="image.pdf")
image(1:10,1:6,log(abs(HOSVD$Z@data[,,1,2])),xlab="l1",ylab="l2",cex.lab=1.5,cex.axis=2)
dev.off()

pdf(file="Uj.pdf")
par(mfrow=c(1,2))
barplot(HOSVD$U[[2]][,1],type="h",ylim=c(-0.5,0.0),cex.lab=1.5,cex.axis=2,ylab="U1j",xlab="j");abline(0,0,col=2,lty=2,lwd=2)
barplot(HOSVD$U[[2]][,2],type="h",cex.lab=1.5,cex.axis=2,ylab="U2j",xlab="j");abline(0,0,col=2,lty=2,lwd=2)
par(mfrow=c(1,1))
dev.off()

k0<-6
th <- function(sd){
  P2<- pchisq(((HOSVD$U[[1]][,k0]-mean(HOSVD$U[[1]][,k0]))/sd)^2,1,lower.tail=F)
  hc<- hist(1-P2,breaks=100,plot=F)
  return(sd(hc$count[1:sum(hc$mid<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
  #hc<- hist(1-P2[p.adjust(P2,"BH")>0.01],breaks=100,plot=F)
  #return(sd(hc$count))
}

pdf(file="optimize.pdf")
par(mfrow=c(1,2))
#for (k0 in c(3:4,6))
k0<-6
  cat(k0," ")
sd <- optim(0.003,th)$par
P1<- pchisq(((HOSVD$U[[1]][,k0]-mean(HOSVD$U[[1]][,k0]))/sd)^2,1,lower.tail=F)
aa <- seq(0.5*sd,2*sd,by=0.05*sd)
bb<-apply(matrix(seq(0.5*sd,2*sd,by=0.05*sd),ncol=1),1,th)
plot(aa,bb,xlab="sigma_l",ylab="sigma_h",type="o",cex.lab=2,cex.axis=2)
arrows(sd,max(bb),sd,min(bb),col=2)
hist(1-P1,breaks=100,xlab="1-Pi",cex.lab=2,cex.axis=2)
par(mfrow=c(1,1))
dev.off()
P1<- pchisq(((HOSVD$U[[1]][,k0]-mean(HOSVD$U[[1]][,k0]))/sd)^2,1,lower.tail=F)
table(p.adjust(P1,"BH")<0.01)
#FALSE  TRUE 
#16394  1781 

genes <- data.frame(x[p.adjust(P1,"BH")<0.01,20])
write.table(file="genes_PC6_prot.csv",genes,sep="\t",row.names=F,col.names=F,quote=F)

pdf(file="image2.pdf")
image(1:10,1:6,log(abs(HOSVD$Z@data[,,1,3])),xlab="l1",ylab="l2",cex.lab=1.5,cex.axis=2)
dev.off()

k0<-3
th <- function(sd){
  P2<- pchisq(((HOSVD$U[[1]][,k0]-mean(HOSVD$U[[1]][,k0]))/sd)^2,1,lower.tail=F)
  hc<- hist(1-P2,breaks=100,plot=F)
  return(sd(hc$count[1:sum(hc$mid<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
  #hc<- hist(1-P2[p.adjust(P2,"BH")>0.01],breaks=100,plot=F)
  #return(sd(hc$count))
}

pdf(file="optimize2.pdf")
par(mfrow=c(1,2))
k0<-3
  cat(k0," ")
sd <- optim(0.002,th)$par
P1<- pchisq(((HOSVD$U[[1]][,k0]-mean(HOSVD$U[[1]][,k0]))/sd)^2,1,lower.tail=F)
#hist(1-P1,breaks=100)
aa <- seq(0.5*sd,2*sd,by=0.05*sd)
bb<-apply(matrix(seq(0.5*sd,2*sd,by=0.05*sd),ncol=1),1,th)
plot(aa,bb,xlab="sigma_l",ylab="sigma_h",type="o",cex.lab=2,cex.axis=2)
arrows(sd,max(bb),sd,min(bb),col=2)
hist(1-P1,breaks=100,xlab="1-Pi",cex.lab=2,cex.axis=2)
par(mfrow=c(1,1))
dev.off()

P1<- pchisq(((HOSVD$U[[1]][,k0]-mean(HOSVD$U[[1]][,k0]))/sd)^2,1,lower.tail=F)
table(p.adjust(P1,"BH")<0.01)

#FALSE  TRUE 
#17948   227 

genes <- data.frame(x[p.adjust(P1,"BH")<0.01,20])
write.table(file="genes_PC3_prot.csv",genes,sep="\t",row.names=F,col.names=F,quote=F)


