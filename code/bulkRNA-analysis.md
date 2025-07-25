# survival analysis
### for example, POSTN, SPP1, and POSTN+SPP1, use `median()` to select cut point.
```
####SPP1+ & POSTN+ survival#####
setwd("/media/gigabyte/Data/DTH/PDAC/RNA/")
# TCGA gene survival
rt<-read.table("./gsva/cli/TCGA-PAAD_cli.txt",sep = '\t',header = 1)
exp<-read.table("./gsva/exp/TCGA-PAAD_tumor_symbol.txt")
colnames(exp)<-gsub('.','-',colnames(exp),fixed = T)
nexp<-as.data.frame(t(exp))
#select survival gene
A<-'POSTN'
B<-'SPP1'
deg.data<-nexp[,c('POSTN','SPP1')]
deg.data$sample<-rownames(deg.data)
deg.data<-merge(rt,deg.data,by = "sample")
colnames(deg.data)[4:5]<-c("A","B")
deg.data$group<-NA
rt.cut<-surv_cutpoint(deg.data,time = "survival.time",event = "survival.status",
                      variables = c("A","B"))
res.cut<-surv_categorize(rt.cut)
#POSTN high && SPP1 high survival
res.cut$C<-NA
res.cut$C[which((res.cut$A=="high" )&(res.cut$B=="high"))] =paste(A,"high",B,"high")
res.cut$C[which((res.cut$A=="low" )&(res.cut$B=="low"))] =paste(A,"low",B,"low")
## survival plot
dir.create("./subtype/Fibro_Macro/Fibro_Macro_subtype_survival/gene")
setwd("./subtype/Fibro_Macro/Fibro_Macro_subtype_survival/gene")
cols<-c("#DC0000CC","#3C5488CC")
#POSTN high + SPP1 high / POSTN low + SPP1 low
{
    fit<-survdiff(as.formula(paste0('Surv(survival.time,survival.status)~',"C")),data = res.cut)
    pValue<-1-pchisq(fit$chisq,df=1)
    pValue<-round(pValue,3)
    if(pValue < 0.001){pValue<-"p < 0.001"}
    fit<-survfit(Surv(survival.time,as.numeric(survival.status))~C,data=res.cut)
    p1<-ggsurvplot(fit,
                   data = res.cut,
                   legend.labs=c("High","Low"),
                   palette = cols,#c("red","blue"),
                   pval = pValue,
                   pval.size=10,
                   risk.table = T,
                   surv.median.line = "hv", 
                   #palette=c("red", "blue"),  
                   legend.title="Exp",
                   title=paste("POSTN+_SPP1+","Overall survival",sep = ' '), #title
                   ylab="Cumulative survival (percentage)",xlab = " Time (Months)", 
                   censor.shape = 124,censor.size = 2,conf.int = FALSE, 
                   break.x.by = 12)
    pdf(file = paste("./","POSTN_SPP1-survival.pdf",sep =''),width = 5,height = 5,onefile = F)
    print(p1)
    dev.off()
}
```
# GSVA of cell types
```
library(GSVA)
#makergenesets
genesets = readLines("/media/gigabyte/Data/DTH/PDAC/RNA/celltype.gmt")
res <- strsplit(genesets, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
genesets <- lapply(res, "[", -c(1:2))
#
expr<-list()
#k=list.files(path="./exp/")[2]
for (i in list.files(path="./gsva/exp/",pattern = ".txt")) {
    expr[[i]] <-read.table(paste("./gsva/exp/",i,sep = ''),sep = '\t')
}
#
all<-matrix(,,6)
for (k in list.files(path="./exp/")) {
gsva.res<-gsva(as.matrix(expr[[k]]),genesets,method="gsva",min.sz > 0)
write.table(gsva.res,paste(strsplit(k,'_')[[1]][1],'_gsva.txt',sep = ''),sep='\t',quote = F)
n=1
cox.res<-matrix(0,length(genesets)*length(genesets),6)
colnames(cox.res)<-c("data","celltp1","celltp2","p","R","color")
for (i in 1:dim(gsva.res)[1]) {
    for (j in 1:dim(gsva.res)[1]) {
        cox.res[n,1]<-strsplit(k,'_')[[1]][1]
        cox.res[n,2]<-names(genesets)[i]
        cox.res[n,3]<-names(genesets)[j]
        a<-cor.test(gsva.res[i,],gsva.res[j,],method = "spearman")
        cox.res[n,4]<-a$p.value
        cox.res[n,5]<-a$estimate
        n=n+1
    }
}    
cox.res[as.numeric(cox.res[,4])<0.05 & 0.3<=as.numeric(cox.res[,5]),6]<-"red"
cox.res[as.numeric(cox.res[,4])<0.05 & -0.3>=as.numeric(cox.res[,5]),6]<-"blue"
cox.res[as.numeric(cox.res[,4])>0.05 ,6]<-"grey"
cox.res[-0.3<as.numeric(cox.res[,5]) & 0.3>as.numeric(cox.res[,5]),6]<-"grey"
all<-rbind(all,cox.res)
write.table(cox.res,paste(strsplit(k,'_')[[1]][1],'_cox.txt',sep = ''),sep='\t',quote = F)
}
all<-all[-1,]
write.table(all,'all_cox.txt',sep='\t',quote = F,row.names = F)
all<-as.data.frame(all)

bulkresult <-  read.table('analysis/Bulk-analysis/bulk_cellytype_all_cox_delete_F_M.txt',sep = "\t",check.names = F,header = T)
bulkresult <- bulkresult[,-1]
bulkresult <- bulkresult[-grep("0",bulkresult$celltp1),]

result<-matrix(0,1*81*3,5)
colnames(result)<-c("Celltype1","Celltype2","Celltype1_celltype2","correlation","count")


cell1<- unique(bulkresult$celltp1)
cell2<- unique(bulkresult$celltp2)


n=1
m=2
l=3
for (i in 1:length(cell1)) {
    tmp <- bulkresult[grep(cell1[i],bulkresult$celltp1),]
     for (j in 1:length(cell2)) {
     tmp2 <- tmp[grep(cell2[j],tmp$celltp2),]
     
     result[n,1] <- cell1[i]
     result[n,2] <- cell2[j]
     result[n,3] <- paste(cell1[i],"_",cell2[j],sep = "")
     result[n,4] <- "red"
     result[n,5] <-sum(tmp2 =="red")
     result[m,1] <- cell1[i]
     result[m,2] <- cell2[j]
     result[m,3] <- paste(cell1[i],"_",cell2[j],sep = "")
     result[m,4] <- "grey"
     result[m,5] <-sum(tmp2 =="grey")
     result[l,1] <- cell1[i]
     result[l,2] <- cell2[j]
     result[l,3] <- paste(cell1[i],"_",cell2[j],sep = "")
     result[l,4] <- "blue"
     result[l,5] <-sum(tmp2 =="blue")
     
      n=n+3
      m=n+1
      l=n+2
     
  }
 
}

result <- as.data.frame( result)
result$count <- as.numeric(result$count)

ggplot(data =result, aes(x = factor(1), y = count, fill = factor(correlation))) +
  geom_bar(stat="identity",position = "stack", color = NA) +
  scale_y_continuous(limits = c(0,12),expand=c(0,0)) +
  coord_polar("y") +
  facet_grid(Celltype2~Celltype1) +
  theme(axis.text=element_blank(),axis.title = element_blank(), 
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle =0,hjust=0,color="black",size=8),strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,angle = 90,vjust = 0),legend.text = element_text(size=10),
        legend.position = "bottom",panel.spacing  = unit(0.01, "lines")) +
  scale_fill_manual(limits=c("blue","grey","red"),values=c( "#377EB8","lightgray","#E41A1C"))


```




