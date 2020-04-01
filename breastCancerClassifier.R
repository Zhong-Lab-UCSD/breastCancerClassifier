library('zoo')
library('e1071')
library('DESeq')
library('ggplot2')

# SVM related functions ####
roc<-function(dis) # ROC curve calculation
{
  temp_dis=dis[order(dis,decreasing = T)]
  return(sapply(temp_dis,function(x){n=dis[1:8];c=dis[9:56];tp=sum(n>=x);fp=sum(c>=x);return(c(tp/8,fp/48))}))
  #tpr=tp/(tp+fn)  fpr=fp/(fp+tn)
}
AUC <- function(x, y) # AUC calculation
{
  sum(diff(x)*rollmean(y,2))
}
gg_color_hue <- function(n) # generate color in plot
{
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# load the TPM file ####
tpm <- read.table('tpm_96_nodup.txt', row.names = 1, stringsAsFactors = F)
colnames(tpm) <- paste0('C',1:96)
tpm_nz = tpm[rowSums(tpm)>0,] # non-zero genes
readcounts <- read.table('readcounts_96_nodup.txt',row.names = 1, stringsAsFactors = F)
readcounts <- readcounts[1:60675,]
colnames(readcounts) <- paste0('C',1:96)

# load meta data file ####
patient_info <- read.csv('patient_info.csv',stringsAsFactors = F) # recurrence status file

preselectedList <- as.vector(unlist(read.table('preselectedList'))) # breast cancer biomarkers

# calculate differential expressed genes among the marker list ####
condition <- as.factor(patient_info$recurStatus)
cds <- newCountDataSet(readcounts[rownames(readcounts) %in% preselectedList,],condition)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "R", "N")

# initialize parameters ####
trainset <- matrix(0,15,3) # initialize trainset
cost_value <- 5 # svm paramter cost
averageAUC <- rep(0,20) # average AUC for all iterations
geneNum <- rep(0,20) # number of biomarkers (genes) of corresponding threshold

# to increase performance, this could be written into apply function
for(j in 1:length(seq(0.05,1,0.05)))
{
  th <- seq(0.05,1,0.05)[j] #threshold
  
  difGene <- res[res$pval<th & !is.na(res$pval),1] # differential gene selected from the threshold
  difGene <- rownames(tpm)[rownames(tpm) %in% difGene]

  temp <- t(tpm_nz[rownames(tpm_nz) %in% difGene,])
  temp <- as.data.frame(temp)
  temp$status[patient_info$recurStatus=='R'] <- 1
  temp$status[patient_info$recurStatus=='N'] <- 0
  temp$status <- as.factor(temp$status)
  
  index <- 1:nrow(temp)
  auc <- rep(0,100)
  set.seed(1234)
  for(i in 1:100) # for each iteration randomly generate the 20 recurrence samples and 20 non-recurrence samples to balance the traning set
  {
    c <- 1
    while(c==1 | length(which(colSums(trainset[,1:(ncol(trainset)-1)])==0))>0)
    {
      testindex <- c(sample(1:28,8),sample(29:96,48))
      testset <- temp[testindex,,drop=F]
      trainset <- temp[-testindex,]
      c <- c + 1
    }
    svm.model <- svm(status ~ .,data=trainset,cost=cost_value,kernel='sigmoid')
    svm.pred <- predict(svm.model,testset[,-ncol(testset)],decision.value=T)
    roc_curve <- roc(attr(svm.pred,'decision.values'))
    auc[i] <- AUC(x=roc_curve[2,],y=roc_curve[1,])
    if(i == 4 & j == 11) # when j == 1, the average AUC is maximum
      df_roc4 <- as.data.frame(t(roc_curve))
    if(i == 5 & j == 11)
      df_roc5 <- as.data.frame(t(roc_curve))
    if(i == 6 & j == 11)
      df_roc6 <- as.data.frame(t(roc_curve))
  }
  averageAUC[j] <- mean(auc)
  geneNum[j] <- length(difGene)
}

# plot the number of biomarker against average AUC plot ####

df <- data.frame(geneNum = geneNum, averageAUC = averageAUC)

ggplot(df, aes(x=geneNum, y=averageAUC)) +
  geom_point(shape=16, color = 'red',size=1.2) + 
  scale_y_continuous(limits = 0:1) +
  theme(plot.background = element_rect(fill = 'white', colour = 'grey',linetype = 'solid'),
        panel.background = element_rect(fill = 'white',color = 'grey',size=0.75,linetype='solid'),
        panel.grid.major = element_line(colour = "grey",size=0.2),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.75)) +
  labs(x = 'Number of genes', y = 'Average AUC',fill = '') + 
  theme(axis.text.x = element_text(size=8,hjust = 1),axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=11),axis.title.x = element_text(size=11),
        legend.text = element_text(size = 8)) 
ggsave('averageAUC.pdf')

# plot roc curve from randomly generated when threshold is 0.55 ####
cols <- gg_color_hue(3)

g1 <- ggplot(df_roc4, aes(x=V2, y=V1)) + geom_line(color = cols[1]) +
  labs(x = '1 - Specificity', y = 'Sensitivity',fill = '') +
  theme(plot.background = element_rect(fill = 'white', colour = 'grey',linetype = 'solid'),
        panel.background = element_rect(fill = 'white',color = 'grey',size=0.75,linetype='solid'),
        panel.grid.major = element_line(colour = "grey",size=0.2),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.75)) +
  theme(axis.text.x = element_text(size=8,hjust = 1),axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=11),axis.title.x = element_text(size=11),
        legend.text = element_text(size = 8)) +
  scale_y_continuous(limits = 0:1) +
  scale_x_continuous(limits = 0:1)
g1 <- g1 + geom_line(data = df_roc5,color = cols[2])
g1 <- g1 + geom_line(data = df_roc6,color = cols[3])
g1 <- g1 + annotate('text',label = paste0('Average AUC = ',format(round(averageAUC[11], 2), nsmall = 3)), 
                  x = 0.65, y=0.2, size = 3)
g1
gsave('ROCplot.pdf')

