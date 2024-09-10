## R script to process Scer alleles from VCF files containing freebayes RO and AO

## Feb 2022
## R code from Yue Hue:
############################################################################################


#figure settings
fig.margin.defaults <- c(1.02, 0.82, 0.82, 0.42)
fig.margin.width <- sum(fig.margin.defaults[c(2,4)])
fig.margin.height <- sum(fig.margin.defaults[c(1,3)])
page.margins <- 2.0
page.width <- 11.693
plot.width <- page.width-page.margins-fig.margin.width
plot.height <- plot.width * 0.7071
fig.width <- plot.width + fig.margin.width
fig.height <- plot.height + fig.margin.height


#founder scer
Scer <- read.csv("Input/founder/HS15_S66.csv", stringsAsFactors=FALSE,row.names =    ## This should be the alternate founder ie HS15 (HS14 is the ref)
                   1)
scer_chr <- unique(Scer$CHROM)
Allele_scer <-list()
for (i in 1:17) {
  Allele_scer[[i]] <- Scer[which(Scer$CHROM == scer_chr[i]),]
}
names(Allele_scer) <- scer_chr

names(Allele_scer) ## check names



#founder skud
skud <- read.csv("Input/founder/HS11_S62.csv", stringsAsFactors=FALSE, row.names =1)  ## This should be the alternate founder ie HS11 (HS10 is the ref)
skud_chr <- unique(skud$CHROM)

Allele_skud <-list()
for (i in 1:length(skud_chr)) {
  Allele_skud[[i]] <- skud[which(skud$CHROM == skud_chr[i]),]
}
names(Allele_skud) <- skud_chr[1:length(skud_chr)]

names(Allele_skud)

#Fit
HY3_list <- list.files('Input/pool/')
#HY3_ploidy2 <- list()
HY3 <- list()
for (i in 1:length(HY3_list)) {
    HY3[[i]] <- read.csv(paste0('Input/pool/',HY3_list[i]),stringsAsFactors = FALSE,row.names = 1)
}
dir_name <- unlist(strsplit(HY3_list,'.csv'))


## check
head( HY3[[i]])
tail( HY3[[i]])


names(HY3) <- dir_name
HY3_scer <- list()
HY3_skud <- list()


HY3_scer <-lapply(HY3,function(x){x[grep('chr',x$CHROM),]})
HY3_skud <- lapply(HY3, function(x){x[grep('Sk',x$CHROM),]})


for ( i in 1:2){ HY3_scer[[i]]$CHROM <- as.factor(HY3_scer[[i]]$CHROM)
#levels(HY3_scer[[i]]$CHROM)<- paste0('chr',c(1,2,3,4,9,5,6,7,8,10,11,12,13,14,15,16))
}
chr <- unique(Scer$CHROM)

Inter_chrom <- list()
for(i in 1:16){
  Inter_chrom[[i]] <- Reduce(intersect,lapply(HY3_scer,function(x){x[x$CHROM==chr[i],2]}))
}
names(Inter_chrom) <-  chr[1:16]

#allele scer
H_scer <- list()

## Choose the bi-allele markers which inheritated from founders.

for (i in 1:2) {
  hybrid_scer<- HY3_scer[[i]]
  #ALT have 1 base different with reference
  hybrid_scer_alt1 <- hybrid_scer[which(nchar(hybrid_scer$ALT)==nchar(hybrid_scer$REF)),]   
  AO <- as.numeric(hybrid_scer_alt1$AO)
  RO <- as.numeric(hybrid_scer_alt1$RO)
  filter_index <- which((AO+RO>10)&(AO+RO<1000) == TRUE)
  hybrid_scer <- data.frame()
  H_scer[[i]] <- hybrid_scer_alt1[filter_index,]
}


HY3scer_af <- lapply(H_scer, function(x){as.numeric(x[,5])/(as.numeric(x[,5])+as.numeric(x[,6]))})
HY3scer_reads <- lapply(H_scer, function(x){as.numeric(x[,5])+as.numeric(x[,6])})
for (i in 1:2) { 
  HY3scer_af[[i]] <- cbind(H_scer[[i]],HY3scer_af[[i]],HY3scer_reads[[i]])
  HY3scer_af[[i]] <- as.data.frame(HY3scer_af[[i]])
  colnames(HY3scer_af[[i]])[7] <- 'Allele_frequency'
  colnames(HY3scer_af[[i]])[8] <- 'Total_reads'
}

HY3.scer <- list()
for (j in 1:2) {
  HY3.scer.fit <- list()
  for (i in 1:16) {
    sub_HY3 <- HY3scer_af[[j]][HY3scer_af[[j]]$CHROM == chr[i], ]
    sub_scer <- Scer[Scer$CHROM == chr[i],]
    HY3.scer.fit[[i]] <- sub_HY3[which((sub_HY3[,2] %in% sub_scer[,2])==TRUE),]
  }
  HY3.scer[[j]] <- HY3.scer.fit
  names(HY3.scer[[j]]) <- chr[1:16]
}
names(HY3.scer) <- dir_name


pdf(file = 'HY3_T_allele_scer_plot.pdf',width = fig.width, height = fig.height, onefile = TRUE)
plot.device <- dev.cur()
for (i in 1:16) {
  plot(HY3.scer[[1]][[i]]$POS, HY3.scer[[1]][[i]]$Allele_frequency, type = 'l',xlim = range(c(0,HY3.scer[[7]][[i]]$POS)),ylim=range(0,1.2),axes = FALSE, ann = FALSE, col = 'green')
  lines(HY3.scer[[2]][[i]]$POS, HY3.scer[[2]][[i]]$Allele_frequency, col = 'red' )
  
  title(main=chr[i], col.main="black", font.main=4)
  axis(1, at = HY3.scer[[1]][[i]]$POS, labels = HY3.scer[[1]][[i]]$POS )
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0))
  title(xlab="Position", col.lab= 'black')
  title(ylab = 'Allele frequency',col.lab = 'black')
  legend('topright', legend = c('High','Low'),col=c('green','red'),lty=1)
}
dev.off()


for (i in 1:16) {
  for (j in 1:2) {
    write.table(HY3.scer[[j]][[i]][,c(2,5,6)],file = paste0('Output_scer/',dir_name[j],'/',chr[i],'.txt'),quote = FALSE,row.names = FALSE,col.names = FALSE)
  }
  
}


pdf(file = 'HY3_MAL_allele_scer_plot.pdf',width = fig.width, height = fig.height, onefile = TRUE)
plot.device <- dev.cur()
for (i in 1:16) {
  plot(HY3.scer[[1]][[i]]$POS, HY3.scer[[1]][[i]]$Allele_frequency, type = 'l',xlim = range(c(0,HY3.scer[[1]][[i]]$POS)),ylim=range(0,1.2),axes = FALSE, ann = FALSE, col = 'green')
  lines(HY3.scer[[2]][[i]]$POS, HY3.scer[[2]][[i]]$Allele_frequency, col = 'red' )
  title(main=chr[i], col.main="black", font.main=4)
  axis(1, at = HY3.scer[[1]][[i]]$POS, labels = HY3.scer[[1]][[i]]$POS )
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0))
  title(xlab="Position", col.lab= 'black')
  title(ylab = 'Allele frequency',col.lab = 'black')
  legend('topright', legend = c('High','Low'),col=c('green','red'),lty=1)
}
dev.off()

