library(FITSio)
# library(DescTools)    # used for Gini function

# define lm-based rescaling algorithm
lm.scale = function(source,target,breaks=seq(.15,.85,.05)){
  return(unname(coef(lm(quantile(target,breaks,na.rm=T)~quantile(source,breaks,na.rm=T)))))
}

# define loess diff peak extraction algorithm
loess.peak = function(y,span=0.2,flip=TRUE){
  x=1:length(y)
  return((2*(.5-flip))*(y-predict(loess(y~x,span=span),data.frame(x=x))))
}

# define low.pass function
low.pass = function(flux,thresh=200,strength=0.03,debug=FALSE){
  f.raw = Re(fft(fft(flux)*exp(-exp(((1+strength)^2-1)*((1:length(flux))-thresh))),inverse=T))
  ab = lm.scale(f.raw,flux)
  return(f.raw*ab[2]+ab[1])
}

# define normalized peak area function
norm.peak.area = function(vector,index){
  fallingOrZeroBound = function(vec){
    return(which( c(F,diff(vec)>0) | vec==0 )[1]-1)
  }
  n = length(vector)
  return(sum(vector[c(
    max(1,index-fallingOrZeroBound(rev(vector[1:(index-1)])),na.rm=T):
      min(n,index+fallingOrZeroBound(vector[(index+1):n]),na.rm=T)
  )])/sum(vector))
}

# # create template to use
# cb58                = readFrameFromFITS("cB58_Lyman_break.fit")
# names(cb58)         = tolower(names(cb58))
# template            = data.frame(flux.pass=low.pass(cb58$flux))
# template$flux.peak  = loess.peak(template$flux.pass)
# save(template, file ="cb58_template.Rdata")

# load template and save # rows for use later
load("cb58_template.Rdata")
N = nrow(template)

# fraction of max-conv-peak used as noise threshold
C_thresh = 1/6

# fraction of template used for scaling (to eliminate quasar spectra with broad high-wavelength emissions)
N_frac   = floor(N*2/3)

# K_thresh is K-stat of 1353 (noisy cB58); above this, do not penalize
K_thresh = 0.8326

ScoreSpec = function(spec.file){
  
  # read in spectrum and check column names
  spec = readFrameFromFITS(spec.file)
  names(spec) = tolower(names(spec))
  spec = spec[,"flux",drop=F]
  
  # smooth spectrum and enhance peaks
  spec$flux.pass = low.pass(spec$flux)
  spec$flux.peak = loess.peak(spec$flux.pass)
  
  # convolve with template, and get offset of max
  conv       = convolve(spec$flux.peak,template$flux.peak,type="filter")
  conv.max   = which.max(conv)
  
  # transform to emphasize peaks, and find normalized area of peak at max-conv offset
  conv.trans = pmax(0,abs(conv)-conv[conv.max]*C_thresh)^2
  A          = norm.peak.area(conv.trans,conv.max)
  
  # subset spectrum at this offset
  newSpec = spec$flux.pass[conv.max:(conv.max+N-1)]
  
  # rescale fromt 2/3 of template to spectrum subset
  # (this attempts to eliminate quasar spectra with broad high-wavelength emissions)
  coefs   = lm.scale(template$flux.pass[1:N_frac],newSpec[1:N_frac])
  newTemp = template$flux.pass*coefs[2]+coefs[1]
  
  # find 1-Kolmogorov-Smirnov statistic to characterize quality of match
  K = unname(1-ks.test(newSpec,newTemp)$statistic)
  
  # generate scores and return
  scores = c(A,K,A*pmin(1,K+(1-K_thresh))^2)
  names(scores) = c("A","K","Composite")
  scores
}

files = list.files("data/",full.names=T)
res = as.data.frame(t(sapply(files,ScoreSpec)))
res = res[order(res$Composite,decreasing=T),]
rownames(res) = gsub(".*/|(-\\d{5}|_temp).*","",rownames(res))
res$names = factor(rownames(res),levels=rev(rownames(res)),ordered=T)
res$types = "Other"
res[grep("spec-1353|8oclockArc",res$names),]$types = "Lyman-break galaxy"
res[grep("spec-(5302|5324|5328)",res$names),]$types = "Quasar"
res[grep("spec-6064",res$names),]$types = "Carbon star"
res$types = factor(res$types,ordered=T,
                   levels=c("Lyman-break galaxy","Quasar","Carbon star","Other"))

# visualize result
library(ggplot)
ggplot(res[1:30,],aes(x=names,y=Composite,fill=types)) + geom_col() + 
  coord_flip() + labs(x=NULL,fill="Type of observation") + 
  theme(legend.position=c(.98,.03),legend.justification=c(1,0)) + 
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1),expand=c(0,0)) + 
  scale_fill_manual(values=rev(c('black','#a1dab4','#41b6c4','#225ea8'))) + 
  ylab(paste0("Normalized convolution peak area  \u00d7  scaled ",
  "1-(K-S statistic of partially-quantile-matched spectra)")) + 
  ggtitle("Top 30 spectra scores in sample of 100")
ggsave("res30.png",width=9,height=6,dpi=600)

# save result
write.table(capture.output(res),file="res2.txt",quote=F,row.names=F,col.names=F)
