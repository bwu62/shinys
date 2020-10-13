library(FITSio)
library(DescTools)

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

# # create template to use
# cb58                = readFrameFromFITS("cB58_Lyman_break.fit")
# names(cb58)         = tolower(names(cb58))
# template            = data.frame(flux.pass=low.pass(cb58$flux))
# template$flux.peak  = loess.peak(template$flux.pass)
# save(template, file ="cb58_template.Rdata")

# load template and save # rows for use later
load("cb58_template.Rdata")
N = nrow(template)

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
  
  # convolve with template, get offset of max, and transform for Gini later
  conv       = convolve(spec$flux.peak,template$flux.peak,type="filter")
  conv.max   = which.max(conv)
  M          = conv[conv.max]
  conv.trans = pmax(0,abs(pmax(conv,-M))-M/8)^2
  
  # subset spectrum at this offset
  newSpec = spec$flux.pass[conv.max:(conv.max+N-1)]
  
  # rescale fromt 2/3 of template to spectrum subset
  # (this attempts to eliminate quasar spectra with broad high-wavelength emissions)
  coefs   = lm.scale(template$flux.pass[1:N_frac],newSpec[1:N_frac])
  newTemp = template$flux.pass*coefs[2]+coefs[1]
  
  # generate scores and return
  scores = c(
    G <- Gini(conv.trans),
    K <- unname(1-ks.test(newSpec,newTemp)$statistic),
    
    # 0.157 selected for N_frac=2/3 so 1353 (noisy cB58) not penalized for K
    # for N_frac=3/4, use 0.168
    G*pmin(1,K+(1-K_thresh))^2
  )
  names(scores) = c("G","K","Composite")
  scores
}

files = list.files("data/",full.names=T)
res = as.data.frame(t(sapply(files,ScoreSpec)))
res = res[order(res$Composite,decreasing=T),]

# save result
write.table(capture.output(res),file="res.txt",quote=F,row.names=F,col.names=F)
