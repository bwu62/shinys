library(FITSio)

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

# load template and save # rows for use later
load("cb58_template.Rdata")
N = nrow(template)

# fraction of max-conv-peak used as noise threshold
C_thresh = 1/6

# fraction of template used for scaling (to eliminate quasars with broad high-wavelength emissions)
N_frac   = floor(N*2/3)

# K_thresh is K-stat of 1353 (noisy cB58); above this, do not penalize
K_thresh = 0.1673544246

# define K-stat rescaling function
K_scale = function(x){1-pmin(pmax(x-K_thresh,0)/(0.8-K_thresh),1)}

ScoreSpec = function(spec.file){
  
  scores = setNames(rep(NA,4),c("Area.peak","KS.raw","KS.trans","Combined"))
  
  tryCatch({ 
    
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
    Area.peak  = norm.peak.area(conv.trans,conv.max)
    
    # subset spectrum at this offset
    newSpec = spec$flux.pass[conv.max:(conv.max+N-1)]
    
    # rescale left 2/3 of template to spectrum subset
    coefs   = lm.scale(template$flux.pass[1:N_frac],newSpec[1:N_frac])
    newTemp = template$flux.pass*coefs[2]+coefs[1]
    
    # compute raw and rescaled Kolmogorov-Smirnov statistics
    KS.raw = unname(ks.test(newSpec,newTemp)$statistic)
    KS.trans = K_scale(KS.raw)
    
    # combine scores and return
    scores[1:4] = c(Area.peak,KS.raw,KS.trans,Area.peak*KS.trans)
    
    # define y-offset for plotting both together
    y.off = 5*IQR(newSpec)
    
    # get spec name
    spec.name = gsub(".*/|\\.fits?","",sub("^.*(spec[-0-9]*).*$","\\1",spec.file))
    
    # create file for plot output 
    png(file=paste0("plots/",spec.name,".png"),res=300,width=5,height=4,units="in")
    
    # make basic plot of output
    par(cex=.6,las=1)
    plot(1:N,newSpec,type='l',col='black',lwd=1,ylim=c(min(newTemp)-y.off,max(newSpec)),
         ylab="Rescaled Flux (template offset for legibility)",xlab="Wavelength index",main=spec.name)
    lines(1:N,newTemp-y.off*.9,type='l',lwd=1,col='blue')
    legend(x=N,y=min(newTemp)-y.off,legend=c("Spectrum","Template"),
           col=c("black","blue"),lty=1,lwd=2,cex=.8,xjust=.9,yjust=.1)
    
    # close plot device
    dev.off()
    
  }, error = function(e) {
    write(toString(e), stderr())
  })
  
  return(scores)
  
}


# get fits files
files = list.files("data/",full.names=T)

# create/clear plots output directory
if(dir.exists("plots/")){
  unlink("plots/*")
} else dir.create("plots/")

# run analysis, sort by score, and change rownames
res = as.data.frame(t(sapply(files,ScoreSpec)))
res = res[order(res$Combined,decreasing=T),]
rownames(res) = sub("data/","",rownames(res))

# get number prefix of spec
spec.num = sub("(\\d+).*","\\1",grep("\\d{4}\\..*",list.files("."),value=T))[1]

# write output
write.table(capture.output(print(res,digits=3)),
            file=paste0(spec.num,".txt"),quote=F,row.names=F,col.names=F)
