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
low.pass = function(flux,thresh=200,strength=0.03){
  f.raw = Re(fft(fft(flux)*exp(-exp(((1+strength)^2-1)*((1:length(flux))-thresh))),inverse=T))
  ab = lm.scale(f.raw,flux)
  return(f.raw*ab[2]+ab[1])
}

# simple find n-th local-min-or-zero function
nth.local.min.or.zero = function(vec,nth=1){
  # if(0 %in% diff(vec[vec>0])) warning("Found adjacent-duplicate; may cause unexpected behavior!")
  return( which( c(0,diff(sign(diff(vec))),0)==2 | vec==0 )[nth] )
}

# define function to get peak indices
peak.indices = function(vec,index=NULL){
  # define function for finding end of peak
  rising.or.zero.bound = function(vec){
    return(which( c(F,diff(vec)>0) | vec==0 )[1]-1)
  }
  n = length(vec)
  if(is.null(index)) index=which(vec==max(vec))
  return(c(max(1,index-rising.or.zero.bound(rev(vec[1:(index-1)])),na.rm=T),
           min(n,index+rising.or.zero.bound(vec[(index+1):n]),na.rm=T)))
}

# define normalized peak area function
norm.peak.area = function(vec,index){
  bounds = peak.indices(vec,index=index)
  return(sum(vec[bounds[1]:bounds[2]])/sum(vec))
}

# load template and save # rows for use later
load("cb58_template.Rdata")
N = nrow(template)

# anti-gibbs window size; 10 indices seem good enough
G_window = 10

# anti-gibbs cut-off trough count
G_cutoff = 5

# fraction of max-conv-peak used as noise threshold
C_thresh = 1/5

# fraction of template used for scaling (to eliminate quasars with broad high-wavelength emissions)
N_frac   = floor(N*2/3)

# K_thresh is K-stat of spec-1353-..-0579 (noisy cB58); above this, do not penalize
K_thresh = 0.1701054562

# define K-stat rescaling function
K_scale = function(x){1-pmin(pmax(x-K_thresh,0)/(0.8-K_thresh),1)}

# get fits files
files = list.files("data",full.names=T)

ScoreSpec = function(spec.file){
  
  # get spec name
  spec.name = gsub(".*/|\\.fits?","",sub("^.*(spec[-0-9]*).*$","\\1",spec.file))
  
  # cat(paste0("Running ",spec.file,"\n"))
  scores = setNames(rep(NA,3),c("Area.peak","KS.trans","Combined"))
  
  tryCatch({
    
    # read in spectrum and check column names
    spec = readFrameFromFITS(spec.file)
    names(spec) = tolower(names(spec))
    spec = spec[,"flux",drop=F]
    
    # if spec is all zero, skip
    if(all(spec$flux==0)) stop("skip")
    
    # store original M (spectrum length) for later
    orig_M = nrow(spec)
    
    # in case the spectrum isn't sufficiently long for this to work
    tryCatch({
      
      # NEW anti-gibbs code, based on detecting very large amplitude oscillations at edges
      # match expected near beginning of spectrum, not end, so we use a more careful left trim,
      # and a more generous right trim
      gibbs.enhanced = convolve(spec$flux^4,rep(1,G_window),type="filter")
      trim.left  = nth.local.min.or.zero(gibbs.enhanced,G_cutoff)
      trim.right = 150
      
      # apply trims
      spec = spec[trim.left:(nrow(spec)-trim.right+1),,drop=F]
      
    },error=function(e){
      stop("skip")
    })
    
    # smooth spectrum and enhance peaks
    spec$flux.pass = low.pass(spec$flux)
    spec$flux.peak = loess.peak(spec$flux.pass)
    
    # convolve with template
    conv = convolve(spec$flux.peak,template$flux.peak,type="filter")
    
    # deprecated following section; replaced with (hopefully) improved anti-gibbs code above^^^
    # # get start and end indices of max peak of absolute value of conv
    # max.indices = peak.indices(abs(conv))
    # # 
    # # if this max peak is directly at beginning or end of conv vector,
    # # probably due to gibbs phenomenon, so we ignore it
    # if( 1 %in% max.indices | (M-N+1) %in% max.indices ){
    #   conv[max.indices[1]:max.indices[2]] = 0
    # }
    
    # get max of conv without absolute value (check only for positive-correlation match)
    conv.max = which.max(conv)
    
    # transform to emphasize peaks, and find normalized area of peak at max-conv offset
    conv.trans = pmax(0,abs(conv)-conv[conv.max]*C_thresh)^2
    Area.peak  = norm.peak.area(conv.trans,conv.max)
    
    # subset spectrum at this offset
    newSpec = spec$flux.pass[conv.max:(conv.max+N-1)]
    
    # rescale left 2/3 of template to spectrum subset
    coefs   = lm.scale(template$flux.pass[1:N_frac],newSpec[1:N_frac])
    newTemp = template$flux.pass*coefs[2]+coefs[1]
    
    # compute raw and rescaled Kolmogorov-Smirnov statistics
    KS.trans = K_scale( unname(ks.test(newSpec,newTemp)$statistic) )
    
    # combine scores and return
    scores[1:3] = c(Area.peak,KS.trans,Area.peak*KS.trans)
    
    if(scores[3] >= 0.8){
      
      # create file for plot output
      # png(file=sprintf("plots/%.6f_%s.png",scores[3],spec.name),width=5.5,height=4.5,res=300,units="in")
      pdf(file=sprintf("plots/%.6f_%s.pdf",scores[3],spec.name),width=5.5,height=4.5)
      
      # rescale template temporarily to match min/max of spectrum
      newTemp = (newTemp-mean(newTemp))/diff(range(newTemp))*diff(range(newSpec))+mean(newSpec)
      
      # get plot lower limit and offset
      y.min = max(newSpec) - 1.05*(diff(range(newSpec))+diff(range(newTemp)))
      y.off = 1.05*(max(newTemp)-min(newSpec))
      
      # plot output
      xbreaks = 1:N+(conv.max-1)+(trim.left-1)
      par(cex=.6,las=1)
      plot(xbreaks,newSpec,type='l',col='black',lwd=1,ylim=c(y.min,max(newSpec)),
           ylab="Rescaled Flux (template offset for legibility)",
           xlab=sprintf("             Wavelength index   (orig spec length: %d; shifts length: %d; L-trim: %d)",
                        orig_M,length(conv.trans),trim.left),
           main=sprintf("%s  (offset: %d (%.0f%%); A/K/C scores: %.2f, %.2f, %.4f)",
                        spec.name,conv.max,100*(conv.max-1)/(orig_M-N),scores[1],scores[2],scores[3]))
      axis(1,at=seq(100*ceiling(min(xbreaks)/100),100*floor(max(xbreaks)/100),100),labels=NA)
      lines(xbreaks,newTemp-y.off,type='l',lwd=1,col='blue')
      legend(x=max(xbreaks),y=y.min,legend=c("Spectrum","Template"),
             col=c("black","blue"),lty=1,lwd=2,cex=.8,xjust=.9,yjust=.1)
      
      # close plot device
      dev.off()
      
    }
    
  }, error = function(e) {
    if(e$message != "skip"){
      write(paste(spec.name,toString(e)), stderr())
    }
  })
  
  return(scores)
  
}

# create/clear plots output directory
if(dir.exists("plots")){
  unlink("plots/*")
} else dir.create("plots")

# run analysis, sort by score, and change rownames
res = as.data.frame(t(sapply(files,ScoreSpec)))
res = res[order(res$Combined,decreasing=T),]
rownames(res) = sub("data/","",rownames(res))

# get number prefix of spec
spec.num = sub("(\\d+).*","\\1",grep("\\d{4}\\..*",list.files("."),value=T))[1]

# use this weird roundabout way to ensure consistency in printing
options(scipen=999,digits=22)
f = file(paste0(spec.num,".txt"))
writeLines(capture.output(trunc(1e7*res)/1e7),f)
close(f)
