library(rdd)
library(dplyr)

# custom DCdensity function 

DCdensity_custom<-function (runvar, cutpoint, bin = NULL, bw = NULL, verbose = FALSE, 
                            plot = TRUE, ext.out = FALSE, htest = FALSE, ylim = NULL, xlim = NULL, xlab = NULL, ylab = NULL, main = NULL, col = "deepskyblue4") 
{
  runvar <- runvar[complete.cases(runvar)]
  rn <- length(runvar)
  rsd <- sd(runvar)
  rmin <- min(runvar)
  rmax <- max(runvar)
  if (missing(cutpoint)) {
    if (verbose) 
      cat("Assuming cutpoint of zero.\n")
    cutpoint <- 0
  }
  if (cutpoint <= rmin | cutpoint >= rmax) {
    stop("Cutpoint must lie within range of runvar")
  }
  if (is.null(bin)) {
    bin <- 2 * rsd * rn^(-1/2)
    if (verbose) 
      cat("Using calculated bin size: ", sprintf("%.3f", 
                                                 bin), "\n")
  }
  l <- floor((rmin - cutpoint)/bin) * bin + bin/2 + cutpoint
  r <- floor((rmax - cutpoint)/bin) * bin + bin/2 + cutpoint
  lc <- cutpoint - (bin/2)
  rc <- cutpoint + (bin/2)
  j <- floor((rmax - rmin)/bin) + 2
  binnum <- round((((floor((runvar - cutpoint)/bin) * bin + 
                       bin/2 + cutpoint) - l)/bin) + 1)
  cellval <- rep(0, j)
  for (i in seq(1, rn)) {
    cnum <- binnum[i]
    cellval[cnum] <- cellval[cnum] + 1
  }
  cellval <- (cellval/rn)/bin
  cellmp <- seq(from = 1, to = j, by = 1)
  cellmp <- floor(((l + (cellmp - 1) * bin) - cutpoint)/bin) * 
    bin + bin/2 + cutpoint
  if (is.null(bw)) {
    leftofc <- round((((floor((lc - cutpoint)/bin) * bin + 
                          bin/2 + cutpoint) - l)/bin) + 1)
    rightofc <- round((((floor((rc - cutpoint)/bin) * bin + 
                           bin/2 + cutpoint) - l)/bin) + 1)
    if (rightofc - leftofc != 1) {
      stop("Error occurred in bandwidth calculation")
    }
    cellmpleft <- cellmp[1:leftofc]
    cellmpright <- cellmp[rightofc:j]
    P.lm <- lm(cellval ~ poly(cellmp, degree = 4, raw = T), 
               subset = cellmp < cutpoint)
    mse4 <- summary(P.lm)$sigma^2
    lcoef <- coef(P.lm)
    fppleft <- 2 * lcoef[3] + 6 * lcoef[4] * cellmpleft + 
      12 * lcoef[5] * cellmpleft * cellmpleft
    hleft <- 3.348 * (mse4 * (cutpoint - l)/sum(fppleft * 
                                                  fppleft))^(1/5)
    P.lm <- lm(cellval ~ poly(cellmp, degree = 4, raw = T), 
               subset = cellmp >= cutpoint)
    mse4 <- summary(P.lm)$sigma^2
    rcoef <- coef(P.lm)
    fppright <- 2 * rcoef[3] + 6 * rcoef[4] * cellmpright + 
      12 * rcoef[5] * cellmpright * cellmpright
    hright <- 3.348 * (mse4 * (r - cutpoint)/sum(fppright * 
                                                   fppright))^(1/5)
    bw = 0.5 * (hleft + hright)
    if (verbose) 
      cat("Using calculated bandwidth: ", sprintf("%.3f", 
                                                  bw), "\n")
  }
  if (sum(runvar > cutpoint - bw & runvar < cutpoint) == 0 | 
      sum(runvar < cutpoint + bw & runvar >= cutpoint) == 0) 
    stop("Insufficient data within the bandwidth.")
  if (plot) {
    d.l <- data.frame(cellmp = cellmp[cellmp < cutpoint], 
                      cellval = cellval[cellmp < cutpoint], dist = NA, 
                      est = NA, lwr = NA, upr = NA)
    pmin <- cutpoint - 2 * rsd
    pmax <- cutpoint + 2 * rsd
    for (i in 1:nrow(d.l)) {
      d.l$dist <- d.l$cellmp - d.l[i, "cellmp"]
      w <- kernelwts(d.l$dist, 0, bw, kernel = "triangular")
      newd <- data.frame(dist = 0)
      pred <- predict(lm(cellval ~ dist, weights = w, data = d.l), 
                      interval = "confidence", newdata = newd)
      d.l$est[i] <- pred[1]
      d.l$lwr[i] <- pred[2]
      d.l$upr[i] <- pred[3]
    }
    d.r <- data.frame(cellmp = cellmp[cellmp >= cutpoint], 
                      cellval = cellval[cellmp >= cutpoint], dist = NA, 
                      est = NA, lwr = NA, upr = NA)
    for (i in 1:nrow(d.r)) {
      d.r$dist <- d.r$cellmp - d.r[i, "cellmp"]
      w <- kernelwts(d.r$dist, 0, bw, kernel = "triangular")
      newd <- data.frame(dist = 0)
      pred <- predict(lm(cellval ~ dist, weights = w, data = d.r), 
                      interval = "confidence", newdata = newd)
      d.r$est[i] <- pred[1]
      d.r$lwr[i] <- pred[2]
      d.r$upr[i] <- pred[3]
    }
    if (is.null(ylim)){
      plot(d.l$cellmp, d.l$est, lty = 1, lwd = 2, col = col, 
           type = "l", xlim = c(pmin, pmax), ylim = c(min(cellval[cellmp <= 
                                                                    pmax & cellmp >= pmin]), max(cellval[cellmp <= 
                                                                                                           pmax & cellmp >= pmin])), xlab = NA, ylab = NA, 
           main = NA)
      lines(d.l$cellmp, d.l$lwr, lty = 2, lwd = 1, col = col, 
            type = "l")
      lines(d.l$cellmp, d.l$upr, lty = 2, lwd = 1, col = col, 
            type = "l")
      lines(d.r$cellmp, d.r$est, lty = 1, lwd = 2, col = col, 
            type = "l")
      lines(d.r$cellmp, d.r$lwr, lty = 2, lwd = 1, col = col, 
            type = "l")
      lines(d.r$cellmp, d.r$upr, lty = 2, lwd = 1, col = col, 
            type = "l")
      points(cellmp, cellval, type = "p", pch = 20)
    } else {
      plot(d.l$cellmp, d.l$est, lty = 1, lwd = 2, col = col, 
           type = "l", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
           main = main, pch=20)
      lines(d.l$cellmp, d.l$lwr, lty = 2, lwd = 1, col = col, 
            type = "l")
      lines(d.l$cellmp, d.l$upr, lty = 2, lwd = 1, col = col, 
            type = "l")
      lines(d.r$cellmp, d.r$est, lty = 1, lwd = 2, col = col, 
            type = "l")
      lines(d.r$cellmp, d.r$lwr, lty = 2, lwd = 1, col = col, 
            type = "l")
      lines(d.r$cellmp, d.r$upr, lty = 2, lwd = 1, col = col, 
            type = "l")
      abline(v=0, col="red", lty=2)
      points(cellmp, cellval, type = "p", pch = 20, col = col)
    }
  }
  cmp <- cellmp
  cval <- cellval
  padzeros <- ceiling(bw/bin)
  jp <- j + 2 * padzeros
  if (padzeros >= 1) {
    cval <- c(rep(0, padzeros), cellval, rep(0, padzeros))
    cmp <- c(seq(l - padzeros * bin, l - bin, bin), cellmp, 
             seq(r + bin, r + padzeros * bin, bin))
  }
  dist <- cmp - cutpoint
  w <- 1 - abs(dist/bw)
  w <- ifelse(w > 0, w * (cmp < cutpoint), 0)
  w <- (w/sum(w)) * jp
  fhatl <- predict(lm(cval ~ dist, weights = w), newdata = data.frame(dist = 0))[[1]]
  w <- 1 - abs(dist/bw)
  w <- ifelse(w > 0, w * (cmp >= cutpoint), 0)
  w <- (w/sum(w)) * jp
  fhatr <- predict(lm(cval ~ dist, weights = w), newdata = data.frame(dist = 0))[[1]]
  thetahat <- log(fhatr) - log(fhatl)
  sethetahat <- sqrt((1/(rn * bw)) * (24/5) * ((1/fhatr) + 
                                                 (1/fhatl)))
  z <- thetahat/sethetahat
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  if (verbose) {
    cat("Log difference in heights is ", sprintf("%.3f", 
                                                 thetahat), " with SE ", sprintf("%.3f", sethetahat), 
        "\n")
    cat("  this gives a z-stat of ", sprintf("%.3f", z), 
        "\n")
    cat("  and a p value of ", sprintf("%.3f", p), "\n")
  }
  if (ext.out) 
    return(list(theta = thetahat, se = sethetahat, z = z, 
                p = p, binsize = bin, bw = bw, cutpoint = cutpoint, 
                data = data.frame(cellmp, cellval)))
  else if (htest) {
    structure(list(statistic = c(z = z), p.value = p, method = "McCrary (2008) sorting test", 
                   parameter = c(binwidth = bin, bandwidth = bw, cutpoint = cutpoint), 
                   alternative = "no apparent sorting"), class = "htest")
  }
  else return(p)
}


# custom plot.RD function

plot.RD <- function(x,gran=400,bins=100,which=1, ylab = NULL, xlab = NULL, main = NULL, col = "deepskyblue4", range,...) {
  frm<-FALSE
  if("frame" %in% names(x$call)) frm<-eval.parent(x$call$frame)
  if(!frm){
    x$call$frame<-TRUE
    x$call$verbose<-FALSE
    x<-eval.parent(x$call)
  }
  d<-as.data.frame(x$frame)
  
  if(length(x$na.action)>0)
    d<-d[-x$na.action,]
  
  if("kernel" %in% names(x$call)) 
    kern<-eval.parent(x$call$kernel)
  else 
    kern<-"triangular"
  
  if("cutpoint" %in% names(x$call)) 
    cut<-eval.parent(x$call$cutpoint)
  else
    cut<-0
  
  bw<-x$bw[1]
  
  if(missing(range)) {
    range<-c(cut-10*bw,cut+10*bw)
    if(range[1]<min(d$X)) range[1]<-min(d$X)
    if(range[2]>max(d$X)) range[2]<-max(d$X)
  }
  
  if(range[1]=="min")
    range[1]<-min(d$X)
  if(range[2]=="max")
    range[2]<-max(d$X)
  range<-as.double(range)
  
  rdplot<-function(d) {
    d.l<-data.frame(X=d$X[d$X<cut],Y=d$Y[d$X<cut])
    lval<-seq(range[1],cut,length.out=(gran%/%2))
    lest<-vector(length=(gran%/%2))
    llwr<-vector(length=(gran%/%2))
    lupr<-vector(length=(gran%/%2))
    for(i in 1:(gran%/%2)) {
      sub<-d.l$X>=(lval[i]-bw) & d.l$X<=(lval[i]+bw)
      w<-kernelwts(X=d.l$X[sub],center=lval[i],bw=bw,kernel=kern)
      ly<-d.l$Y[sub]
      lx<-d.l$X[sub]
      if(length(lx)<=2)
        pred<-rep(NA,3)
      else
        pred<-predict(lm(ly~lx,weights=w),interval="confidence",newdata=data.frame(lx=lval[i]))
      lest[i]<-pred[1]
      llwr[i]<-pred[2]
      lupr[i]<-pred[3]
    }
    
    d.r<-data.frame(X=d$X[d$X>=cut],Y=d$Y[d$X>=cut])
    rval<-seq(cut,range[2],length.out=(gran%/%2))
    rest<-vector(length=(gran%/%2))
    rlwr<-vector(length=(gran%/%2))
    rupr<-vector(length=(gran%/%2))
    for(i in 1:(gran%/%2)) {
      sub<-d.r$X>=(rval[i]-bw) & d.r$X<=(rval[i]+bw)
      w<-kernelwts(X=d.r$X[sub],center=rval[i],bw=bw,kernel=kern)
      ry<-d.r$Y[sub]
      rx<-d.r$X[sub]
      if(length(rx)<=2)
        pred<-rep(NA,3)
      else
        pred<-predict(lm(ry~rx,weights=w),interval="confidence",newdata=data.frame(rx=rval[i]))
      rest[i]<-pred[1]
      rlwr[i]<-pred[2]
      rupr[i]<-pred[3]
    }
    
    #plot to the left
    if(length(unique(d$Y))==2) {
      #DO THIS for when the outcome is dichotomous
      ep<-(max(d$X)-min(d$X))/(2*bins)
      nX<-seq(min(d$X)-ep,max(d$X)+ep,length=bins+1)
      nY<-rep(NA,length(nX))
      for(i in (1:(length(nX)-1))){
        if(sum(!is.na(d$Y[d$X>nX[i] & d$X<=nX[i+1]]))==0)
          next
        nY[i]<-sum(d$Y[d$X>nX[i] & d$X<=nX[i+1]],na.rm=TRUE)/sum(!is.na(d$Y[d$X>nX[i] & d$X<=nX[i+1]]))
      }
      sub<-nX>=range[1] & nX<=range[2]
      subl<-lval>=range[1] & lval<=range[2]
      subr<-rval>=range[1] & rval<=range[2]
      plot(nX,nY,
           type="p",pch=20,cex=.5,col=col,
           xlim=c(range[1],range[2]),
           ylim=c(min(c(llwr[subl],rlwr[subr]),na.rm=T),
                  max(c(lupr[subl],rupr[subr]),na.rm=T)),
           xlab=xlab,
           ylab=ylab,
           main=main
      )
    } else {
      subl<-lval>=range[1] & lval<=range[2]
      subr<-rval>=range[1] & rval<=range[2]
      plot(d$X,d$Y,
           type="p",pch=20,cex=.5,col=col,
           xlim=c(range[1],range[2]),
           ylim=c(min(c(llwr[subl],rlwr[subr]),na.rm=T),
                  max(c(lupr[subl],rupr[subr]),na.rm=T)),
           xlab=xlab,
           ylab=ylab,
           main=main
      )
    } 
    #plot to the left
    lines(lval,lest,
          lty=1,lwd=2,col=col,type="l"
    )
    
    lines(lval,llwr,
          lty=2,lwd=1,col=col,type="l"
    )
    lines(lval,lupr,
          lty=2,lwd=1,col=col,type="l"
    )
    
    #plot to the right
    lines(rval,rest,
          lty=1,lwd=2,col=col,type="l"
    )
    lines(rval,rlwr,
          lty=2,lwd=1,col=col,type="l"
    )
    lines(rval,rupr,
          lty=2,lwd=1,col=col,type="l"
    )
  }
  if(x$type=="sharp" | 1%in%which){
    rdplot(d)
    dev.flush()
  }
  if(x$type=="fuzzy" & 2%in%which){
    d$Y<-d$Z
    rdplot(d)
    dev.flush()
  }
}

#Read data from csv
data <- read.csv("All_States_GE.csv")

##Prepare data for McCrary Test 

#Take data for LS2019 Elections and combine Constituency Name with State Name since Constituency Names may be the same
data_f <- data[data$Assembly_No == 17 & data$Poll_No == 0, ]
data_f$Constituency_Name <- paste(data_f$Constituency_Name, "_", data_f$State_Name, sep = "")

#Calculate the vote share of the winner and the runner up of each seat
data_w <- data_f[data_f$Position == 1, c('Constituency_Name', 'Vote_Share_Percentage')]
colnames(data_w) <- c('Constituency_Name', 'Vote_Share_Percentage_W')
data_r <- data_f[data_f$Position == 2, c('Constituency_Name', 'Vote_Share_Percentage')]
colnames(data_r) <- c('Constituency_Name', 'Vote_Share_Percentage_RU')
data_f$Vote_Share_Percentage_W <- data_w$Vote_Share_Percentage_W[match(data_f$Constituency_Name, data_w$Constituency_Name)]
data_f$Vote_Share_Percentage_RU <- data_r$Vote_Share_Percentage_RU[match(data_f$Constituency_Name, data_w$Constituency_Name)]

#Calculate Win Margin for each candidate
data_f$Win_Margin <- ifelse(data_f$Position == 1, data_f$Vote_Share_Percentage - data_f$Vote_Share_Percentage_RU, data_f$Vote_Share_Percentage - data_f$Vote_Share_Percentage_W)
data_f$Win_Margin <- data_f$Win_Margin/100

#Filter for BJP
data_f <- data_f[data_f$Party == "BJP", ]

BJP_ruled_UT <- c("Assam", "Bihar", "Goa", "Gujarat", "Haryana", "Himachal_Pradesh", "Jharkhand", "Maharashtra",
                  "Manipur", "Nagaland", "Tripura", "Uttar_Pradesh", "Uttarakhand", "Dadra_&_Nagar_Haveli", "Daman_&_Diu", "Andaman_&_Nicobar_Islands", "Lakshadweep", "Chandigarh")

BJP_ruled <- c("Assam", "Bihar", "Goa", "Gujarat", "Haryana", "Himachal_Pradesh", "Jharkhand", "Maharashtra",
               "Manipur", "Nagaland", "Tripura", "Uttar_Pradesh", "Uttarakhand")


data_fb <- data_f[data_f$State_Name %in% BJP_ruled, ]
data_fnb <- data_f[!(data_f$State_Name %in% BJP_ruled), ]

#Run McCrary test. (Bandwidth (bw) can be changed here)
DCdensity_custom(data_f$Win_Margin, ylim=c(0,4), xlim=c(-0.5, 0.5), xlab = "BJP Win Margin", ylab = "Density", main = "General Election 2019")
DCdensity_custom(data_fb$Win_Margin, ylim=c(0,4), xlim=c(-0.5, 0.5), xlab = "BJP Win Margin", ylab = "Density", main = "General Election 2019: BJP Ruled States")
DCdensity_custom(data_fnb$Win_Margin, ylim=c(0,4), xlim=c(-0.5, 0.5), xlab = "BJP Win Margin", ylab = "Density", main = "General Election 2019: Non-BJP Ruled States")

##Prepare data for RDD 
data_f17e <- data %>%
  filter(Assembly_No == 17, Poll_No == 0, Position == 1) %>%
  mutate(Constituency_Name = paste(Constituency_Name, State_Name, sep = "_")) %>%
  select(Constituency_Name, Electors) %>%
  rename(Electors_19 = Electors)

data_f16e <- data %>%
  filter(Assembly_No == 16, Poll_No == 0, Position == 1) %>%
  mutate(Constituency_Name = paste(Constituency_Name, State_Name, sep = "_")) %>%
  select(Constituency_Name, Electors) %>%
  rename(Electors_14 = Electors)

data_g <- data_f %>%
  select(Constituency_Name, Win_Margin)


# Fix discrepancies in Constituency Names between elections
repl_col <- c("BARDHAMAN DURGAPUR_West_Bengal" = "BURDWAN - DURGAPUR_West_Bengal",
              "BIKANER (SC)_Rajasthan" = "BIKANER_Rajasthan",
              "CHEVELLA_Andhra_Pradesh" = "CHELVELLA_Andhra_Pradesh",
              "DADRA AND NAGAR HAVELI_Dadra_&_Nagar_Haveli" = "DADAR & NAGAR HAVELI_Dadra_&_Nagar_Haveli",
              "JAYNAGAR_West_Bengal" = "JOYNAGAR_West_Bengal")

data_f17e$Constituency_Name <- gsub("^(.*?)_Telangana", "\\1_Andhra_Pradesh", data_f17e$Constituency_Name)
data_f17e$Constituency_Name <- ifelse(data_f17e$Constituency_Name %in% names(repl_col), repl_col[data_f17e$Constituency_Name], data_f17e$Constituency_Name)

data_g$Constituency_Name <- gsub("^(.*?)_Telangana", "\\1_Andhra_Pradesh", data_g$Constituency_Name)
data_g$Constituency_Name <- ifelse(data_g$Constituency_Name %in% names(repl_col), repl_col[data_g$Constituency_Name], data_g$Constituency_Name)

data_e <- inner_join(data_f16e, data_f17e, by = "Constituency_Name") %>%
  mutate(g = (Electors_19 - Electors_14) / Electors_14)

data_e <- inner_join(data_e, data_g, by = 'Constituency_Name')

x<-data_e$Win_Margin
y<-data_e$g*100
res<-RDestimate(y~x, bw=0.107)
print(res$est)
print(res$se)
print(res$p)
plot(res, xlab = "BJP Win Margin", ylab = "Electorate Growth (%)", col = "deepskyblue4", main = "Mean Electorate Growth 2014-19 = 9.2% " )
abline(v=0, col="red", lty=2)

data_ex = data_e %>% filter(!(Win_Margin < 0) | !(Win_Margin > -0.015) | !(g > 0.15))

x<-data_ex$Win_Margin
y<-data_ex$g*100
res<-RDestimate(y~x, bw=0.107)
print(res$est)
res$se
res$p
plot(res, xlab = "BJP Win Margin", ylab = "Electorate Growth (%)", col = "deepskyblue4", main = "Mean Electorate Growth 2014-19 = 9.2% " )
abline(v=0, col="red", lty=2)