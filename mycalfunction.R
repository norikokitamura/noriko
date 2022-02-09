
myfun2=function (x, xlab, ylab, xlim, ylim,legend = TRUE, subtitles = TRUE,
          cex.subtitles = 0.75, riskdist = FALSE, scat1d.opts = list(nhistSpike = 200),group=TRUE, 
          ...) 
{
  at <- attributes(x)
  if (missing(ylab)) 
    ylab <- if (at$model == "lr") 
      "Actual Probability"
  else paste("Observed", at$yvar.name)
  if (missing(xlab)) {
    if (at$model == "lr") {
      xlab <- paste("Predicted Pr{", at$yvar.name, 
                    sep = "")
      if (at$non.slopes == 1) {
        xlab <- if (at$lev.name == "TRUE") 
          paste(xlab, "}", sep = "")
        else paste(xlab, "=", at$lev.name, "}", 
                   sep = "")
      }
      else xlab <- paste(xlab, ">=", at$lev.name, 
                         "}", sep = "")
    }
    else xlab <- paste("Predicted", at$yvar.name)
  }
  p <- x[, "predy"]
  p.app <- x[, "calibrated.orig"]
  p.cal <- x[, "calibrated.corrected"]
  if (missing(xlim) & missing(ylim)) 
    xlim <- ylim <- range(c(p, p.app, p.cal), na.rm = TRUE)
    #xlim <- ylim <- range(c(0, 1, 0.1), na.rm = TRUE)
  else {
    if (missing(xlim)) 
      xlim <- range(p)
    if (missing(ylim)) 
       ylim <- range(c(p.app, p.cal, na.rm = TRUE))
  }
  plot(p, p.app, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
       type = "n", ...)
  predicted <- at$predicted
  err <- NULL
  if (length(predicted)) {
    s <- !is.na(p + p.cal)
    err <- predicted - approx(p[s], p.cal[s], xout = predicted, 
                              ties = mean)$y
    cat("\nn=", n <- length(err), "   Mean absolute error=", 
        round(mae <- mean(abs(err), na.rm = TRUE), 3), "   Mean squared error=", 
        round(mean(err^2, na.rm = TRUE), 5), "\n0.9 Quantile of absolute error=", 
        round(quantile(abs(err), 0.9, na.rm = TRUE), 3), 
        "\n\n", sep = "")
    if (subtitles) 
      title(sub = paste("Mean absolute error=", round(mae, 
                                                      3), " n=", n, sep = ""), cex.sub = cex.subtitles, 
            adj = 1)
    if (riskdist) 
      do.call("scat1d", c(list(x = predicted), scat1d.opts))
  }
  lines(p, p.app, lty = 3,lwd=1) #caroline used lwd=3
  lines(p, p.cal, lty = 1,lwd=1)
  abline(a = 0, b = 1, lty = 2,lwd=1)
  if (subtitles) 
    title(sub = paste("B=", at$B, "repetitions,", 
                      at$method), cex.sub = cex.subtitles, adj = 0)
  if (!(is.logical(legend) && !legend)) {
    if (is.logical(legend)) 
      legend <- list(x = 0.80, y = 0.8)
    legend(legend, c("Apparent", "Bias-corrected", 
                     "Ideal"), lty = c(3, 1, 2), lwd=c(1,1,1),bty = "n",cex=1.2)#0.7 to 1.2
  }
  invisible(err)
}
