###################################################################################
##' Function for visualizing the clustering results
##'
##' @param res object return by function clusvis
##' @param add.centers boolean. If TRUE, centers are plotted.
##' @param add.obs   boolean. If TRUE, coordinnates of the observations are plotted.
##' @param threshold   numeric. It contains the thersholds used for computing the level curves.
##' @param dim   numeric. This vector of size two choose the axes to represent.
##' @param col.dots   numeric/character. It specifies the color of the dots.
##' @param col.centers  numeric/character. It specifies the color of the centers.
##' @param cex.dots   numeric. It specifies the size of the observations.
##' @param ...   other parameters of the plot function.
##'
##'
##' @return NULL
##' @examples
##' set.seed(123)
##' @export
##'
plotDensityClusVisu <- function(res,
                                dim=c(1,2),
                                threshold=0.95,
                                add.obs=FALSE,
                                positionlegend="topright",
                                xlim=NULL,
                                ylim=NULL){
  if (add.obs){
    input <- list()
    if (is.null(xlim)) xlim <- c(min(c(res$centers[,dim[1]], res$y[,dim[1]]))-1, max(c(res$centers[,dim[1]], res$y[,dim[1]])) + 1)
    if (is.null(ylim)) ylim <- c(min(c(res$centers[,dim[2]], res$y[,dim[2]]))-1, max(c(res$centers[,dim[2]], res$y[,dim[2]])) + 1)
    input$xval <- seq(xlim[1], xlim[2], length.out = 400)
    input$yval <- seq(ylim[1], ylim[2], length.out = 400)
    tmp <- sapply(1:length(res$prop), function(k) as.numeric(outer(input$xval, input$yval, dmixtmvnorm, mu=res$centers[k,dim], prop=res$prop[k])))
    input$z  <- t(matrix(rowSums(tmp), length(input$xval), length(input$yval)))
    input$z <- (matrix(apply(sweep(tmp, 1, rowSums(tmp), "/"), 1, max), length(input$xval), length(input$yval)))
    input$z[which(input$z<0.55)] <- 0.4
    contour(input,
            main = paste0("Difference between entropies: ", round(res$EM - res$EV,2)),
            levels=c(0.95,0.8,0.5),
            col=c("gray30","gray30",1),
            lwd=c(1,1,2),
            lty=c(2,2,1),
            labcex = 0.8,
            xlab=paste0("Dim1 (", 100*round(res$inertia[dim[1]]/sum(res$inertia), 4), "%)"),
            ylab=paste0("Dim2 (", 100*round(res$inertia[dim[2]]/sum(res$inertia), 4), "%)"))
    colset <- c("darkorange1", "dodgerblue2", "black", "chartreuse2", "darkorchid2", "gold2", "deeppink2", "deepskyblue1", "firebrick2", "cyan1")
    if (length(colset) < ncol(res$logtik.obs)) colset <- rep("black", ncol(res$logtik.obs))
    # text(res$y[,dim[1]],
    #      res$y[,dim[2]],
    #      as.character(apply(res$logtik.obs, 1, which.max)),
    #      cex=0.6,
    #      lwd=1.2,
    #      col=colset[apply(res$logtik.obs, 1, which.max)])
    points(res$y[,dim[1]],
         res$y[,dim[2]],
         pch=20,
         cex=0.7,
         col=colset[apply(res$logtik.obs, 1, which.max)])
    if (length(colset) >= ncol(res$logtik.obs)) legend(x = positionlegend, legend = paste0("Compo.", 1:ncol(res$logtik.obs)), col = colset[1:ncol(res$logtik.obs)], pch = 20)


  }else{
    if (is.null(xlim)) xlim <- c(min(res$centers[,dim[1]])-4, max(res$centers[,dim[1]]) + 4)
    if (is.null(ylim)) ylim <- c(min(res$centers[,dim[2]])-4, max(res$centers[,dim[2]]) + 4)
    xval <- seq(xlim[1], xlim[2], length.out = 400)
    yval <- seq(ylim[1], ylim[2], length.out = 400)
    tmp <- sapply(1:length(res$prop), function(k) as.numeric(outer(xval, yval, dmixtmvnorm, mu=res$centers[k,dim], prop=res$prop[k])))
    z <- t(matrix(rowSums(tmp), length(xval), length(yval)))
    tikmax <- t(matrix(apply(sweep(tmp, 1, rowSums(tmp), "/"), 1, max), length(xval), length(yval)))
    class <- t(matrix(apply(sweep(tmp, 1, rowSums(tmp), "/"), 1, which.max), length(xval), length(yval)))
    tmp <- optimize(function(alpha, z, threshold) abs(sum(z*(z>alpha))/sum(z) - threshold),
                    interval = c(0, max(z)),
                    z=z,
                    threshold=threshold)
    bound <- min(as.numeric(tikmax)[which(as.numeric(z)>tmp$minimum)])

    tikmax <- tikmax *(z>tmp$minimum)
    image(xval,yval,t(tikmax),
          main =  paste0("Difference between entropies: ", round(res$EM - res$EV,2)),
          col=c("white","gray30","gray60","gray80"), breaks = c(0,0.001,0.8,0.95,1),
          xlab=paste0("Dim1 (", 100*round(res$inertia[dim[1]]/sum(res$inertia), 4), "%)"),
          ylab=paste0("Dim2 (", 100*round(res$inertia[dim[2]]/sum(res$inertia), 4), "%)")
    )
    legend(legend = c("0.95<Pr. Classif.", "0.8<Pr. Classif.<0.95", "Pr. Classif.<0.8", "outside the conf. level."), x = positionlegend, fill = c("gray80","gray60","gray30", "white"), cex=0.7)
    input <- list(x=xval, y=yval, z=t(tikmax))

    contour(input, add=TRUE, levels=c(0.95,0.8,0.001), drawlabels = F,
            lwd=c(1,1,2),
            lty=c(2,2,1),
            labcex = 0.8)
    text(res$centers[,dim[1]], res$centers[,dim[2]], as.character(1:nrow(res$centers)), lwd=1.2)
  }
}
  #
  # plotDensityClusVisu <- function(res,
  #                                 dim=c(1,2),
  #                                 threshold=0.95,
  #                                 add.obs=FALSE,
  #                                 add.iso=1-add.obs){
  #   xval <- seq(min(c(res$centers[,dim[1]], res$y[,dim[1]]))-4, max(c(res$centers[,dim[1]], res$y[,dim[1]])) + 4, length.out = 200)
  #   yval <- seq(min(c(res$centers[,dim[2]], res$y[,dim[2]]))-4, max(c(res$centers[,dim[2]], res$y[,dim[2]])) + 4, length.out = 200)
  #
  #   tmp <- sapply(1:length(res$prop), function(k) as.numeric(outer(xval, yval, dmixtmvnorm, mu=res$centers[k,dim], prop=res$prop[k])))
  #   z <- t(matrix(rowSums(tmp), length(xval), length(yval)))
  #   tikmax <- t(matrix(apply(sweep(tmp, 1, rowSums(tmp), "/"), 1, max), length(xval), length(yval)))
  #   class <- t(matrix(apply(sweep(tmp, 1, rowSums(tmp), "/"), 1, which.max), length(xval), length(yval)))
  #   tmp <- optimize(function(alpha, z, threshold) abs(sum(z*(z>alpha))/sum(z) - threshold),
  #                   interval = c(0, max(z)),
  #                   z=z,
  #                   threshold=threshold)
  #   bound <- min(as.numeric(tikmax)[which(as.numeric(z)>tmp$minimum)])
  #
  #   areatik <- (tikmax<threshold) * 0.5
  #
  #   f <- list(
  #     family = "sans serif",
  #     size = 14,
  #     color = "black"
  #   )
  #
  #   xaxis <- list(
  #     title = paste0("Dim1 (", 100*round(res$inertia[dim[1]]/sum(res$inertia), 4), "%)"),
  #     titlefont = f
  #   )
  #   yaxis <- list(
  #     title = paste0("Dim2 (", 100*round(res$inertia[dim[2]]/sum(res$inertia), 4), "%)"),
  #     titlefont = f
  #   )
  #   p <- plot_ly(data=as.data.frame(res$y), showlegend=FALSE) %>%
  #     add_contour(z = areatik,
  #                 x=xval,
  #                 y=yval,
  #                 colors = colorRamp(c("white", "grey60")),
  #                 hoverinfo = 'text',
  #                 ncontours = 20,
  #                 showscale=FALSE ,
  #                 line=list(color="black"),
  #                 text = matrix(paste0("<br>Classif:",
  #                                      class,
  #                                      "<br>Prob. error:",
  #                                      round(1-tikmax, 3)), nrow(z), ncol(z))) %>%
  #     layout( xaxis= xaxis, yaxis=yaxis)
  #   p
  #   if (add.iso)
  #     p <-  add_contour(p,
  #                       z = 10*(z>tmp$minimum),
  #                       x=xval,
  #                       y=yval,
  #                       hoverinfo = 'none',
  #                       colors=colorRamp(c("black", "black")),
  #                       ncontours = 20,
  #                       showscale=FALSE ,
  #                       contours=list(coloring="heatmap"),
  #                       line=list(color="black", dash="dot"))%>%
  #     add_trace(p,
  #               x=res$centers[,dim[1]],
  #               y=res$centers[,dim[2]],
  #               type="scatter",
  #               mode="text",
  #               hoverinfo = 'none',
  #               textfont = list(color = 'red',
  #                               family = 'sans serif',
  #                               size = 16,
  #                               font = "white"),
  #               text=as.character(1:nrow(res$centers)))%>%
  #     add_trace(p,
  #               x=res$centers[,dim[1]],
  #               y=res$centers[,dim[2]],
  #               type="scatter",
  #               opacity=0,
  #               mode="markers",
  #               hoverinfo = 'text',
  #               text = paste0("Center of Class ",
  #                             1:nrow(res$centers),
  #                             "<br>Dim", dim[1], ": ", round(res$centers[,dim[1]], 2),
  #                             "<br>Dim", dim[2], ": ", round(res$centers[,dim[2]], 2))
  #     )
  #
  #
  #   if (add.obs)
  #     p <- add_trace(p,
  #                    x=res$y[,dim[1]],
  #                    y=res$y[,dim[2]],
  #                    type="scatter",
  #                    mode="text",
  #                    hoverinfo = 'none',
  #                    textfont = list(color = apply(res$logtik.obs, 1, which.max),
  #                                    family = 'sans serif',
  #                                    size = 10,
  #                                    font = "white"),
  #                    text=as.character(apply(res$logtik.obs, 1, which.max))) %>%
  #     add_trace(p,
  #               x=res$y[,dim[1]],
  #               y=res$y[,dim[2]],
  #               type="scatter",
  #               marker=list(color="blue"),
  #               opacity=0,
  #               mode="markers",
  #               hoverinfo = 'text',
  #               text = paste0("Obs:     ",
  #                             1:nrow(res$y),
  #                             "<br>Classif.:",
  #                             apply(res$logtik.obs, 1, which.max),
  #                             "<br>Prob. error:",
  #                             round(1-apply(exp(res$logtik.obs), 1, max), 3)))
  #   p
  # }
  #


  # plotDensityClusVisu <- function(res,
  #                                 dim=c(1,2),
  #                                 threshold=0.99,
  #                                 add.obs=FALSE,
  #                                 add.modes=TRUE){
  # xval <- seq(min(c(res$centers[,dim[1]], res$y[,dim[1]]))-2, max(c(res$centers[,dim[1]], res$y[,dim[1]])) + 2, length.out = 200)
  # yval <- seq(min(c(res$centers[,dim[2]], res$y[,dim[2]]))-2, max(c(res$centers[,dim[2]], res$y[,dim[2]])) + 2, length.out = 200)
  #
  # tmp <- sapply(1:length(res$prop), function(k) as.numeric(outer(xval, yval, dmixtmvnorm, mu=res$centers[k,dim], prop=res$prop[k])))
  # z <- t(matrix(rowSums(tmp), length(xval), length(yval)))
  # tikmax <- t(matrix(apply(sweep(tmp, 1, rowSums(tmp), "/"), 1, max), length(xval), length(yval)))
  # class <- t(matrix(apply(sweep(tmp, 1, rowSums(tmp), "/"), 1, which.max), length(xval), length(yval)))
  # tmp <- optimize(function(alpha, z, threshold) abs(sum(z*(z>alpha))/sum(z) - threshold),
  #                 interval = c(0, max(z)),
  #                 z=z,
  #                 threshold=threshold)
  # bound <- min(as.numeric(tikmax)[which(as.numeric(z)>tmp$minimum)])
  #
  # area <- (z>tmp$minimum) * (tikmax - bound) + (z>tmp$minimum) * 0.1
  #
  # f <- list(
  #   family = "sans serif",
  #   size = 14,
  #   color = "black"
  # )
  #
  # a <- list(
  #   x = res$centers[,dim[1]],
  #   y = res$centers[,dim[2]],
  #   text = 1:nrow(res$centers),
  #   xref = "x",
  #   yref = "y",
  #   showarrow = F,
  #   ax = 0,
  #   ay = 0,
  #   font = list(color = 'red',
  #               family = 'sans serif',
  #               size = 16,
  #               font = "white")
  # )
  #
  # xaxis <- list(
  #   title = paste0("Dim1 (", 100*round(res$inertia[dim[1]]/sum(res$inertia), 4), "%)"),
  #   titlefont = f
  # )
  # yaxis <- list(
  #   title = paste0("Dim2 (", 100*round(res$inertia[dim[2]]/sum(res$inertia), 4), "%)"),
  #   titlefont = f
  # )
  # p <- plot_ly(data=as.data.frame(res$y)) %>%
  #   add_contour(z = area,
  #               x=xval,
  #               y=yval,
  #               type = "contour",
  #               colors = colorRamp(c("white", "grey30")),
  #               hoverinfo = 'text',
  #               ncontours = 20,
  #               showscale=FALSE ,
  #               text = matrix(paste0("log-density:     ",
  #                                    round(log(z),3),
  #                                    "<br>Classif:",
  #                                    class,
  #                                    "<br>Proportion:",
  #                                    round(res$prop[as.numeric(class)],2),
  #                                    "<br>Prob. error:",
  #                                    round(1-tikmax, 3)), nrow(z), ncol(z))) %>%
  #   layout( annotations = a, xaxis= xaxis, yaxis=yaxis)
  #   if (add.obs)
  #     p <- add_trace(p,
  #                    x=res$y[,dim[1]],
  #                    y=res$y[,dim[2]],
  #                    type="scatter",
  #                    marker=list(color="blue"),
  #                    mode="markers",
  #                    hoverinfo = 'text',
  #                    text = paste0("Obs:     ",
  #                                  1:nrow(res$y),
  #                                  "<br>Classif.:",
  #                                  apply(res$logtik.obs, 1, which.max),
  #                                  "<br>Prob. error:",
  #                                  round(1-apply(exp(res$logtik.obs), 1, max), 3)))
  #
  #   if (add.modes)
  #     p <- add_trace(p,
  #                    x=res$modes[,dim[1]],
  #                    y=res$modes[,dim[2]],
  #                    type="scatter",
  #                    marker=list(color="red"),
  #                    mode="markers")
  #
  #   p
  #
  # }

#
#   plotDensityClusVisu <- function(res,
#                                   add.centers=FALSE,
#                                   threshold=.8,
#                                   dim=c(1,2),
#                                   col.dots=res$partition$hard,
#                                   col.centers=1,
#                                   cex.dots=0.8,
#                                   lwd.dots=0.9,...){
#     n <- 4 * (10**2)
#     rangey1 <- range(c(res$centers[,dim[1]], res$y[,dim[1]])) + c(-2, 2)
#     rangey2 <- range(c(res$centers[,dim[2]], res$y[,dim[2]])) + c(-2, 2)
#     input <- list(x=(seq(rangey1[1], rangey1[2], length.res = n)), y=(seq(rangey2[1], rangey2[2], length.res = n)))
#     input$z <- reser(input$x,
#                      input$y,
#                      function(x, y, centers, prop){
#                        res <- matrix(sapply(1:nrow(centers),
#                                             function(k) dnorm(x, centers[k,dim[1]], log=TRUE) + dnorm(y, centers[k,dim[2]], log=TRUE) + log(prop[k])), ncol=nrow(centers))
#                        res <- exp(sweep(res, 1, apply(res, 1, max), "-"))
#                        res <- sweep(res, 1, rowSums(res), "/")
#
#                        return(apply(res, 1, max))
#                      },
#                      centers=res$centers, prop=res$prop)
#     threshold <- sort(threshold, decreasing = F)
#     input$zthreshold <- input$z * 0
#     for (h in 1:length(threshold))  input$zthreshold <- input$zthreshold + (input$z>threshold[h])
#     image(input$x, input$y, input$zthreshold, col=c("white",paste0("gray", floor(50 + (1-threshold)*50))),
#           xlab=paste0("Dim1 (", 100*round(res$inertia[dim[1]]/sum(res$inertia), 4), "%)") ,
#           ylab=paste0("Dim2 (", 100*round(res$inertia[dim[2]]/sum(res$inertia), 4), "%)", ...))
#     xlab="dim1",
#     ylab="dim2", ...)
# maxtik <- exp(apply(res$logtik.obs, 1, max))
# maxtik <- (maxtik - 1/ncol(res$logtik.obs)) / (1-1/ncol(res$logtik.obs))
# colo <- rep("gray10", length(maxtik))
# if (any(maxtik>threshold)) colo[which(maxtik > threshold)] <- "gray70"
# plot(res$y, pch=20,#14+apply(res$logtik.obs, 1, which.max),
#      cex=cex.dots,
#      lwd=lwd.dots,
#      col=  colo,
#      xlab=paste0("Dim1 (", 100*round(res$inertia[dim[1]]/sum(res$inertia), 4), "%)") ,
#      ylab=paste0("Dim2 (", 100*round(res$inertia[dim[2]]/sum(res$inertia), 4), "%)"), ...)
# res$centersy <- res$centers * 0
# for (k in 1:nrow(res$centers)){
#   res$centersy[k,] <- t(res$y) %*% as.numeric(exp(res$logtik.obs[,k])) / sum(exp(res$logtik.obs[,k]))
# }
# points(res$centersy[,1], res$centersy[,2], pch=22, col="black", bg="white", cex=1.8)
# text(res$centersy[,1], res$centersy[,2], 1:nrow(res$centers), cex=0.7)
#contour(input, levels = threshold, add=TRUE, labcex=.8)
#if (add.centers) points(res$centers, pch=14 + (1:ncol(res$logtik.obs)), cex=2, col=col.centers, lwd=2)

# if (add.obs) points(res$y, pch=14+apply(res$logtik.obs, 1, which.max), cex=cex.dots,
#                     col=  paste0("gray", ceiling( .8*(100 - maxtik*100))))


