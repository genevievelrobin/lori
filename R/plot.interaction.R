#' Plot the rows and columns of the contingency table in the Euclidean space defined by two of the
#' first principal directions of the interaction matrix Theta. Theta corresponds to the interaction
#' remaining after discarding the effects of the covariates. The interpretation is the following.
#' A row and a column that are close in Euclidean distance interact highly. Two rows or two columns that
#' are close in Euclidean distance have similar profiles.
#' @param x a lori object resulting from the lori function
#' @param axes a vector of integers of length two, indicating which axes to plot (c(1,2) means that the first two directions are plotted).
#' @param xlim a vector of length two indicating the limits of the plot abscissa.
#' @param ylim a vector of length two indicating the limits of the plot ordinate.
#' @param invisible string indicating if some points should be unlabelled ("row", "col", "row.sup", "col.sup","quali.sup")
#' @param choix the graph to plot ("CA" for the CA map, "quanti.sup" for the supplementary quantitative variables)
#' @param col.row a color for the rows points
#' @param col.col a color for columns points
#' @param col.row.sup a color for the supplementary rows points
#' @param col.col.sup a color for supplementary columns points
#' @param col.quali.sup a color for the supplementary quanlitative variables
#' @param col.quanti.sup a color for the supplementary quantitative variables
#' @param label a list of character for the elements which are labelled (by default, all the elements are labelled ("row", "row.sup", "col", "col.sup","quali.sup")
#' @param title string corresponding to the title of the graph you draw (by default NULL and a title is chosen)
#' @param palette the color palette used to draw the points. By default colors are chosen. If you want to define the colors : palette=palette(c("black","red","blue")); or you can use: palette=palette(rainbow(30)), or in black and white for example: palette=palette(gray(seq(0,.9,len=25)))
#' @param autoLab if autoLab="auto", autoLab is equal to "yes" if there are less than 50 elements and "no" otherwise; if "yes", the labels of the drawn elements are placed in a "good" way (can be time-consuming if many elements), and if "no" the elements are placed quickly but may overlap
#' @param new.plot boolean, if TRUE, a new graphical device is created
#' @param selectRow a selection of the rows that are drawn; see the details section
#' @param selectCol a selection of the columns that are drawn; see the details section
#' @param unselect may be either a value between 0 and 1 that gives the transparency of the unselected objects (if unselect=1 the transparceny is total and the elements are not drawn, if unselect=0 the elements are drawn as usual but without any label) or may be a color (for example unselect="grey60")
#' @param shadowtext boolean; if true put a shadow on the labels (rectangles are written under the labels which may lead to difficulties to modify the graph with another program)
#' @param habillage color the individuals among a categorical variable (give the number of the categorical supplementary variable or its name)
#' @param legend a list of arguments that defines the legend if needed (when individuals are drawn according to a variable); see the arguments of the function legend
#' @param ... further arguments passed to or from other methods, such as cex, cex.main, ...
#' @method plot lori
#' @S3method plot lori
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y <- matrix(rpois(length(c(X)), exp(c(X))), 5)
#' res_lori <- lori(Y, cov=cbind(c(X),c(X)), lambda1=1, lambda2=1)
#' p <- plot(res_lori)
plot.lori <- function(x, axes = c(1, 2),
                      xlim = NULL, ylim = NULL, invisible = c("none","row", "col", "row.sup", "col.sup","quali.sup"), choix = c("CA","quanti.sup"), col.row = "blue",
                      col.col = "red", col.row.sup = "darkblue", col.col.sup = "darkred",col.quali.sup ="magenta",
                      col.quanti.sup="blue",label = c("all","none","row", "row.sup", "col","col.sup", "quali.sup"), title = NULL, palette=NULL,
                      autoLab = c("auto","yes","no"),new.plot=FALSE, selectRow = NULL, selectCol = NULL,
                      unselect = 0.7,shadowtext = FALSE, habillage = "none", legend = list(bty = "y", x = "topleft"),...){
  theta <- x$theta
  duv <- svd(theta)
  ## Multiply interaction by singular values
  u <- duv$u%*%diag(sqrt(duv$d))
  v <- duv$v%*%diag(sqrt(duv$d))
  res_ca <- FactoMineR::CA(x$imputed, graph = F)
  res_ca$row$coord <- u
  rownames(res_ca$row$coord) <- rownames(theta)
  res_ca$col$coord <- v
  rownames(res_ca$col$coord) <- colnames(theta)
  res_ca$eig[,1] <- duv$d[1:length(res_ca$eig[,1])]
  denom <- sum(res_ca$eig[,1]^2)
  if(denom==0) denom <- 1
  res_ca$eig[,2] <- res_ca$eig[,1]^2/denom
  res_ca$eig[,3] <- cumsum(res_ca$eig[,2])
  if(is.null(title)){
    title <- "2D Display plot of interaction directions"
  }
  FactoMineR::plot.CA(res_ca, axes = axes, xlim = xlim, ylim = ylim, invisible = invisible,
                      choix = choix, col.row = col.row, col.col = col.col, col.row.sup = "darkblue",
                      col.col.sup = col.col.sup, col.quanti.sup = col.quanti.sup, label = label, title = title,
                      palette = palette, autoLab = autoLab, new.plot = new.plot, selectRow = selectRow,
                      selectCol = selectCol, unselect = unselect, shadowtext = shadowtext, habillage = habillage,
                      legend = legend)
}
