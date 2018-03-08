#' Plot the rows and columns of the contingency table in the Euclidean space defined by two of the
#' first principal directions of the interaction matrix Theta. Theta corresponds to the interaction
#' remaining after discarding the effects of the covariates. The interpretation is the following.
#' A row and a column that are close in Euclidean distance interact highly. Two rows or two columns that
#' are close in Euclidean distance have similar profiles.
#' @param theta An interaction matrix ($Theta in output from the application of the lori function).  Theta can have extra columns corresponding to supplementary quantitative variables
#' @param quanti.sup a vector of column indices giving the column indices of the supplementary quantitative variables if there are some.
#' @param axes a vector of integers of length two, indicating which axes to plot (c(1,2) means that the first two directions are plotted).
#' @param xlim a vector of length two indicating the limits of the plot abscissa.
#' @param ylim a vector of length two indicating the limits of the plot ordinate.
#' @param invisible string indicating if some points should be unlabelled ("row", "col", "row.sup", "col.sup","quali.sup")
#' @param choix the graph to plot ("CA" for the CA map, "quanti.sup" for the supplementary quantitative variables)
#' @param col.row a color for the rows points
#' @param col.col a color for columns points
#' @param col.row.sup a color for the supplementary rows points
#' @param col.col.sup a color for supplementary columns points
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
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y <- matrix(rpois(length(c(X)), exp(c(X))), 5)
#' res_lori <- lori(Y)
#' p <- plot_interaction(res_lori$Theta)
plot_interaction <- function(theta, quanti.sup = NULL, axes = c(1, 2), xlim = NULL, ylim = NULL,
                             invisible = c("none","row", "col", "row.sup", "col.sup", "quali.sup"),
                             choix = "CA", col.row = "blue", col.col = "red",
                             col.row.sup = "darkblue", col.col.sup = "darkred",
                             col.quanti.sup = "blue",label = c("all", "none", "row", "row.sup", "col",
                                                               "col.sup", "quali.sup"), title = NULL,
                             palette = NULL, autoLab = c("auto", "yes", "no"), new.plot = FALSE,
                             selectRow = NULL, selectCol = NULL, unselect = 0.7, shadowtext = FALSE,
                             habillage = "none", legend = list(bty = "y",  x = "topleft"), ...){
  duv <- svd(theta[, setdiff(1:ncol(theta), quanti.sup)])
  ## Multiply interaction by singular values
  u <- duv$u%*%diag(sqrt(duv$d))
  v <- duv$v%*%diag(sqrt(duv$d))
  theta[,setdiff(1:ncol(theta), quanti.sup)] <- theta[,setdiff(1:ncol(theta), quanti.sup)]+2
  res_ca <- FactoMineR::CA(theta, graph = F, quanti.sup = quanti.sup)
  res_ca$row$coord <- u
  rownames(res_ca$row$coord) <- rownames(theta)
  res_ca$col$coord <- v
  rownames(res_ca$col$coord) <- colnames(theta)[setdiff(1:ncol(theta), quanti.sup)]
  if(!is.null(quanti.sup)){
    res_ca$quanti.sup$coord <- t(theta[, quanti.sup])%*%duv$u
    norms <- sqrt(rowSums(res_ca$quanti.sup$coord^2))
    res_ca$quanti.sup$coord <- sweep(res_ca$quanti.sup$coord, 1, norms, "/")

  }
  if(is.null(title)){
    if(choix == "CA") {
      title <- "2D Display plot of interaction directions"
    } else if (choix == "quanti.sup") title <- "Supplementary variables on the interaction map"
  }
  FactoMineR::plot.CA(res_ca, axes = axes, xlim = xlim, ylim = ylim, invisible = invisible,
                      choix = choix, col.row = col.row, col.col = col.col, col.row.sup = "darkblue",
       col.col.sup = col.col.sup, col.quanti.sup = col.quanti.sup, label = label, title = title,
       palette = palette, autoLab = autoLab, new.plot = new.plot, selectRow = selectRow,
       selectCol = selectCol, unselect = unselect, shadowtext = shadowtext, habillage = habillage,
       legend = legend)
}
