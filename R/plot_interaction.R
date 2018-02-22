plot_interaction <- function(theta, quanti.sup = NULL, axes = c(1, 2), xlim = NULL, ylim = NULL,
                             invisible = c("none","row", "col", "row.sup", "col.sup", "quali.sup"),
                             choix = "CA", col.row = "blue", col.col = "red",
                             col.row.sup = "darkblue", col.col.sup = "darkred", col.quali.sup = "magenta",
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
  res_ca=CA(theta, graph = F, quanti.sup = quanti.sup)
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
  plot(res_ca, axes, xlim, ylim, invisible, choix, col.row , col.col, col.row.sup = "darkblue",
       col.col.sup, col.quali.sup, col.quanti.sup, label, title, palette, autoLab, new.plot,
       selectRow, selectCol, unselect, shadowtext, habillage, legend)
}
