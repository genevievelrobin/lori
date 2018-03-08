\name{aravo}
\alias{aravo}
\docType{data}
\title{Aravo
}
\description{Count table of the abundances of 82 plant species across 75 environments in the French Alps}
\usage{data("aravo")}
\format{
  The format is:
List of 4
 $ spe      :'data.frame':	75 obs. of  82 variables:
 $ env      :'data.frame':	75 obs. of  6 variables:
 $ traits   :'data.frame':	82 obs. of  8 variables:
 $ spe.names: Named chr [1:82] "Agrostis rupestris All." "Alopecurus alpinus Vill." "Anthoxanthum odoratum L.  subsp. nipponicum (Honda) Tzvelev" "Helictotrichon sedenense (DC.) Holub" ...
  ..- attr(*, "names")= chr [1:82] "Agro.rupe" "Alop.alpi" "Anth.nipp" "Heli.sede" ...
}

\source{Dray, S. and Dufour, A.B. (2007): The ade4 package: implementing the duality diagram for ecologists. Journal of Statistical Software. 22(4): 1-20.
}
\references{
 Chessel, D. and Dufour, A.B. and Thioulouse, J. (2004): The ade4 package-I- One-table methods. R News. 4: 5-10.,
 Dray, S. and Dufour, A.B. and Chessel, D. (2007): The ade4 package-II: Two-table and K-table methods. R News. 7(2): 47-52.
}
\examples{
data(aravo)
}
\keyword{datasets}
