#' @title
#' Multi-Criteria Method for Decision (TOPSIS)
#'
#' @description
#' Assists in the TOPSIS analysis process, designed to return a ranking table with the best option as the analysis proposes. TOPSIS is a technique developed by Hwang and Yoon (1981), assuming the best alternative is closest to the positive ideal solution and farthest from the negative one. (LIU, H. et al., 2019) <doi:10.1016/j.agwat.2019.105787>.
#'
#' @param TabTo A table or matrix containing row and column names, with only numeric values. Default is \code{NULL}.
#' @param P.Pesos A vector of weights for each criterion. Must have a weight for each column of the table. The sum of weights must be 1. Default is \code{NULL}.
#' @param solpo A vector indicating if each criterion is to be maximized ("+") or minimized ("-"). Must match the number of columns. Default is \code{NULL}.
#' @param sort.rank Logical. If \code{TRUE}, the result is ordered by rank. If \code{FALSE}, keeps original input order. Default is \code{TRUE}.
#'
#' @return A matrix with the results of the TOPSIS analysis, including the distance to the positive ideal solution (d+), the distance to the negative ideal solution (d-), the relative closeness (eu), and the final rank.
#'
#'
#' @export
#'
#' @examples
#' data(Topex1)
#' Tabela <- Topex1
#' Pesos <- c(0.5,0.3,0.2)
#' solucao <- c("-","-","+")
#' TopSisWM(Tabela, Pesos, solucao)
#'
#' # Important: use row.names = 1 when importing local data
#' # Tabela <- read.csv("file.csv", sep=";", h=TRUE, row.names=1)
#' # Tabela <- read.table("Tabela.txt", h=TRUE, row.names=1)
#'
#' @references
#' LIU, H. et al. (2019). Optimizing irrigation frequency and amount to balance yield, fruit quality and water use efficiency of greenhouse tomato. Agricultural Water Management, 226, 1–11.




TopSisWM <- function (TabTo = NULL, P.Pesos = NULL, solpo = NULL,sort.rank = TRUE)
{
  e <- NULL
  TabTo <- as.matrix(TabTo)
  if (length(solpo) != ncol(TabTo)) {
    stop("The number of positive solutions differs from the number of criteria")
  }
  if (length(P.Pesos) == 0) {
    P.Pesos <- rep(round(1/ncol(TabTo), 2), ncol(TabTo))
  }
  if (length(P.Pesos) < ncol(TabTo)) {
    stop("The number of weights differs from the number of criteria")
  }
  if (round(sum(P.Pesos), 1) < 1) {
    stop("The total sum of the weights must be equal to 1")
  }
  if (sum(P.Pesos) > 1) {
    stop("The total sum of the weights must be equal to 1")
  }
  VerSolP <- c("+", "-")
  if (length(setdiff(solpo, VerSolP)) != 0) {
    stop("The arguments allowed in the positive solution are:\n       - when the ideal value is the lowest\n       + when the ideal value is the highest")
  }
  tryCatch({
    norm <- sqrt(colSums(TabTo^2))
    TaN <- t(t(TabTo)/norm)
    TaNP <- t(t(TaN) * P.Pesos)
    solne <- solpo
    TaSP <- matrix(solpo, nrow = 1)
    TaSN <- matrix(solne, nrow = 1)
    SolP = c()
    SolN = c()
    x = 0
    for (i in solpo) {
      x = x + 1
      colg <- TaNP[, x]
      SolP[x][TaSP[x] == "+"] <- max(colg)
      SolP[x][TaSP[x] == "-"] <- min(colg)
    }
    x = 0
    for (i in solne) {
      x = x + 1
      colg <- TaNP[, x]
      SolN[x][TaSN[x] == "+"] <- min(colg)
      SolN[x][TaSN[x] == "-"] <- max(colg)
    }
    x = 0
    dip = c()
    for (i in 1:nrow(TaNP)) {
      x = x + 1
      try(dip[x] <- (sqrt(sum((TaNP[x, ] - SolP)^2))))
    }
    x = 0
    din = c()
    for (i in 1:nrow(TaNP)) {
      x = x + 1
      try(din[x] <- (sqrt(sum((TaNP[x, ] - SolN)^2))))
    }
    MD <- matrix(c(dip, din), ncol = 2)
    x = 0
    Prox = c()
    for (i in rep(nrow(MD), nrow(MD))) {
      x = x + 1
      Prox[x] <- (MD[x, 2]/(MD[x, 1] + MD[x, 2]))
    }
    rank = rank(-Prox)
    Prox <- matrix(c(dip, din, Prox, rank), ncol = 4)
    row.names(Prox) <- row.names(TabTo)
    colnames(Prox) <- c("d+", "d-", "eu", "rank")
    lista = list(TaN, TaNP, Prox)
    names(lista) <- c("Normalized Table", "Weighted Normalized Table",
                      "Decision Table - TOPSIS")

    if(sort.rank){Prox <- Prox[order(Prox[, 4]), ]}

    return(Prox)

  }, warning = function(war) {
    print(paste("WARNING: ", war))
  }, error = function(err) {
    print(paste("ERROR: ", err))
  }, finally = function(f) {
    return(print(paste("e: ", e)))
  })
}



#'@title
#'Topex1 - TopSisWM
#'
#'@description
#'Example Data
#'
#'@format Matrix
#'
#'@source Topex1
"Topex1"
