#' @title
#' Multi-Criteria Method for Decision (TOPSIS)
#'
#' @description
#' Assists in the TOPSIS analysis process, designed to return at the end of the answer of the TOPSIS multicriteria analysis,  a ranking table with the best option as the analysis proposes. TOPSIS is basically a technique developed by Hwang and Yoon in 1981, starting from the point that the best alternative should be closest to the positive ideal solution and farthest from the negative one, based on several criteria to result in the best benefit. (LIU, H. et al., 2019) <doi:10.1016/j.agwat.2019.105787>.
#'
#'
#'
#' @return returns a table with the results of TOPSIS.
#'
#' @author Wagner Martins dos Santos
#'
#' @usage
#' TopSisWM(TabTo,P.Pesos,solpo)
#'
#' @param TabTo A table or matrix containing row and column names, containing only numeric values.
#' @param P.Pesos A vector containing the weights assigned to each criterion, where there must be a weight for each column of the table and the sum of the weights must be equal to 1.
#' @param solpo A vector containing the information if the respective criterion must be ideal, a maximum or minimum value, with ("+") for maximum and ("-") for minimum, must have a definition for each column.
#'
#' @export
#' TopSisWM
#' @examples
#'data(Topex1)
#'Tabela <- Topex1
#'Pesos <- c(0.5,0.3,0.2)
#'solucao <- c("-","-","+")
#'TopSisWM(Tabela,Pesos,solucao)
#'
#' #Important to use the parameter row.names = 1 when importing local data
#' #Tabela <- read.csv("file.csv",sep=";", h=T,row.names=1)
#' #Tabela <- read.table("Tabela.txt", h=T,row.names=1)
#' #row.names = 1 ensures import of data with row names as fact names.
#'
#'@references
#'LIU, H. et al.(2019). Optimizing irrigation frequency and amount to balance yield, fruit quality and water use efficiency of greenhouse tomato. Agricultural Water Management, v. 226, 1-11.
#'
#'



TopSisWM <- function(TabTo=NULL,P.Pesos=NULL,solpo=NULL){
  e <- NULL
  as.matrix(TabTo)
  if (length(solpo) != ncol(TabTo)){
    stop("The number of positive solutions differs from the number of criteria")
  }
  if (length(P.Pesos) == 0){
    P.Pesos <- rep(round(1/ncol(TabTo),2),ncol(TabTo))
  }
  if (length(P.Pesos)<ncol(TabTo)){
    stop("The number of weights differs from the number of criteria")
  }
  if (sum(P.Pesos)<1){
    stop("The total sum of the weights must be equal to 1")
  }
  if (sum(P.Pesos)>1){
    stop("The total sum of the weights must be equal to 1")
  }


  VerSolP <- c("\u02B","\u02D")

  if (length(setdiff(solpo,VerSolP)) != 0){
    stop("The arguments allowed in the positive solution are:
       \u02D when the ideal value is the lowest
       \u02B when the ideal value is the highest")
  }

  tryCatch({

    norm<-sqrt(colSums(TabTo^2))
    TaN <- t(t(TabTo) / norm)
    TaNP <- t(t(TaN) * P.Pesos)
    solne <- solpo
    TaSP <- matrix(solpo,nrow=1)
    TaSN <- matrix(solne,nrow=1)
    SolP = c()
    SolN = c()
    x=0

    for (i in solpo){
      x = x+1
      colg <- TaNP[,x]
      SolP[x][TaSP[x] == "\u02B"] <- max(colg)
      SolP[x][TaSP[x] == "\u02D"] <- min(colg)
    }

    x = 0

    for (i in solne){
      x = x+1
      colg <- TaNP[,x]
      SolN[x][TaSN[x] == "\u02B"] <- min(colg)
      SolN[x][TaSN[x] == "\u02D"] <- max(colg)
    }

    x=0
    dip = c()
    for (i in SolP){
      x = x+1
      try(dip[x] <- (sqrt(sum((TaNP[x,]-SolP)^2))))
    }


    x=0
    din = c()
    for (i in SolN){
      x = x+1
      try(din[x] <- (sqrt(sum((TaNP[x,]-SolN)^2))))
    }

    MD <- matrix(c(dip,din),ncol=2)

    x=0
    Prox = c()
    for (i in rep(nrow(MD),nrow(MD))){
      x = x+1
      Prox[x] <- (MD[x,2]/(MD[x,1]+MD[x,2]))
    }


    rank = rank(-Prox)
    Prox <- matrix(c(dip,din,Prox,rank),ncol=4)

    row.names(Prox) <- row.names(TabTo)
    colnames(Prox) <- c("d+","d-","eu","rank")
    lista = list(TaN,TaNP,Prox)
    names(lista) <- c("Normalized Table",
                      "Weighted Normalized Table",
                      "Decision Table - TOPSIS")

    return(Prox)

  },
  warning=function(war)
  {
    print(paste("WARNING: ", war))
  },
  error=function(err)
  {
    print(paste("ERROR: ", err))
  },
  finally=function(f)
  {
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
