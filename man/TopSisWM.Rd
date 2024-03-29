% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopSisWM.R
\name{TopSisWM}
\alias{TopSisWM}
\title{Multi-Criteria Method for Decision (TOPSIS)}
\usage{
TopSisWM(TabTo,P.Pesos,solpo)
}
\arguments{
\item{TabTo}{A table or matrix containing row and column names, containing only numeric values.}

\item{P.Pesos}{A vector containing the weights assigned to each criterion, where there must be a weight for each column of the table and the sum of the weights must be equal to 1.}

\item{solpo}{A vector containing the information if the respective criterion must be ideal, a maximum or minimum value, with ("+") for maximum and ("-") for minimum, must have a definition for each column.}
}
\value{
returns a table with the results of TOPSIS.
}
\description{
Assists in the TOPSIS analysis process, designed to return at the end of the answer of the TOPSIS multicriteria analysis,  a ranking table with the best option as the analysis proposes. TOPSIS is basically a technique developed by Hwang and Yoon in 1981, starting from the point that the best alternative should be closest to the positive ideal solution and farthest from the negative one, based on several criteria to result in the best benefit. (LIU, H. et al., 2019) \url{doi:10.1016/j.agwat.2019.105787}.
}
\examples{
data(Topex1)
Tabela <- Topex1
Pesos <- c(0.5,0.3,0.2)
solucao <- c("-","-","+")
TopSisWM(Tabela,Pesos,solucao)

#Important to use the parameter row.names = 1 when importing local data
#Tabela <- read.csv("file.csv",sep=";", h=T,row.names=1)
#Tabela <- read.table("Tabela.txt", h=T,row.names=1)
#row.names = 1 ensures import of data with row names as fact names.

}
\references{
LIU, H. et al.(2019). Optimizing irrigation frequency and amount to balance yield, fruit quality and water use efficiency of greenhouse tomato. Agricultural Water Management, v. 226, 1-11.
}
\author{
Wagner Martins dos Santos
}
