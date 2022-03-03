#' Plot function to show graph of mcem class
#'
#'
#' Copyright (C) 2021 The ORF Project Team
#'
#' This program is free software; you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation; either version 3 of the License, or
#' (at your option) any later version.
#
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#
#' A copy of the GNU General Public License is available at
#' http://www.gnu.org/licenses/
#'
#' @param object = mcem object
#' @param X = dataset to plot
#' @param ... = parameter
#'
#' @import psych
#'
#' @method plot mcem
#' @export
plot.mcem <- function(object,X){
  data <- X
  plot.it <- data.frame(Score=colMeans(data$Y, na.rm=TRUE)/data$N*10,
                        Time.secs=colMeans(exp(data$logT), na.rm=TRUE)/data$N*10)
  pairs.panels(plot.it,
               method = "pearson",
               hist.col = "#00AFBB",
               cex.cor=0.5)
}
