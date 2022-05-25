#' summary the information of mcem class
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
#' @param object = object
#' @param digits = print out numeric with specific digits
#' @param ... = parameter
#'
#' @import mvtnorm
#' @import tidyverse
#' @return printing information
#'
#' @method summary mcem
#' @export
summary.mcem <- function(object, digits=4,...) {
  z <- object
  tab <- cbind(z$pass.param$passage_id,
               as.vector(z$pass.param$a),
               as.vector(z$pass.param$b),
               z$pass.param$alpha,
               z$pass.param$beta)
  colnames(tab)=c("Passage_id","a","      b","   alpha","   beta")
  rownames(tab)=paste0(c(rep(1:length(tab[,1]))),".")
  #  print(tab[, 1:3]) # only print a part of columns
  #  print(tab)
  print(tab, digits = digits, print.gap = 2L) # specific minimum digits
  cat("\n====== Hyper Parameters ======\n")
  cat(paste(paste0("Variance of ", greek$tau), ":     "))
  cat(paste(format(z$hyper.param$vartau,digits=6,nsmall=digits), "\n"))
  cat(paste(greek$rho), "            :     ")
  cat(paste(format(z$hyper.param$rho,digits=6,nsmall=digits), "\n"))
}
#' summary the information of wcpm class
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
#' @param object = object
#' @param digits = print out numeric with specific digits
#' @param verbose - boolean, if TRUE, shows the summary, default is TRUE
#' @param factor.scores - theta and tau output flag, default is FALSE
#'
#' @return wcpm dataset with passage information and estimated score
#' @method summary wcpm
#' @export
summary.wcpm <- function(object, digits=4,verbose=TRUE,factor.scores=FALSE) {
  z <- object

  tb <- as.data.frame(t(do.call(rbind, z)))

  # don't output theta and tau, if FALSE
  if (factor.scores==FALSE) {
    tb <- tb %>% select(-contains(c("tau", "theta")))
  }

  getNames <- colnames(tb)

  cols_num <- ncol(tb)
  #set screen print out to be short decimal
  tt <- as.matrix(unlist(lapply(as.double(unlist((tb[,c(6:cols_num)]))),
                                sprintf, fmt = "%6.2f")))
  dim(tt) <- c(dim(tb)[1],(cols_num-5))
  tt <- cbind(tb[,c(1:5)], tt)
  colnames(tt) <- getNames
  if (verbose == TRUE) {
    # only verbose TRUE will print out on screen
    print(tt)
    return(invisible(tb))
  } else {
    return(invisible(tb))
  }

}
#' summary the information of bootstrap class
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
#' @param object = bootstrap object
#' @param digits = print out numeric with specific digits
#' @param geterror, summary error case, default FALSE
#' @param verbose show summary on screen, default TRUE
#'
#' @return
#' @export
summary.bootstrap <- function(object, digits=4, geterror=FALSE,verbose=TRUE) {
  z <- object

  tb <- z$bootstrap.out
  if (geterror == TRUE) {
    if (length(z$error_case) != 0) {
      print(z$error_case)
      return(invisible(z$error_case))
    } else {
      print("Bootstrap has no error cases.")
    }
  } else {
    if (ncol(tb) != 0) {
      getNames <- colnames(tb)
      cols_num <- ncol(tb)

      #set screen print out to be short decimal
      tt <- as.matrix(unlist(lapply(as.double(unlist((tb[,c(6:cols_num)]))),
                                    sprintf, fmt = "%6.2f")))
      dim(tt) <- c(nrow(tb),(cols_num-5))
      tt <- cbind(tb[,c(1:5)], tt)
      colnames(tt) <- getNames
      if (verbose == TRUE) {
        # only verbose TRUE will print out on screen
        print.noquote(tt)
        return(invisible(tb))
      } else {
        return(invisible(tb))
      }
    } else {
      print("Bootstrap has 0 obs.")
    }
  }
}
