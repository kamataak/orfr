#' This file includes the log functions of orfr package.
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
#' This is log.initiating function, which is used to log output.
#'
#' @import tryCatchLog
#' @import futile.logger
#'
#'
#'
log.initiating <- function() {

options(keep.source = TRUE) # source code file name and line number tracking
options("tryCatchLog.write.error.dump.file" = FALSE)
options("tryCatchLog.include.full.call.stack" = FALSE) # reduce the ouput for demo purposes
options("tryCatchLog.include.compact.call.stack" = FALSE) # reduce the ouput

# Set log level
flog.threshold(INFO) # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

# Define logger, write to both console and file
flog.logger("orfrlog", INFO, appender=appender.tee('orfrlog.log'))

}
