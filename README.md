# orfr
## Model-Based Calibration and Scoring for Oral Reading Fluency Assessment Data with R

`orfr` is an R package that allows model-based calibration and scoring for oral reading fluency (ORF) assessment data.

**Installation:**  
To use `orfr`, you need a `MultiGHQuad` package installed. As of May 31, 2022, `MultiGHQuad` package has been removed from CRAN, so you need to download a binary package file from the archive https://cran.r-project.org/src/contrib/Archive/MultiGHQuad/ and install.

To install `orf` package, follow the steps.
1. Install `remotes` package, if you have not.
if(!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
2. Install `orfr`
remotes::install_github("kamataak/orfr")

**Usage:**  



Copyright (C) 2022 The ORF Project Team

The orfr package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
 any later version.

The orfr package is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this package.  If not, see <http://www.gnu.org/licenses/>.
