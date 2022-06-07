orfr
================

## Model-Based Calibration and Scoring for Oral Reading Fluency Assessment Data with R

`orfr` is an R package that allows model-based calibration and scoring
for oral reading fluency (ORF) assessment data.

### Installation:

To install `orf` package, follow the steps below.

1.  The `MultiGHQuad` package needs to be installed. As of May 31, 2022,
    `MultiGHQuad` package has been removed from CRAN, so download an
    archive package file from
    <https://cran.r-project.org/src/contrib/Archive/MultiGHQuad/> and
    manually install.

2.  Install `remotes` package, by running the following R code.

``` r
if(!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
```

3.  Install `orfr` by running the following R code.

``` r
remotes::install_github("kamataak/orfr")
```

### Basic Usage:

It is recommended that the data are prepared as a long-format data
frame, where each row is data for a unique case, namely, a specific
passage from a specific student in a specific testing occasion. The data
set should minimally contain the following variables: (1) student ID,
(2) grade level, (3) testing occasion ID, (4) passage ID, (5) the number
of words in the passage, (6) the number of words correctly read for the
passage, and (7) time that took to read the passage.

Variable names and the order of the variables can be flexible in the
data frame. When running functions in this package, variable names for
required variables need to be specified.

To demonstrate some basic usage of key functions in the package, the
`orfr` package contains an example passage-level student data set
`passage`.

Load required packages, and load/view the example data set `passage`.

``` r
library(tidyverse)
library(orfr)
data(passage)
View(passage)
```

#### Passage Calibration

Calibrate the passages using the `mcem()` function.

``` r
test_MCEM <- mcem(passage,
                  studentid = "id.student",
                  passageid = "id.passage",
                  numwords.p = "numwords.pass",
                  wrc = "wrc",
                  time = "sec",
                  k.in = 5,
                  reps.in = 2,
                  est = "mcem")
test_MCEM
```

#### Estimating WCPM scores

To estimate WCPM scores, we can do in three steps.

**Step 1:** Prepare the data using the `preplong()` function, where
required variables for the `wcpm()` function are prepared, including
changing variable names and a generation of the natural-logarithm of the
time data. In addition, we can use another utility function
`get.cases()` to generate a list of unique cases with the output of the
`preplong()` function.

``` r
datalong <- preplong(data = passage,
                     studentid = "id.student",
                     passageid = "id.passage",
                     season = "occasion",
                     grade = "grade",
                     numwords.p = "numwords.pass",
                     wrc = "wrc",
                     time = "sec")
```

Generate a list of unique cases:

``` r
get.cases(datalong)
```

**Step 2:** Specify a list of cases for which WCPM scores are estimated.
It has to be a one-variable data frame with a variable name `cases`.

The format of case values should be: `studentid_season`.

``` r
sample.cases <- data.frame(cases = c("2056_fall", "2056_winter", "2056_spring"))
```

If cases are not specified, the `wcpm()` function will estimate WCPM
scores for all cases in the data in the next step.

**Step 3:** Run the `wcpm()` function to estimate WCPM scores for
selected cases. Note that we pass the output object `test_MCEM` from the
passage calibration step, as well as the manipulated data `datalong`
from Step 1. Also, if the `cases =` argument is omitted, WCPM scores
will be estimated for all cases in the data.

``` r
test_WCPMEAP <- wcpm(test_MCEM, 
                     stu.data = datalong,
                     cases = sample.cases, 
                     est = "eap", 
                     se = "analytic")
summary(test_WCPMEAP)
```

Also, we can specify a set of passages to scale the WCPM scores. If WCPM
scores are scaled with a set of passages that is different from the set
of passages the student read, the set of passages is referred to as an
**external passage set**.

The use of an external passage set is particularly important to make the
estimated WCPM scores to be comparable between students for
cross-sectional data, as well as within students for longitudinal data.

``` r
test_WCPMEAP_EXT <- wcpm(test_MCEM, 
                     stu.data = datalong,
                     cases = sample.cases, 
                     external = c("22007","22013","22036","22043","22048","22079"),
                     est = "eap", 
                     se = "analytic")
summary(test_WCPMEAP_EXT)
```

**\[*Not yet implemented*\]** Alternatively, we can run the `wcpm()`
function without **Step 1** above, by entering the original data
`passage` directly as follows.

``` r
test_WCPMEAP <- wcpm(test_MCEM, 
                     stu.data = passage,
                     studentid = "id.student",
                     passageid = "id.passage",
                     season = "occasion",
                     grade = "grade",
                     numwords.p = "numwords.pass",
                     wrc = "wrc",
                     time = "sec",
                     cases = sample.cases, 
                     external = c("22007","22013","22036","22043","22048","22079"),
                     est = "eap", 
                     se = "analytic")
summary(test_WCPMEAP)
```

The package also contains calibrated passage data `MCME`, in which
passages were calibrated with much larger sample.

Please see the [package website](https://kamataak.github.io/orfr/) for
more detailed usage of the package.

### Citation

### Copyright Statement

Copyright (C) 2022 The ORF Project Team

The orfr package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later
version.

The orfr package is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details. You should have received a copy of the
GNU General Public License along with this package. If not, see
<http://www.gnu.org/licenses/>.
