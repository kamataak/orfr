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

The `orfr` package contains an example passage-level student data set
`passage`, which is stored in the `passage.rda` file, to demonstrate
some basic usage of the package.

#### Passage Calibration

Load the example data set `passage`.

``` r
load("passage.rda")
```

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
                  est="mcem")
summary(test_MCEM)
```

#### Estimating WCPM scores

To estimate WCPM scores, we can do in three steps.

**Step 1:** Prepare the data using the `preplong()` function, where
required variables to run the `wcpm()` function are prepared, including
the change of variable names and a generation of the natural-logarithm
of the time data. In addition, we can use another utility function
`get.cases()` to generate a list of unique cases.

``` r
datalong <- preplong(stu.data = passage,
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

**Step 2:** Specify a list of cases for which WCPM scores are estimated:

``` r
sample.cases <- data.frame(cases = c("2056_fall", "2056_winter", "2056_spring"))
```

**Step 3:** Run the `wcpm()` function to estimate WCPM scores for
selected cases. Here, note that we enter the manipulated data `datalong`
from the Step 1.

``` r
test_WCPMEAP <- wcpm(test_MCEM, 
                     stu.data=datalong,
                     cases=sample.cases, 
                     est="eap", 
                     se="analytic"))
summary(test_WCPMEAP)
```

Alternatively, we can run the `wcpm()` function without **Step 1**
above, by entering the original data `passage` directly as follows.

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
                     cases=sample.cases, 
                     est="eap", 
                     se="analytic")
summary(test_WCPMEAP)
```

The package also contains calibrated passage data `MCME.rda`, in which
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
