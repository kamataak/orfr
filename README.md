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
set should minimally contain the following 7 variables: (1) student ID,
(2) grade level, (3) testing occasion ID, (4) passage ID, (5) the number
of words in the passage, (6) the number of words correctly read for the
passage, and (7) time that took to read the passage.

Variable names and the order of the variables can be flexible in the
data frame. When running functions in this package, variable names for
required variables need to be specified.

The `orfr` package comes with several data sets. To demonstrate some
basic usage of key functions in the package, an example passage-level
student data set `passage2` is used here. The data set `passage2` is
consisted of reading accuracy and time data for 12 passages from 85
students. Although the 85 students were assigned to all 12 passages, the
number of passages read by the 85 students varied from 2 to 12 passages.
The number of students per passage were between 59 to 79.

Load required packages, and load/view the example data set `passage`.

``` r
library(tidyverse)
library(orfr)
View(passage2)
```

#### Passage Calibration

Calibrate the passages using the `mcem()` function.

``` r
system.time(
MCEM_run <- mcem(stu.data=passage2,
                 studentid = "id.student",
                 passageid = "id.passage",
                 numwords.p = "numwords.pass",
                 wrc = "wrc",
                 time = "sec",
                 k.in = 5,
                 reps.in = 50,
                 est = "mcem")
)
MCEM_run
```

By default, the standard errors for the model parameters are not
estimated. This will allow one to increase the number of Monte-Carlo
iteration `reps.in` to improve the quality of the estimates of the model
parameters, while minimizing the computation time. The number of
`reps.in` should be 50-100 in realistic calibrations. SE’s for model
parameters are not required for running `wcpm()` function to estimate
WCPM scores in the next step. In order to compute the standard errors
for the model parameters, an additional argument `se = "analytical"` or
`se = "bootstrap"` needs to be added to the `mcem()` function.

#### Estimating WCPM scores 1

To estimate WCPM scores, we can do in two steps.

**Step 1:** Prepare the data using the `prep()` function, where required
data for the `wcpm()` function are prepared, including changing variable
names and a generation of the natural-logarithm of the time data.

The output from the `prep()` function is a list of two components. The
`data.long` component is a data frame, which is a long format of student
response data, and the `data.wide` is list that contains four
components, including a wide format of the data, as well as other
information such as the number of passages and the number of words for
each passage.

One benefit of this two-step approach is that we can use another utility
function `get.cases()` to generate a list of unique cases with the
output of the `prep()` function. This list can be useful when our
interest is to estimate WCPM scores only for selected cases.

``` r
data <- prep(data = passage2,
             studentid = "id.student",
             season = "occasion",
             grade = "grade",
             passageid = "id.passage",
             numwords.p = "numwords.pass",
             wrc = "wrc",
             time = "sec")
```

Generate a list of unique cases:

``` r
get.cases(data$data.long)
```

**Step 2:** Run the `wcpm()` function to estimate WCPM scores. Note that
we pass the output object `MCEM_run` from the passage calibration phase,
as well as the manipulated data `data.long` from Step 1. By default,
WCPM scores will be estimated for all cases in the data. Additionally,
there are several estimator options and standard error estimation
options.

``` r
WCPM_MAP_all <- wcpm(calib.data=MCEM_run, 
                          stu.data = data$data.long,
                          est = "map", 
                          se = "analytic")
summary(WCPM_MAP_all)
```

If the computations of WCPM scores for only selected cases are desired,
we can create a list of cases and provide the list by the `cases =`
argument. The list of cases has to be a one-variable data frame with a
variable name `cases`. The format of case values should be:
`studentid_season`, just like the output of the `get.cases()` function.

``` r
sample.cases <- data.frame(cases = c("2033_fall", "2043_fall", "2089_fall"))
WCPM_MAP_sample <- wcpm(calib.data=MCEM_run, 
                             stu.data = data$data.long,
                             cases = sample.cases,
                             est = "map", 
                             se = "analytic")
summary(WCPM_MAP_sample)
```

Also, we can specify a set of passages to scale the WCPM scores. If WCPM
scores are scaled with a set of passages that is different from the set
of passages the student read, the set of passages is referred to as an
**external passage set**.

The use of an external passage set is particularly important to make the
estimated WCPM scores to be comparable between students for
cross-sectional data, as well as within students for longitudinal data.

``` r
WCPM_MAP_ext <- wcpm(calib.data=MCEM_run, 
                     stu.data = data$data.long,
                     cases = sample.cases, 
                     external = c("32004","32010","32015","32016","33003","33037"),
                     est = "map", 
                     se = "analytic")
summary(WCPM_MAP_ext)
```

#### Estimating WCPM scores 2

**Alternatively, we can run the `wcpm()` function without Step 1**
above, by entering the original data `passage2` directly as follows.

``` r
test_WCPM_MAP_ALT <- wcpm(calib.data=MCEM_run, 
                     stu.data = passage2,
                     studentid = "id.student",
                     passageid = "id.passage",
                     season = "occasion",
                     grade = "grade",
                     numwords.p = "numwords.pass",
                     wrc = "wrc",
                     time = "sec",
                     cases = sample.cases, 
                     external = c("32004","32010","32015","32016","33003","33037"),
                     est = "map", 
                     se = "analytic")
summary(test_WCPM_MAP_ALT)
```

The package also contains calibrated passage data `MCEM`, in which
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
