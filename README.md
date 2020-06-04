# EITrans: Empirical Inverse Transform Function

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

<!-- vim-markdown-toc GitLab -->

* [Overview](#overview)
* [Installation](#installation)

<!-- vim-markdown-toc -->

## Overview

This package proposes an algorithm that use [weather analogs](https://weiming-hu.github.io/AnalogsEnsemble/) for ensemble forecast calibration.

The package is currently under-developing and relevant research is continuing.

Please visit [Analog Ensembles](https://weiming-hu.github.io/AnalogsEnsemble/) to learn more about weather analogs.

Please visit [my website](https://weiming-hu.github.io/) to learn more about my research.


## Installation

To install `EITrans`, first install the following dependent packages:

```
# The repos argument is just to skip the repo selection.
# You can simply omit the argument if you want.
#
install.packages(c('abind', 'foreach', 'progress'), repos = 'http://cran.us.r-project.org')

# More installation on the following package, please refer to
# https://weiming-hu.github.io/AnalogsEnsemble/
#
install.packages("https://github.com/Weiming-Hu/AnalogsEnsemble/raw/master/RAnalogs/releases/RAnEn_latest.tar.gz", repos = NULL)

devtools::install_github('Weiming-Hu/RAnEnExtra')
```

Then, you can install `EITrans` by running

```
devtools::install_github('Weiming-Hu/EITrans')
```
