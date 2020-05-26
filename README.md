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

Please visit [Analog Ensembles](weiming-hu.github.io/) to learn more about weather analogs.

Please visit [my website](weiming-hu.github.io/) to learn more about my research.


## Installation

To install `EITrans`, first install the following dependent packages:

```
install.packages("https://github.com/Weiming-Hu/AnalogsEnsemble/raw/master/RAnalogs/releases/RAnEn_latest.tar.gz", repos = NULL)
devtools::install_github('Weiming-Hu/RAnEnExtra')

install.packages(c('progress', 'abind', 'doSNOW', 'foreach'))
```

The following packages are optional. They are mainly for high performance computing.

```
install.packages('Rmpi')
```

Then, you can install `EITrans` by running

```
devtools::install_github('Weiming-Hu/EITrans')
```

```
# "`-''-/").___..--''"`-._
#  (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
#  (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#    _ ..`--'_..-_/  /--'_.' ,'
#  (il),-''  (li),'  ((!.-'
# 
# Author: 
#     Weiming Hu <weiming@psu.edu>
#         
# Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
# Department of Geography and Institute for CyberScience
# The Pennsylvania State University
```
