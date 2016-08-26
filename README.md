# Binding_Affinity_Prediction
This requires R 3.2.3 at least.

It also depends on the following packages:

- glmnet >=1.8-12
- caret >=6.0-70
- mgcv >=1.8-12
- KRLS >=0.3-7
- hydroGOF >=0.3-8
- optparse >=1.3.2
- ggplot2 >=2.1.0

To ease dependency management, this project uses packrat (see https://rstudio.github.io/packrat/). For this, you need the following project files:

- .Rproject
- packrat/init.R
- packrat/packrat.lock
- packrat/packrat.opts

To install the dependencies, set your working directory to the project top directory, and run (you need packrat >=0.4.7 installed):

```R
> library(packrat)
> packrat::init()
```

At the question of removing the dependencies, reply "N" (since you want them to be installed, instead).

To run using the packrat installed dependencies, you cannot use Rscript --vanilla, since it does not reads the .Rproject file which configures the system to use packrat. Instead, run

```
$ Rscript --no-save --no-restore --no-site-file --no-environ myscript.R
```
