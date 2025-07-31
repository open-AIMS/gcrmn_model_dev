## Get R version 4.4.3
FROM rocker/r-ver:4.4.3

## Install packages
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libudunits2-dev \
    libssl-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    cmake \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    glpk-utils \
    libglpk-dev \ 
    git \ 
  && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  install.packages(\"pak\"); \
"

RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  pak::pkg_install(c('rmarkdown', 'quarto')); \
  pak::pkg_install(c('tidyverse')); \
  pak::pkg_install(c('sf', 'knitr', 'patchwork')); \
"

RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  pak::pkg_install(c('testthat', 'usethis', 'remotes')); \
"
RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  pak::pkg_install(c('glmmTMB', 'emmeans', 'DHARMa', 'performance', 'see')); \
  pak::pkg_install(c('gbm', 'dbarts')); \
  pak::pkg_install(c('simstudy')); \
"
RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  pak::pkg_install(c('targets', 'tarchetypes')); \
"
RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  pak::pkg_install(c('graph', 'Rgraphviz', 'gridGraphics')); \
"

RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  pak::pkg_install(c('stan-dev/cmdstanr')); \
  pak::pkg_install(c('brms')); \
  pak::pkg_install(c('tidybayes', 'posterior', 'bayesplot', 'HDInterval', 'jmgirard/standist', 'bayestestR', 'rstan')); \
"

# Install CmdStan
RUN R -e "cmdstanr::check_cmdstan_toolchain(fix = TRUE); \
  cmdstanr::install_cmdstan(cores = parallel::detectCores()); \
"

# Install quarto
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    curl \
    libcurl4-openssl-dev \
    gdebi-core \
  && rm -rf /var/lib/apt/lists/*

ARG QUARTO_VERSION="1.7.26"
RUN curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb


## Python

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  python3 \
  python3-pip \
  python3-dev \
  && rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages --ignore-installed packaging
## RUN pip3 install --break-system-packages pandas
## RUN pip3 install --break-system-packages numpy
## RUN pip3 install --break-system-packages arviz
## RUN pip3 install --break-system-packages seaborn
RUN pip3 install --break-system-packages pymc
RUN pip3 install --break-system-packages pymc-bart
RUN pip3 install --break-system-packages preliz


RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2024-04-11/\")); \
  pak::pkg_install(c('reticulate', 'styler')); \
"



## INLA

RUN R -e "install.packages('INLA',repos=c(getOption('repos'),INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)"

RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2025-04-11/\")); \
  pak::pkg_install(c('open-AIMS/synthos')); \
"

## ## A selection of tidyverse packages
## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('dplyr'); \
##   install.packages('lubridate'); \
##   install.packages('ggplot2'); \
##   install.packages('readr'); \
##   install.packages('stringr'); \
##   install.packages('tidyr'); \
##   install.packages('tidyverse'); \
## "  

## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('crayon'); \
##   install.packages('cli'); \
##   install.packages('validate'); \
## "  

## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('remotes'); \
## "

## ## Project specific packages
## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   remotes::install_github('open-AIMS/status'); \
## "

## ## Other packages
## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('markdown'); \
##   install.packages('bookdown'); \
##   install.packages('rmarkdown'); \
##   install.packages('quarto'); \
## "  

## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('testthat'); \
##   install.packages('assertthat'); \
## "  

## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('plotly'); \
## "  

## ## Install extra packages required for quarto 
## RUN apt-get update \
##   && apt-get install -y --no-install-recommends \
##     curl \
##     gdebi-core \
##   && rm -rf /var/lib/apt/lists/*

## ARG QUARTO_VERSION="1.3.450"
## RUN curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb
## RUN gdebi --non-interactive quarto-linux-amd64.deb


## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('sf'); \
## "  

## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('emmeans');   \
##   install.packages('DHARMa');   \
##   install.packages('patchwork');   \
## "  

## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('future');   \
##   install.packages('purrr');   \
##   install.packages('insight');   \
##   install.packages('gridGraphics');   \
## "  

## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages("fmesher"); \
##   remotes::install_version("INLA", version="24.05.10",repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE); \
## "
  
## RUN R -e "options(repos = \
##   list(CRAN = 'https://packagemanager.posit.co/cran/2025-04-11/')); \
##   install.packages('glmmTMB');   \
##   install.packages('brms');   \
##   install.packages('emdstanr');   \
##   install.packages('cmdstanr');   \
## "  

RUN pip3 install --break-system-packages seaborn

RUN apt-get clean

RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2024-04-11/\")); \
  pak::pkg_install(c('gridGraphics')); \
"

RUN R -e "options(repos = \
    list(CRAN = \"https://packagemanager.posit.co/cran/2024-04-11/\")); \
  pak::pkg_install(c('gbm-developers/gbm')); \
"

RUN pip3 install --break-system-packages ploomber

RUN mkdir /home/Project
WORKDIR /home/Project
