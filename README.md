GCRMN model development
===========================

## Building the environment

Most of the dependencies can be inferred by examination of the
`Dockerfile`. In fact, the safest way of ensuring that the codes will
run is the build a docker image from the `Dockerfile` and run within a
container.

## Running in the container

To run these models inside the docker container (after building the
container), use the following steps:

1. `make R_container`: this will run the the following analyses in R:
  - `synthetic_data()`: to generate the synthetic data
  - `site_replacement()`: R based models for exploring site
    replacement scenarios
  - `missing_years()`: R based models for exploring missing years
    scenarios
  - `incomplete_spatial()`: R based models for exploring incomplete
    spatial scenarios
2. `make python_container`: this will run the the following analyses
   in R:
  - `site_replacement()`: python based models for exploring site
    replacement scenarios
  - `missing_years()`: python based models for exploring missing years
    scenarios
  - `incomplete_spatial()`: python based models for exploring
    incomplete spatial scenarios
3. `make R_container`: a second time (after `make python_container`) in order to 
   combine the previous R and python models together for the purpose
   of comparing the models
4. `make docs_container`: this will render each of the quarto documents

## Running the codes

Nevertheless, if you are running on bare metal, then ensure that both
`R` and `python` are installed and that their respective packages
indicated in the Dockerfile are installed and available.

Ideally, it should be possible to run in the following order:

1. `make run_R`: this will run all the R based analyses using the
   `targets` package to ensure all steps are performed in the correct
   order
2. `make run_python`: this will run all the python based analyes using
   the `ploomber` library to ensure all steps are performed in the
   correct order
3. `make run_R`: a second time (after `make run_python`) in order to 
   combine the previous R and python models together for the purpose
   of comparing the models
4. `make render_docs`: this will render each of the quarto documents

I say ideally, because in my haste and to quickly knit together the
python analyses and R analyses so that they can be compared, I
included R code to read in the results of the python analyses. Of
course the first time that this all gets run, those python analyses
wont occur and so it will error. Unfortunately, it is not as simple as
just running all the python analyses first (e.g. swapping steps 1 and
2 above) because the python analyses rely on the simulated landscapes
that are created early on in the R analyses.

I should re-write this to be a bit more granular so that it runs the
first part of the R analyses, then python, then the later R steps -
however I ran out of time!  I will hopefully do this step soon.

## Github actions to publish the documents as github pages

I have written github workflows to automate the whole process on a git
commit. Unfortunately, the workflow hits the hard walltime of 6hrs and
therefore fails. I am currently exploring options for engaging a
runner to allow full continuous integration.


