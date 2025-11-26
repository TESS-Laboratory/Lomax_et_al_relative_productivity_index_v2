# The Relative Productivity Index V2

This repository contains R and Google Earth Engine code for the analysis in Lomax, Guy A., Thomas W.R. Powell, Timothy M. Lenton, and Andrew M. Cunliffe, (in prep) ‘The relative productivity index for assessing local impacts  on rangeland condition: refinement and evaluation’.

## Getting started

To get a local copy of the project up and running, follow the steps below.

1. **Clone the repository from GitHub**

*From RStudio Projects*
- In RStudio, click File > New Project
- Click "Version Control > Git"
- In the Repository URL box, paste `https://github.com/TESS-Laboratory/Lomax_et_al_relative_productivity_index_v2`
- Edit the folder location if desired.
- Click "Create Project".

*From the terminal*

- In an open terminal, navigate to the folder in which you would like to clone the project.
- Run the following command:
`git clone https://github.com/TESS-Laboratory/Lomax_et_al_relative_productivity_index_v2`

2. **Set up project environment using renv**
- Install the renv package using `install.packages("renv")`
- Active renv using `renv::activate()`
- Load the project environment from the lockfile with `renv::restore()`

3. **Download and extract input datasets**
- Download "data.zip" from XXXX.
- Extract data.zip to project directory.

4. **Run the scripts in numeric order**
- Scripts 1a and 1b calculate two spatial covariates - topographic wetness index and distance to river.
- Script 2 tunes a quantile regression forest model on a sample of the input covariates.
- Script 3a predicts the tuned model to calculate RPI across the study area; Script 3b calculates RESTREND residuals for the same area.
- Script 4 calculates performance metrics for RPI and RESTREND methods.
- Scripts 5a, 5b and 5c conduct analysis for the three evaluation case studies.
