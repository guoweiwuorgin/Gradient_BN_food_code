# Gradient_BN_food_code

This code is used to run FC gradient in BN and HC group, and use R run statistic, the gradient calculation code mainly come from Xia Mingrui (2023, Moculer Psychiatry)

# Research Code Documentation: Motion Filtering, Gradient Analysis, and Statistical Processing

## 1. Introduction

This document outlines the code structure and functionality used in our research project. Our analysis pipeline consists of two main steps: 
1. Subject screening and motion filtering
2. Gradient computation, statistical analysis, and visualization

## 2. Code Structure

Our code repository contains the following main files:

### Step 1:
1. `gradient_motion_filter.R`: Motion filtering script
2. `a_analysis_pipeline.m`: Main script for individual gradient computation
3. `x_GRF.m`: Auxiliary function
4. `x_compute_diffusion_map.m`: Diffusion map computation function
5. `x_gen_matrix_voxel.m`: Voxel matrix generation function
6. `x_reslice.m`: Image reslicing function

### Step 2:
1. `Gradient_statistical_andplotting.R`: R script for statistical analysis and plotting
2. `StatForGradient_range_and_voxel.m`: MATLAB script for gradient range and voxel-wise analysis
3. `gradient_plot.Rmd`: R Markdown file for generating gradient plots and reports

## 3. Analysis Pipeline

### 3.1 Step 1: Motion Filtering and Gradient Computation

[Content remains the same as in the previous version]

### 3.2 Step 2: Statistical Analysis and Visualization

This step focuses on statistical processing of the gradients, range analysis, voxel-wise computations, and generating visualizations.

#### 3.2.1 Gradient Statistics and Plotting

**File:** `Gradient_statistical_andplotting.R`

**Usage:**
1. Open R or RStudio.
2. Set the working directory to the location of your R script.
3. Run the script:
   ```R
   source("Gradient_statistical_andplotting.R")
   ```

**Key Functions:**
- Statistical tests on gradient data
- Generation of various plots (e.g., box plots, scatter plots)
- Group-level analyses

**Output:**
- Statistical test results
- Visualization plots
- Summary statistics

#### 3.2.2 Gradient Range and Voxel-wise Analysis

**File:** `StatForGradient_range_and_voxel.m`

**Usage:**
1. Open MATLAB.
2. Navigate to the directory containing the script.
3. Run the script:
   ```matlab
   run('StatForGradient_range_and_voxel.m')
   ```

**Key Functions:**
- Computation of gradient ranges across subjects or conditions
- Voxel-wise statistical analyses
- Generation of statistical maps

**Output:**
- Gradient range statistics
- Voxel-wise statistical maps
- Data files for further analysis or visualization

#### 3.2.3 Gradient Plotting and Reporting

**File:** `gradient_plot.Rmd`

**Usage:**
1. Open RStudio.
2. Open the `gradient_plot.Rmd` file.
3. Knit the document to generate the report.

**Key Functions:**
- Creation of publication-quality gradient plots
- Generation of comprehensive reports
- Integration of statistical results with visualizations

**Output:**
- HTML, PDF, or Word document containing gradient plots and analysis results
- Interactive visualizations (if applicable)

## 4. Data Requirements

- Processed gradient data from Step 1
- Subject metadata (e.g., group assignments, demographic information)
- ROI or atlas information for region-specific analyses

## 5. Dependencies

- R (version 4.0.0 or later) with packages:
  - tidyverse
  - ggplot2
  - rmarkdown
  - [Add any additional R packages used in your scripts]
- MATLAB (R2019b or later) with toolboxes:
  - Statistics and Machine Learning Toolbox
  - [Add any additional MATLAB toolboxes required]
- RStudio (for R Markdown processing)

## 6. Notes for Researchers

- Ensure that the output from Step 1 is properly organized and accessible for Step 2 analyses.
- Modify file paths in the scripts to match your directory structure.
- The `StatForGradient_range_and_voxel.m` script may require adjustment based on your specific gradient computation method.
- When using `gradient_plot.Rmd`, make sure to update any parameters or file paths at the beginning of the document.

## 7. Citation

If you use this code in your research, please cite our paper:

[Insert citation details here]

## 8. Contact

For questions or issues regarding the code, please contact:
wuguowei@psych.ac.cn
We welcome contributions and suggestions for improving this pipeline. Please submit issues or pull requests to our project repository.
