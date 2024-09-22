# Gradient_BN_food_code

This code is used to run FC gradient in BN and HC group, and use R run statistic, the gradient calculation code mainly come from Xia Mingrui (2023, Moculer Psychiatry)

# Research Code Documentation: Motion Filtering and Individual Gradient Analysis

## 1. Introduction

This document aims to elucidate the code structure and functionality used in our research project. Our analysis pipeline primarily consists of two key steps: subject screening (based on motion filtering) and individual gradient computation with metric analysis.

## 2. Code Structure

Our code repository contains the following main files:

1. `gradient_motion_filter.R`: Motion filtering script
2. `a_analysis_pipeline.m`: Main script for individual gradient computation and metric analysis
3. `x_GRF.m`: Auxiliary function
4. `x_compute_diffusion_map.m`: Diffusion map computation function
5. `x_gen_matrix_voxel.m`: Voxel matrix generation function
6. `x_reslice.m`: Image reslicing function

## 3. Analysis Pipeline

### 3.1 Step 1: Motion Filtering (Subject Screening)

Utilize the `gradient_motion_filter.R` script for subject screening. This step evaluates subject head motion during scanning based on the Framewise Displacement (FD) metric.

#### Usage:
1. Ensure the R environment is properly configured with all necessary dependencies installed.
2. In the R environment, run the following command:
   ```R
   source("gradient_motion_filter.R")
   ```
3. Input required parameters as prompted by the script, such as FD threshold and data path.

#### Output:
- A text file containing a list of subjects who passed the motion screening.
- Statistical summary of motion parameters.

### 3.2 Step 2: Individual Gradient Computation and Metric Analysis

Use `a_analysis_pipeline.m` as the main script, which calls other auxiliary functions to perform individual gradient computation and analysis of related metrics.

#### Usage:
1. Open `a_analysis_pipeline.m` in the MATLAB environment.
2. Set necessary parameters such as data path, output directory, etc.
3. Run the script:
   ```matlab
   run('a_analysis_pipeline.m')
   ```

#### Key Functions:
- `x_GRF.m`: Computes the Gradient Response Function
- `x_compute_diffusion_map.m`: Constructs and computes the diffusion map
- `x_gen_matrix_voxel.m`: Generates voxel-level analysis matrices
- `x_reslice.m`: Performs image reslicing to ensure spatial consistency

#### Output:
- Individual gradient maps
- Derived metrics (e.g., gradient magnitude, gradient orientation)
- Statistical results and visualizations

## 4. Data Requirements

- Raw fMRI data (preprocessed according to standard pipelines)
- Motion parameters file for each subject
- Anatomical images (for registration purposes)

## 5. Dependencies

- R (version 3.5.0 or later) with packages:
  - tidyverse
  - stringr
- MATLAB (R2019b or later) with toolboxes:
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox

## 6. Notes for Researchers

- Ensure all data paths are correctly specified before running the scripts.
- The FD threshold in `gradient_motion_filter.R` may need adjustment based on your specific study requirements.
- Intermediate results are saved at various stages of the pipeline to facilitate debugging and allow for checkpoint resumption.
