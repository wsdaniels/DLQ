# Detection, localization, and quantification (DLQ) using continuous monitoring systems (CMS)

Repository contains code used to estimate methane emission start and end time (detection), source location (localization), and emission rate (quantification) using concentration observations from a network of point-in-space continuous monitoring systems. The code is separated into two main scripts located in the "code" directory: 1) MAIN_1_simulate runs the Gaussian puff atmospheric dispersion model, and 2) MAIN_2_DLQ uses output from the Gaussian puff model to perform DLQ.

Inputs to the MAIN_1 and MAIN_2 files are controlled using two configuration files found in the "input_data" directory. The "simulation_config.txt" file controls input for the MAIN_1 script and the "DLQ_config.txt" file controls input for the MAIN_2 script. A README file is provided for each config file in the "input_data" directory.

The MAIN_3 script generates all results and figures for the accompanying manuscript (Daniels et al. 2023): https://doi.org/10.26434/chemrxiv-2022-xxkk8-v3

Note that the "input_data" directory also contains the raw concentration data, sensor locations, and source locations from the ADED experiment discussed in Daniels et al. 2023. The "output_data" directory is where output from the MAIN_1 and MAIN_2 scripts is saved. Output files have been pre-generated and are saved in these directories.

