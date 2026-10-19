# DLQ: Detection, localization, and quantification of single-source methane emissions

This repository contains the R code and data for the detection, localization, and quantification (DLQ) framework described in:

> Daniels, W. S., Jia, M., and Hammerling, D. M. (2024). Detection, localization, and quantification of single-source methane emissions on oil and gas production sites using point-in-space continuous monitoring systems. *Elementa: Science of the Anthropocene*, 12(1), 00110. https://doi.org/10.1525/elementa.2023.00110

## Contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Data](#data)
- [Requirements](#requirements)
- [Usage](#usage)
- [Applying DLQ to your own data](#applying-dlq-to-your-own-data)
- [Computational notes](#computational-notes)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Overview

DLQ estimates when a methane emission starts and ends (detection), which potential source it comes from (localization), and how large it is (quantification), using methane concentration and wind data from a network of point-in-space continuous monitoring system (CMS) sensors. The framework is modular and is not tied to a particular sensor technology: it needs methane concentrations from the sensors and wind speed and direction from at least one of them. It requires no training data and no prior knowledge of emission start or end times, so it can be run continuously on incoming data, and a localization and quantification estimate can be produced from as little as 15 minutes of data.

The framework has four steps, all described in Section 2.3 of the paper. First, the background is removed from each sensor's concentration time series using a gradient-based spike detection algorithm (Section 2.2), so that only concentration enhancements remain. Second, the enhancements are combined across sensors and grouped into emission events with estimated start and end times. Third, a Gaussian puff atmospheric transport model simulates the concentrations that each potential source would produce at each sensor, given the measured wind. Fourth, for each event, the source whose simulated concentrations best match the observations (highest correlation) is taken as the localization estimate, and the emission rate is estimated by scaling the simulated concentrations to the observations, with a bootstrapped 90% confidence interval.

The code here applies DLQ to one month (April 17 to May 16, 2022) of non-blinded single-source controlled releases at the Methane Emissions Technology Evaluation Center (METEC) in Fort Collins, Colorado, collected as part of the Advancing Development of Emissions Detection (ADED) program. The releases came from five potential sources (two wellheads, two separators, and a set of tanks) and were measured by eight CMS sensors, three of which also measured wind. In that evaluation, all controlled releases were detected, 82% were localized correctly, and 5.5% of the estimated events were false positives. See the paper for the quantification results.

The framework is designed for single-source emissions on production sites of roughly 150 m in diameter with spatially distinct equipment groups, surrounded by enough sensors (approximately 4 to 8) that a plume is likely to reach at least one of them regardless of wind direction. Multisource and off-site emissions are outside its scope; the seven multisource releases in the METEC month are excluded from the evaluation. `MAIN_1_simulate.R` and `MAIN_2_DLQ.R` work for any number of sources and sensors, but `MAIN_3_analyze_results.R` is specific to the METEC experiment.

**Note:** this README was generated with Claude Sonnet 5.5.

## Repository structure

```
DLQ/
├── code/
│   ├── MAIN_1_simulate.R                   # runs the Gaussian puff model at the sensor locations
│   ├── MAIN_2_DLQ.R                        # background removal, event detection, localization, quantification
│   ├── MAIN_3_analyze_results.R            # evaluates results against the controlled releases, makes figures
│   ├── HELPER_gpuff_function.R             # Gaussian puff model: stability classes, dispersion parameters, concentrations
│   ├── HELPER_distance_conversions.R       # latitude/longitude to UTM conversion
│   └── HELPER_spike_detection_algorithm.R  # find.spikes(), the spike detection algorithm
├── input_data/
│   ├── ADED_data_clean_SAMPLE.csv          # one week of CMS observations (sample of the full data)
│   ├── ADED_data_clean_README.txt          # column descriptions for the CMS data
│   ├── sensor_locations.csv
│   ├── source_locations.csv
│   ├── leak_data.csv                       # METEC controlled release records (ground truth)
│   ├── leak_data_README.xlsx               # column descriptions for leak_data.csv
│   ├── simulation_config.txt               # settings for MAIN_1_simulate.R
│   ├── simulation_config_README.txt
│   ├── DLQ_config.txt                      # settings for MAIN_2_DLQ.R
│   └── DLQ_config_README.txt
├── output_data/
│   ├── simulation_output_SAMPLE.RData      # MAIN_1 output for the sample week
│   ├── DLQ_output_SAMPLE.RData             # MAIN_2 output for the sample week
│   └── DLQ_output.RData                    # MAIN_2 output for the full month (input to MAIN_3)
├── figures/                                # populated by MAIN_3_analyze_results.R
├── LICENSE
└── README.md
```

The scripts are run in order, with each one reading the output of the previous one:

```
CMS observations --MAIN_1--> simulation_output --MAIN_2--> DLQ_output --MAIN_3--> summary statistics and figures
```

`HELPER_spike_detection_algorithm.R` is also loaded directly from this repository by the [MDLQ](https://github.com/wsdaniels/MDLQ) code, so please do not move or rename it.

## Data

### Input files

| File | Contents |
|---|---|
| `input_data/ADED_data_clean_SAMPLE.csv` | Minute-frequency observations from the eight CMS sensors for one week (2022-05-09 00:01 to 2022-05-15 23:59 UTC). Columns used by the code: `time` (UTC), `name` (sensor name), `methane` (ppm), `wind.speed` (m/s), and `wind.direction` (degrees; the direction the wind is coming from, clockwise from north). Wind is `NA` for sensors that do not measure it (in this file, only `E`, `NW`, and `SW` do). The remaining columns are ignored. |
| `input_data/sensor_locations.csv` | `lon`, `lat`, `height` (m above ground), and `name` of the eight CMS sensors. The names must match the `name` column of the observation file. |
| `input_data/source_locations.csv` | `lon`, `lat`, `height` (m above ground, the release height), and `name` of the five potential sources: `Tank` (4.5 m), `Separator.East`, `Wellhead.East`, `Wellhead.West`, and `Separator.West` (2 m). |
| `input_data/leak_data.csv` | METEC records of the controlled releases from February 1 to May 16, 2022, including equipment group, start and end time (UTC), metered release rate (`tc_C1MassFlow`, g/h), and its uncertainty. Used only by `MAIN_3_analyze_results.R`. Column definitions are in `leak_data_README.xlsx`. |
| `input_data/simulation_config.txt`, `input_data/DLQ_config.txt` | Settings for `MAIN_1_simulate.R` and `MAIN_2_DLQ.R` (see [Usage](#usage)). |

Due to the proprietary nature of the CMS sensor technology used in this study, only a sample of the concentration and wind measurements can be shared publicly, which is why `ADED_data_clean_SAMPLE.csv` covers one week rather than the whole month. The full data can be made available for research purposes upon request (see [Contact](#contact)). Without them, you can run `MAIN_1_simulate.R` and `MAIN_2_DLQ.R` on the sample week, and `MAIN_3_analyze_results.R` on the pre-generated full-month results in `output_data/DLQ_output.RData`, which reproduces the results in the paper.

### Output files

All files in `output_data/` have the extension `.RData`, but they were written with `saveRDS()` and must be read with `readRDS()`, not `load()`.

**`simulation_output_SAMPLE.RData`** (written by `MAIN_1_simulate.R`) is a named list:

| Element | Description |
|---|---|
| `times` | Minute time stamps (POSIXct, in the time zone set by `tz`). Each time stamp marks the end of its minute. |
| `obs` | Data frame of observed methane concentrations (ppm), one row per minute and one column per sensor. Columns are sorted alphabetically by sensor name. |
| `WD`, `WS` | Wind direction and speed for each minute, taken as the median across the sensors that measure wind. `WD` is in radians, in the standard mathematical convention (counterclockwise from east) and points in the direction the wind is blowing toward. `WS` is in m/s. |
| one element per source (`Tank`, `Separator.East`, ...) | Data frame with the same dimensions as `obs`: the simulated concentrations (ppm) at each sensor for an emission of 1 g/s from that source. |

**`DLQ_output_SAMPLE.RData`** (written by `MAIN_2_DLQ.R`) is a named list:

| Element | Description |
|---|---|
| `event.mask` | Data frame with a `time` column and an `events` column. `events` is `NA` outside of events and an integer identifying the event otherwise. Only meaningful when `do.event.detection=T`. |
| `max.obs` | Maximum background-removed concentration across sensors (ppm) for each minute. |
| `localization.estimates` | Estimated source name for each event (or interval). |
| `rate.estimates`, `error.lower`, `error.upper` | Estimated emission rate (kg/h) for each event, with the 5th and 95th percentiles of the bootstrapped rates (a 90% interval). |
| `source.names`, `WD`, `WS` | Source names and the wind data, trimmed to match `event.mask`. |
| `out` | List with one entry per event. Each entry has `q.hat`, `q.hat.lower`, and `q.hat.upper` (one value per source, in the order of `source.names`) and `rates` (the bootstrapped rate estimates, kg/h). |
| one element per source | The simulated concentrations from the input, trimmed to match `event.mask`. |

Estimates are produced only for the most likely source; the rate for every other source is set to 0 (the single-source assumption). Where an estimate cannot be produced, the value is a text label instead of a number: `"no info"` means that the simulated concentrations from that source are too small during the event to be informative (typically because no sensor is downwind of it), and `"no overlap"` means that the event has concentration enhancements but there are too few time steps in which the observations and simulations are both in a spike to estimate a rate. Because of these labels, the `q.hat`, `q.hat.lower`, and `q.hat.upper` entries are stored as character vectors and should be converted with `as.numeric()`.

**`DLQ_output.RData`** holds the results for the full month and contains only the first nine elements above (`event.mask` through `WS`), without `out` or the simulated concentrations. In this file, events without a rate estimate are stored as `NA`.

## Requirements

The code is written in R and was tested with R 4.3.3. It uses the packages `lubridate`, `zoo`, `foreach`, `doParallel`, `scales`, and `rstudioapi`:

```r
install.packages(c("lubridate", "zoo", "foreach", "doParallel", "scales", "rstudioapi"))
```

`rstudioapi` is loaded by every main script and so must be installed even if you do not use RStudio. `MAIN_1_simulate.R` uses `zoo`, `lubridate`, `foreach`, and `doParallel`; `MAIN_2_DLQ.R` uses `lubridate` and `zoo`; and `MAIN_3_analyze_results.R` uses `lubridate` and `scales`. The scripts have no other dependencies and do not need an internet connection.

## Usage

All paths in the scripts and config files are relative to the `code` directory, so run the scripts with `code` as the working directory. Inside RStudio the scripts set this automatically; otherwise use `setwd("code")` or run from a terminal:

```bash
cd code
Rscript MAIN_1_simulate.R
Rscript MAIN_2_DLQ.R
Rscript MAIN_3_analyze_results.R
```

The pre-generated files in `output_data/` let you start at any step. Re-running a step overwrites the pre-generated files that it writes (MAIN_1 and MAIN_2 write to the paths set in their config files, and MAIN_3 overwrites the images in `figures/`), so change the output paths first if you want to keep them.

Both config files contain one `setting=value` pair per line, with no spaces around the equals sign and no quotation marks.

### 1. Simulate concentrations: `MAIN_1_simulate.R`

This script reads the CMS observations and the sensor and source locations, and simulates the concentration that a constant emission from each potential source would produce at each sensor, using the Gaussian puff model. Edit `input_data/simulation_config.txt`:

| Setting | Default value | Meaning |
|---|---|---|
| `num.cores.to.use` | 7 | Number of cores to parallelize over. Should not exceed the number of cores on your machine. |
| `dt` | 1 | Simulation time step in seconds: one puff is released every `dt` seconds. `60/dt` must be an even number (for example, `dt` = 1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, or 60). |
| `cutoff.t` | 20 | Time in minutes after which a puff is dropped. A good rule of thumb is the diameter of the study area divided by the 5th percentile of the wind speed. |
| `ignore.dist` | 80 | Distance in meters. A puff farther than this from every sensor is skipped, which speeds up the simulation. |
| `chunk.size` | 54 | Length in minutes of the time chunks that the simulation is split into for parallel processing. Puffs that cross a chunk boundary are carried into the next chunk. Chunks are the unit of parallelization, so a simulation shorter than `chunk.size` runs on one core. |
| `run.mode` | constant | `constant` fixes each puff's dispersion parameters, direction, and speed at the conditions when it is released (used in the paper, as it can be derived from the advection-diffusion equation). `dynamic` updates them as the puff moves. |
| `emission.rate` | 0.001 | Emission rate used in the simulation, in kg/s. Leave this at 0.001 (1 g/s); see the note below. |
| `start.time`, `end.time` | 2022-05-9T00:00:00, 2022-05-16T00:00:00 | Start and end of the simulation, in the format `YYYY-MM-DDTHH:MM:00`, interpreted in the time zone `tz`. They must fall on a whole minute. |
| `tz` | America/Denver | Local time zone, in the "area/location" convention of the [tz database](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). |
| `raw.sensor.observations.path` | `../input_data/ADED_data_clean_SAMPLE.csv` | CMS observations. |
| `source.locations.path` | `../input_data/source_locations.csv` | Potential sources. A separate simulation is run for each. |
| `sensor.locations.path` | `../input_data/sensor_locations.csv` | CMS sensors. |
| `output.file.path` | `../output_data/simulation_output_SAMPLE.RData` | Where the simulation output is saved. |
| `helper.distance.conversions.path`, `helper.gpuff.function.path` | `../code/HELPER_distance_conversions.R`, `../code/HELPER_gpuff_function.R` | Helper files. |

Two points about the settings are easy to miss. First, the `time` column of the observation file is in UTC, whereas `start.time` and `end.time` are interpreted in `tz`. With the default settings, the simulation therefore starts at 06:00 UTC on May 9 (midnight local time), not at the beginning of the file, which is why the sample simulation has 9,720 minutes. Second, `MAIN_2_DLQ.R` converts its fitted rates to kg/h assuming that the simulation used 1 g/s. If you change `emission.rate`, multiply the rate estimates from `MAIN_2_DLQ.R` (`rate.estimates`, `error.lower`, `error.upper`, and the `rates` in `out`) by `emission.rate / 0.001`.

`MAIN_1_simulate.R` is the slowest step in the workflow (see [Computational notes](#computational-notes)). For a quick test, shorten `start.time` and `end.time` to a few minutes and point `output.file.path` somewhere else.

### 2. Detect, localize, and quantify: `MAIN_2_DLQ.R`

This script reads the output of `MAIN_1_simulate.R`, removes the background from the observations, detects emission events, and then localizes and quantifies each one. Edit `input_data/DLQ_config.txt`:

| Setting | Default value | Meaning |
|---|---|---|
| `gap.time` | 30 | Events separated by `gap.time` minutes or less are combined into one event. |
| `length.threshold` | 15 | Events shorter than `length.threshold` minutes are discarded. |
| `do.event.detection` | T | `T` estimates emission start and end times and then localizes and quantifies each resulting event (the method in the paper). `F` skips event detection and instead localizes and quantifies each non-overlapping interval of `gap.time` minutes. |
| `first.sim.ind` | 5 | Index (starting at 1) of the first source in the list saved by `MAIN_1_simulate.R`, i.e., the number of elements before the simulations (`times`, `obs`, `WD`, `WS`) plus one. |
| `n.samples` | 1000 | Number of bootstrap samples used for the rate estimate and its interval. |
| `forward.model.path` | `../output_data/simulation_output_SAMPLE.RData` | Output of `MAIN_1_simulate.R`. |
| `output.file.path` | `../output_data/DLQ_output_SAMPLE.RData` | Where the results are saved. |
| `helper.spike.detection.alg.path` | `../code/HELPER_spike_detection_algorithm.R` | Helper file with the spike detection algorithm. |

On the sample week this step takes a few seconds and finds 25 events. The bootstrap in `MAIN_2_DLQ.R` does not set a random seed, so the rate estimates and intervals change slightly from run to run, including relative to the provided `DLQ_output_SAMPLE.RData`. In a repeat run on the sample week, the event times and localization estimates were identical to the provided ones.

### 3. Evaluate the results: `MAIN_3_analyze_results.R`

This script compares the DLQ results to the METEC controlled releases, prints the summary statistics, and saves figures. The three paths at the top of the script (`save.dir`, `data.path`, and `leak.path`) are set to `../figures/`, `../output_data/DLQ_output.RData`, and `../input_data/leak_data.csv`, so by default it analyzes the full-month results and reproduces the paper. It takes a few seconds to run. Controlled releases that overlap in time (multisource releases) are removed from the evaluation, along with any estimated events that overlap with them.

To look at the sample week instead, point `data.path` at `../output_data/DLQ_output_SAMPLE.RData` (and `save.dir` at a different folder). The script runs, but the sample week is too short to reproduce the numbers in the paper.

The figures correspond to the paper as follows:

| File in `figures/` | Contents |
|---|---|
| `result_time_series.png` | Figure 2: controlled releases (solid) and estimated events (transparent), with estimated rates and intervals, over the month |
| `detection_localization_results.png` | Figure 3: event-level detection, localization, and rate availability, by true source |
| `parity.png` | Figure 4: parity plot of estimated and true emission rates |
| `over_under_est.png` | Figure 5: percent difference between estimated and true rates, for releases of 1 kg/h or less and above 1 kg/h |
| `timing_error.png` | Histograms of the start and end time errors (Section 3) |
| `cumulative.png` | Cumulative true and estimated emissions as a function of true emission rate (Section 3; Section S10 of the Supporting Information) |
| `leaks_1.png` to `leaks_6.png`, `events_1.png` to `events_6.png` | Diagnostic plots of the true releases and of the detected events, in 5-day windows |

The summary statistics are printed to the console, without labels, in the order in which they are computed: the event-level false positive rate, quantiles of the start and end time errors, the minute-by-minute confusion matrix, the fractions of estimates within a factor of 2 and 3 of the true rate, the percent differences for small and large emissions, and the errors in cumulative emissions.

## Applying DLQ to your own data

`MAIN_1_simulate.R` and `MAIN_2_DLQ.R` can be used on other sites and sensor networks. To do so:

1. Create an observation file with the columns `time` (UTC, rounded to the minute), `name`, `methane` (ppm), `wind.speed` (m/s), and `wind.direction` (degrees the wind is coming from, clockwise from north), with one row per sensor per minute. At least one sensor must report wind. Missing values (`NA`) are allowed. The simulation requires minute-frequency data.
2. Create `sensor_locations.csv` and `source_locations.csv` with the columns `lon`, `lat`, `height`, and `name`, using the same sensor names as in the observation file. The potential sources should include every piece of equipment that could emit, because DLQ can only assign an event to a source in this list.
3. Edit the paths, `start.time`, `end.time`, and `tz` in `simulation_config.txt`, and the paths in `DLQ_config.txt`. Set `first.sim.ind` to 5 unless you change the format of the simulation output.
4. Run `MAIN_1_simulate.R` and then `MAIN_2_DLQ.R`.

A few details of the forward model are specific to the way it was set up for the paper. The Pasquill stability class used for the dispersion parameters is chosen from the wind speed and the time of day, with daytime hard-coded as 07:00 to 18:59 local time, which was set up for summer conditions (`is.day()` in `HELPER_gpuff_function.R`). The time zone `tz` therefore matters. `MAIN_3_analyze_results.R` maps the equipment group IDs in `leak_data.csv` to source names according to the order in which they first appear in that file, so it works only with the METEC data as provided.

## Computational notes

- **`MAIN_1_simulate.R` is by far the slowest step.** The cost grows with the length of the simulation, the number of sources, `cutoff.t`, and `1/dt`, and the simulation for each source is parallelized over time chunks of `chunk.size` minutes. Simulating the full sample week for five sources is expected to take hours rather than minutes, even in parallel; a ten-minute simulation takes on the order of a minute. If you only want to reproduce the paper's results, use the pre-generated files in `output_data/`.
- `MAIN_2_DLQ.R` and `MAIN_3_analyze_results.R` each take seconds on the sample and full-month data.
- `MAIN_1_simulate.R` is deterministic. `MAIN_2_DLQ.R` is not (see above), so re-running it gives rate estimates that differ slightly from the provided files. `MAIN_3_analyze_results.R` is deterministic given its input.

## Troubleshooting

- **`argument is of length zero` in `latlon.to.utm()`, or a column named `X...lon`:** the two location files begin with a UTF-8 byte-order mark, which R discards only when it runs in a UTF-8 locale. Start R in a UTF-8 locale (for example, `LC_ALL=en_US.UTF-8` or `LC_ALL=C.UTF-8`) or re-save the CSV files without the byte-order mark.
- **`cannot open the connection` or `cannot open file`:** check that the working directory is `code` and that the paths in the config files are correct. Output folders must already exist.
- **`bad restore file magic number` when reading an output file:** use `readRDS()` instead of `load()`.
- **`there is no package called 'rstudioapi'`:** install `rstudioapi` even if you do not use RStudio.
- **Times in the output are shifted by several hours:** `start.time` and `end.time` are in `tz`, whereas the `time` column of the observations is in UTC.

## Citation

If you use this code in your research, please cite the accompanying paper:

Daniels, W. S., Jia, M., and Hammerling, D. M. (2024). Detection, localization, and quantification of single-source methane emissions on oil and gas production sites using point-in-space continuous monitoring systems. *Elementa: Science of the Anthropocene*, 12(1), 00110. https://doi.org/10.1525/elementa.2023.00110

```bibtex
@article{daniels_detection_2024,
	title = {Detection, localization, and quantification of single-source methane emissions on oil and gas production sites using point-in-space continuous monitoring systems},
	volume = {12},
	issn = {2325-1026},
	doi = {10.1525/elementa.2023.00110},
	number = {1},
	urldate = {2024-07-31},
	journal = {Elementa: Science of the Anthropocene},
	author = {Daniels, William S. and Jia, Meng and Hammerling, Dorit M.},
	month = mar,
	year = {2024},
	pages = {00110},
}
```

The Gaussian puff atmospheric transport model used for the forward simulations is described in Jia, M., Fish, R., Daniels, W. S., Sprinkle, B., and Hammerling, D. (2025). A fast and lightweight implementation of the Gaussian puff model for near-field atmospheric transport of trace gasses. *Scientific Reports*, 15, 18710. https://www.nature.com/articles/s41598-025-99491-x. That paper presents a faster implementation of the model, available in [FastGaussianPuff](https://github.com/Hammerling-Research-Group/FastGaussianPuff); the code in this repository (`HELPER_gpuff_function.R`) is the R implementation used in the DLQ paper. The [MDLQ](https://github.com/wsdaniels/MDLQ) repository contains a related Bayesian model for multisource emissions, which reuses the background removal from this repository.

## License

This project is released under the MIT License; see [LICENSE](LICENSE).

## Contact

William Daniels (wdanie16@jh.edu). Please also contact the author to request access to the full CMS concentration and wind data from the METEC experiment for research purposes.
